#include "FitCATS3D.h"
#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Histo.h"

#include "TFile.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TAxis.h"

#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

struct ModelCacheEntry3D {
    double Rout = 0.0;
    double Rside = 0.0;
    double Rlong = 0.0;
    double Alpha = 2.0;
    unsigned LastUsed = 0;
    std::unique_ptr<TH3D> Model;
};

struct FitProfile3D {
    double BuildDataSec = 0.0;
    double BuildResonanceSourceSec = 0.0;
    double ConfigureKittySec = 0.0;
    double LoadExternalKernelSec = 0.0;
    double LoadLambdaSec = 0.0;
    double MinuitSec = 0.0;
    double FinalChi2Sec = 0.0;
    double SaveOutputSec = 0.0;

    double CacheLookupSec = 0.0;
    double SetAnaSourceSec = 0.0;
    double KillTheCatSec = 0.0;
    double BuildModelSec = 0.0;
    double Chi2GetModelSec = 0.0;
    double Chi2LoopSec = 0.0;

    int CacheLookupCalls = 0;
    int SetAnaSourceCalls = 0;
    int KillTheCatCalls = 0;
    int BuildModelCalls = 0;
    int Chi2Calls = 0;
};

struct FitState3D {
    FitCATS3DConfig cfg;
    TH3D* data = nullptr;
    TH3D* same = nullptr;
    TH3D* mixed = nullptr;
    CATS3D* kitty = nullptr;
    int nFitBins = 0;
    int nFcnCalls = 0;
    int nCatsEvals = 0;
    int nCacheHits = 0;
    unsigned CacheClock = 0;
    std::vector<ModelCacheEntry3D> ModelCache;
    std::unique_ptr<TH3D> ScratchModel;
    std::unique_ptr<DLM_Histo<double>> LambdaFeedDown;
    DLM_CleverReso3D* ResonanceSource = nullptr;
    FitProfile3D Profile;
    double BkgNormInit = 1.0;
    double BkgParInit[10] = {1.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0};
    double DataNormFactor = 1.0;
};

FitState3D* gFitState3D = nullptr;
double gBkgRejectQminMeV = 100.0;
double gBkgQAxisScaleToMeV = 1000.0;
double gBkgSliceCutMeV = 60.0;

using Clock3D = std::chrono::steady_clock;

double SecondsSince(const Clock3D::time_point& t0)
{
    return std::chrono::duration<double>(Clock3D::now() - t0).count();
}

TH3D* GetTH3DClone(TFile& file, const TString& name, const TString& cloneName)
{
    TH3D* h = dynamic_cast<TH3D*>(file.Get(name));
    if(!h){
        throw std::runtime_error(Form("GetTH3DClone: missing TH3D '%s'", name.Data()));
    }
    gROOT->cd();
    TH3D* out = dynamic_cast<TH3D*>(h->Clone(cloneName));
    out->SetDirectory(nullptr);
    if(out->GetSumw2N() == 0) out->Sumw2();
    return out;
}

TH3D* RebinClone3D(TH3D* h, const int rebinFactor, const TString& name)
{
    if(!h){
        throw std::runtime_error("RebinClone3D: input histogram is null.");
    }

    TH3D* out = dynamic_cast<TH3D*>(h->Clone(name));
    out->SetDirectory(nullptr);
    if(out->GetSumw2N() == 0) out->Sumw2();

    if(rebinFactor > 1){
        out->Rebin3D(rebinFactor, rebinFactor, rebinFactor);
    }

    return out;
}

double RelativeK3MeV(const TH3D* h,
                     const int ix,
                     const int iy,
                     const int iz,
                     const FitCATS3DConfig& cfg)
{
    const double qOut = h->GetXaxis()->GetBinCenter(ix) * cfg.qAxisScaleToMeV;
    const double qSide = h->GetYaxis()->GetBinCenter(iy) * cfg.qAxisScaleToMeV;
    const double qLong = h->GetZaxis()->GetBinCenter(iz) * cfg.qAxisScaleToMeV;
    return cfg.relativeMomentumFactor *
           std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);
}

bool IsInsideFitBoxMeV(const TH3D* h,
                       const int ix,
                       const int iy,
                       const int iz,
                       const FitCATS3DConfig& cfg)
{
    const double qOut = h->GetXaxis()->GetBinCenter(ix) * cfg.qAxisScaleToMeV;
    const double qSide = h->GetYaxis()->GetBinCenter(iy) * cfg.qAxisScaleToMeV;
    const double qLong = h->GetZaxis()->GetBinCenter(iz) * cfg.qAxisScaleToMeV;
    return std::fabs(qOut) < cfg.fitQmaxMeV &&
           std::fabs(qSide) < cfg.fitQmaxMeV &&
           std::fabs(qLong) < cfg.fitQmaxMeV;
}

double Background3DValue(const double qOut,
                         const double qSide,
                         const double qLong,
                         const double* par)
{
    return par[0]
         + par[1] * qOut
         + par[2] * qSide
         + par[3] * qLong
         + par[4] * qOut * qOut
         + par[5] * qSide * qSide
         + par[6] * qLong * qLong
         + par[7] * qOut * qSide
         + par[8] * qOut * qLong
         + par[9] * qSide * qLong;
}

double Background3DValue(const TH3D& h,
                         const int ix,
                         const int iy,
                         const int iz,
                         const double* par)
{
    const double qOut = h.GetXaxis()->GetBinCenter(ix);
    const double qSide = h.GetYaxis()->GetBinCenter(iy);
    const double qLong = h.GetZaxis()->GetBinCenter(iz);
    return Background3DValue(qOut, qSide, qLong, par + 4);
}

double Background3DRejectFemtoTF3(double* x, double* par)
{
    const double qOut = std::fabs(x[0]) * gBkgQAxisScaleToMeV;
    const double qSide = std::fabs(x[1]) * gBkgQAxisScaleToMeV;
    const double qLong = std::fabs(x[2]) * gBkgQAxisScaleToMeV;

    const bool inPure3DCorner =
        qOut >= gBkgRejectQminMeV &&
        qSide >= gBkgRejectQminMeV &&
        qLong >= gBkgRejectQminMeV;
    const bool inOutSlice =
        qOut >= gBkgRejectQminMeV &&
        qSide <= gBkgSliceCutMeV &&
        qLong <= gBkgSliceCutMeV;
    const bool inSideSlice =
        qSide >= gBkgRejectQminMeV &&
        qOut <= gBkgSliceCutMeV &&
        qLong <= gBkgSliceCutMeV;
    const bool inLongSlice =
        qLong >= gBkgRejectQminMeV &&
        qOut <= gBkgSliceCutMeV &&
        qSide <= gBkgSliceCutMeV;

    if(!inPure3DCorner && !inOutSlice && !inSideSlice && !inLongSlice){
        TF1::RejectPoint();
    }
    return Background3DValue(x[0], x[1], x[2], par);
}

DLM_Histo<double>* MakeDlmFeedDownFractionFromLambdaTH1(const TH1* hLambda,
                                                        const FitCATS3DConfig& cfg)
{
    if(!hLambda) return nullptr;

    const unsigned nBins = hLambda->GetNbinsX();
    std::vector<double> edges(nBins+1);
    for(unsigned iBin=0; iBin<=nBins; ++iBin){
        edges[iBin] = hLambda->GetXaxis()->GetBinLowEdge(iBin+1) *
                      cfg.lambdaAxisScaleToMeV;
    }

    DLM_Histo<double>* hFeedDown = new DLM_Histo<double>();
    hFeedDown->SetUp(1);
    hFeedDown->SetUp(0, nBins, edges.data());
    hFeedDown->Initialize(true);

    unsigned whichBin[1];
    const double weakScale = (1.0 - cfg.weakDecayFraction) *
                             (1.0 - cfg.weakDecayFraction);

    for(unsigned iBin=0; iBin<nBins; ++iBin){
        const double lambdaMain = hLambda->GetBinContent(iBin+1) * weakScale;
        const double lambdaErr = hLambda->GetBinError(iBin+1) * weakScale;

        if(!std::isfinite(lambdaMain) || lambdaMain<0.0 || lambdaMain>1.0){
            std::cout << "ERROR: lambda histogram has invalid value in bin "
                      << iBin+1 << ": " << lambdaMain << "\n";
            delete hFeedDown;
            return nullptr;
        }

        whichBin[0] = iBin;
        hFeedDown->SetBinContent(whichBin, 1.0 - lambdaMain);
        hFeedDown->SetBinError(whichBin, lambdaErr);
    }

    return hFeedDown;
}

DLM_Histo<double>* LoadDlmFeedDownFraction(const FitCATS3DConfig& cfg)
{
    if(!cfg.UseLambdaFeedDown) return nullptr;

    TFile fin(cfg.LambdaFileName, "READ");
    if(fin.IsZombie()){
        throw std::runtime_error(Form("LoadLambdaHistogram: cannot open '%s'",
                                      cfg.LambdaFileName.Data()));
    }

    TH1* hIn = dynamic_cast<TH1*>(fin.Get(cfg.LambdaHistoName));
    if(!hIn){
        throw std::runtime_error(Form("LoadLambdaHistogram: missing TH1 '%s'",
                                      cfg.LambdaHistoName.Data()));
    }

    return MakeDlmFeedDownFractionFromLambdaTH1(hIn, cfg);
}

double FeedDownAtKstar(const FitState3D& st, const double kstarMeV)
{
    if(!st.cfg.UseLambdaFeedDown || !st.LambdaFeedDown) return 0.0;

    double evalAt[1] = {kstarMeV};
    double feedDown = st.LambdaFeedDown->Eval(evalAt);

    if(!std::isfinite(feedDown)) feedDown = 0.0;
    if(feedDown < 0.0) feedDown = 0.0;
    if(feedDown > 1.0) feedDown = 1.0;

    return feedDown;
}

double LambdaMainAtKstar(const FitState3D& st, const double kstarMeV)
{
    return 1.0 - FeedDownAtKstar(st, kstarMeV);
}

void PrintLambdaFeedDownSummary(const FitState3D& st)
{
    if(!st.cfg.UseLambdaFeedDown || !st.LambdaFeedDown) return;

    const double weakScale = (1.0 - st.cfg.weakDecayFraction) *
                             (1.0 - st.cfg.weakDecayFraction);
    const double sampleKstarMeV[] = {
        0.0, 20.0, 40.0, 60.0, 80.0, 100.0,
        120.0, 150.0, 180.0, 200.0, 250.0, 300.0
    };

    std::cout << "CATS3D lambda feeddown summary\n"
              << "  input histogram: " << st.cfg.LambdaFileName
              << " :: " << st.cfg.LambdaHistoName << "\n"
              << "  weak primary scale: (1 - "
              << st.cfg.weakDecayFraction << ")^2 = "
              << weakScale << "\n"
              << "  k* [MeV]   lambdaMain   feedDown   inputPurityApprox\n";

    for(const double kstar : sampleKstarMeV){
        const double feedDown = FeedDownAtKstar(st, kstar);
        const double lambdaMain = 1.0 - feedDown;
        const double inputPurityApprox =
            weakScale > 0.0 ? lambdaMain / weakScale : 0.0;
        std::cout << "  " << std::setw(8) << kstar
                  << "   " << std::setw(10) << lambdaMain
                  << "   " << std::setw(8) << feedDown
                  << "   " << std::setw(16) << inputPurityApprox
                  << "\n";
    }
}

double NormalizeCk3D(TH3D* h, const FitCATS3DConfig& cfg)
{
    double sumW = 0.0;
    double sumWC = 0.0;

    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                const double kstar = RelativeK3MeV(h, ix, iy, iz, cfg);
                if(kstar < cfg.normQminMeV || kstar > cfg.normQmaxMeV) continue;

                const double c = h->GetBinContent(ix, iy, iz);
                const double e = h->GetBinError(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;

                const double w = (e > 0.0 && std::isfinite(e)) ? 1.0/(e*e) : 1.0;
                sumW += w;
                sumWC += w*c;
            }
        }
    }

    if(sumW <= 0.0){
        std::cout << "WARNING NormalizeCk3D: no bins in normalization window; data is left unscaled.\n";
        return 1.0;
    }

    const double norm = sumWC / sumW;
    if(norm > 0.0 && std::isfinite(norm)){
        h->Scale(1.0 / norm);
        return norm;
    }
    return 1.0;
}

double MeanInNormWindow3D(const TH3D& h,
                          const TH3D* weightRef,
                          const FitCATS3DConfig& cfg)
{
    double sumW = 0.0;
    double sumWC = 0.0;

    for(int ix=1; ix<=h.GetNbinsX(); ++ix){
        for(int iy=1; iy<=h.GetNbinsY(); ++iy){
            for(int iz=1; iz<=h.GetNbinsZ(); ++iz){
                const double kstar = RelativeK3MeV(&h, ix, iy, iz, cfg);
                if(kstar < cfg.normQminMeV || kstar > cfg.normQmaxMeV) continue;

                const double c = h.GetBinContent(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;

                double w = 1.0;
                if(weightRef){
                    const double e = weightRef->GetBinError(ix, iy, iz);
                    if(e > 0.0 && std::isfinite(e)) w = 1.0/(e*e);
                }

                sumW += w;
                sumWC += w*c;
            }
        }
    }

    if(sumW <= 0.0) return 0.0;
    return sumWC / sumW;
}

void NormalizeModelCk3D(TH3D* model,
                        const TH3D& data,
                        const FitCATS3DConfig& cfg)
{
    if(!model || !cfg.NormalizeModelInNormWindow) return;

    const double norm = MeanInNormWindow3D(*model, &data, cfg);
    if(norm > 0.0 && std::isfinite(norm)){
        model->Scale(1.0 / norm);
    }
    else {
        std::cout << "WARNING NormalizeModelCk3D: no bins in normalization window; model is left unscaled.\n";
    }
}

void PrintAxisMirrorCheck(const TAxis* axis,
                          const TString& label,
                          const double scaleToMeV)
{
    if(!axis) return;

    double maxAbsPairSum = 0.0;
    int worstLow = -1;
    int worstHigh = -1;
    const int nBins = axis->GetNbins();
    const int nPairs = nBins / 2;

    std::cout << "CATS3D q-axis mirror check " << label
              << ": nbins=" << nBins
              << " range=[" << axis->GetXmin() << ", " << axis->GetXmax()
              << "] centers in MeV\n";

    for(int iPair=1; iPair<=nPairs; ++iPair){
        const int iLow = iPair;
        const int iHigh = nBins + 1 - iPair;
        const double low = axis->GetBinCenter(iLow) * scaleToMeV;
        const double high = axis->GetBinCenter(iHigh) * scaleToMeV;
        const double pairSum = low + high;
        if(std::fabs(pairSum) > maxAbsPairSum){
            maxAbsPairSum = std::fabs(pairSum);
            worstLow = iLow;
            worstHigh = iHigh;
        }
        if(iPair <= 3 || iPair > nPairs - 3){
            std::cout << "  pair bins " << iLow << "/" << iHigh
                      << ": " << low << " + " << high
                      << " = " << pairSum << " MeV\n";
        }
    }

    if(nBins % 2 == 1){
        const int centerBin = (nBins + 1) / 2;
        std::cout << "  center bin " << centerBin
                  << ": " << axis->GetBinCenter(centerBin) * scaleToMeV
                  << " MeV\n";
    }

    std::cout << "  worst mirror mismatch: "
              << maxAbsPairSum << " MeV";
    if(worstLow > 0){
        std::cout << " at bins " << worstLow << "/" << worstHigh;
    }
    std::cout << "\n";
}

int MirrorBin(const TAxis* axis, const int bin)
{
    if(!axis) return bin;
    const double center = axis->GetBinCenter(bin);
    return axis->FindBin(-center);
}

TH3D* BuildParityDiff3D(const TH3D& h,
                        const TString& name,
                        const TString& mode)
{
    TH3D* diff = dynamic_cast<TH3D*>(h.Clone(name));
    diff->SetDirectory(nullptr);
    diff->Reset("ICES");

    for(int ix=1; ix<=h.GetNbinsX(); ++ix){
        const int mx = mode.Contains("out") || mode.EqualTo("full")
                     ? MirrorBin(h.GetXaxis(), ix) : ix;
        for(int iy=1; iy<=h.GetNbinsY(); ++iy){
            const int my = mode.Contains("side") || mode.EqualTo("full")
                         ? MirrorBin(h.GetYaxis(), iy) : iy;
            for(int iz=1; iz<=h.GetNbinsZ(); ++iz){
                const int mz = mode.Contains("long") || mode.EqualTo("full")
                             ? MirrorBin(h.GetZaxis(), iz) : iz;
                const double value = h.GetBinContent(ix, iy, iz) -
                                     h.GetBinContent(mx, my, mz);
                diff->SetBinContent(ix, iy, iz, value);
            }
        }
    }

    return diff;
}

void PrintModelParityCheck(const TH3D& h,
                           const TString& label,
                           const FitCATS3DConfig& cfg)
{
    if(!cfg.DebugQSymmetry) return;

    struct ParityMode {
        const char* name;
        const char* mode;
    };
    const ParityMode modes[] = {
        {"out", "out"},
        {"side", "side"},
        {"long", "long"},
        {"full", "full"}
    };

    for(const ParityMode& parityMode : modes){
        double maxAbsDiff = 0.0;
        double maxRelDiff = 0.0;
        int wx = -1;
        int wy = -1;
        int wz = -1;
        int wmx = -1;
        int wmy = -1;
        int wmz = -1;

        for(int ix=1; ix<=h.GetNbinsX(); ++ix){
            const int mx = TString(parityMode.mode).Contains("out") ||
                           TString(parityMode.mode).EqualTo("full")
                         ? MirrorBin(h.GetXaxis(), ix) : ix;
            for(int iy=1; iy<=h.GetNbinsY(); ++iy){
                const int my = TString(parityMode.mode).Contains("side") ||
                               TString(parityMode.mode).EqualTo("full")
                             ? MirrorBin(h.GetYaxis(), iy) : iy;
                for(int iz=1; iz<=h.GetNbinsZ(); ++iz){
                    const int mz = TString(parityMode.mode).Contains("long") ||
                                   TString(parityMode.mode).EqualTo("full")
                                 ? MirrorBin(h.GetZaxis(), iz) : iz;
                    const double c = h.GetBinContent(ix, iy, iz);
                    const double cm = h.GetBinContent(mx, my, mz);
                    const double diff = c - cm;
                    const double absDiff = std::fabs(diff);
                    const double denom = 0.5 * (std::fabs(c) + std::fabs(cm));
                    const double relDiff = denom > 0.0 ? absDiff / denom : 0.0;
                    if(absDiff > maxAbsDiff){
                        maxAbsDiff = absDiff;
                        maxRelDiff = relDiff;
                        wx = ix;
                        wy = iy;
                        wz = iz;
                        wmx = mx;
                        wmy = my;
                        wmz = mz;
                    }
                }
            }
        }

        std::cout << "CATS3D model parity check " << label
                  << " mode=" << parityMode.name
                  << ": maxAbsDiff=" << maxAbsDiff
                  << " maxRelDiff=" << maxRelDiff;
        if(wx > 0){
            std::cout << " at bin (" << wx << "," << wy << "," << wz << ")"
                      << " q=(" << h.GetXaxis()->GetBinCenter(wx)
                      << "," << h.GetYaxis()->GetBinCenter(wy)
                      << "," << h.GetZaxis()->GetBinCenter(wz) << ")"
                      << " mirror bin (" << wmx << "," << wmy << "," << wmz << ")"
                      << " mirror q=(" << h.GetXaxis()->GetBinCenter(wmx)
                      << "," << h.GetYaxis()->GetBinCenter(wmy)
                      << "," << h.GetZaxis()->GetBinCenter(wmz) << ")";
        }
        std::cout << "\n";
    }
}

void PrefitBackground3D(FitState3D& state)
{
    const FitCATS3DConfig& cfg = state.cfg;
    if(!cfg.UseBackground || !cfg.PrefitBackground || !state.data) return;

    TH3D* data = state.data;
    gBkgRejectQminMeV = cfg.BackgroundNonFemtoQminMeV;
    gBkgQAxisScaleToMeV = cfg.qAxisScaleToMeV;
    gBkgSliceCutMeV = cfg.sliceCutAbs * cfg.qAxisScaleToMeV;

    TF3 baseline("fCATS3D_background_prefit",
                 Background3DRejectFemtoTF3,
                 data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(),
                 data->GetYaxis()->GetXmin(), data->GetYaxis()->GetXmax(),
                 data->GetZaxis()->GetXmin(), data->GetZaxis()->GetXmax(),
                 10);
    const char* parNames[10] = {
        "BkgP0", "BkgOut", "BkgSide", "BkgLong",
        "BkgOut2", "BkgSide2", "BkgLong2",
        "BkgOutSide", "BkgOutLong", "BkgSideLong"
    };
    for(int iPar=0; iPar<10; ++iPar){
        baseline.SetParName(iPar, parNames[iPar]);
        baseline.SetParameter(iPar, state.BkgParInit[iPar]);
    }
    baseline.SetParLimits(0, 0.5, 1.5);

    const int fitStatus = data->Fit(&baseline, "RQ0");
    for(int iPar=0; iPar<10; ++iPar){
        state.BkgParInit[iPar] = baseline.GetParameter(iPar);
    }
    state.BkgNormInit = state.BkgParInit[0];

    std::cout << "CATS3D background prefit: "
              << "pol2 p0=" << state.BkgParInit[0]
              << " p1=" << state.BkgParInit[1]
              << " p2=" << state.BkgParInit[2]
              << " p3=" << state.BkgParInit[3]
              << " status=" << fitStatus
              << " region hybrid non-femto: pure 3D corner |q_all| >= "
              << cfg.BackgroundNonFemtoQminMeV
              << " MeV OR projection slice |q_target| >= "
              << cfg.BackgroundNonFemtoQminMeV
              << " MeV with other two |q| <= "
              << gBkgSliceCutMeV << " MeV\n";
}

void ConfigureKittyForData(CATS3D& kitty,
                           const TH3D& data,
                           const FitCATS3DConfig& cfg,
                           DLM_CleverReso3D* resonanceSource)
{
    if(cfg.DebugQSymmetry){
        PrintAxisMirrorCheck(data.GetXaxis(), "data q_out", cfg.qAxisScaleToMeV);
        PrintAxisMirrorCheck(data.GetYaxis(), "data q_side", cfg.qAxisScaleToMeV);
        PrintAxisMirrorCheck(data.GetZaxis(), "data q_long", cfg.qAxisScaleToMeV);
    }

    kitty.SetMomBinsOut(data.GetNbinsX(),
                        data.GetXaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                        data.GetXaxis()->GetXmax()*cfg.qAxisScaleToMeV);
    kitty.SetMomBinsSide(data.GetNbinsY(),
                         data.GetYaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                         data.GetYaxis()->GetXmax()*cfg.qAxisScaleToMeV);
    kitty.SetMomBinsLong(data.GetNbinsZ(),
                         data.GetZaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                         data.GetZaxis()->GetXmax()*cfg.qAxisScaleToMeV);

    kitty.SetGridX(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridY(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridZ(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetFastKernelMode(cfg.UseFastKernel);
    kitty.SetFastKernelSubdivisions(cfg.FastKernelSubdivisions);

    CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
    sourcePars.SetParameter(0, cfg.RoutInit);
    sourcePars.SetParameter(1, cfg.RsideInit);
    sourcePars.SetParameter(2, cfg.RlongInit);
    if(cfg.UseResonanceSource){
        if(!resonanceSource){
            throw std::runtime_error("ConfigureKittyForData: resonance source requested but not built.");
        }
        kitty.SetAnaSource(CatsSourceForwarder, resonanceSource, 4);
        kitty.SetAnaSource(0, cfg.RoutInit);
        kitty.SetAnaSource(1, cfg.RsideInit);
        kitty.SetAnaSource(2, cfg.RlongInit);
        kitty.SetAnaSource(3, cfg.ResonanceAlphaInit);
        std::cout << "CATS3D source mode: resonance"
                  << " weights=(" << cfg.ResonanceWeight1
                  << ", " << cfg.ResonanceWeight2 << ")"
                  << " alphaInit=" << cfg.ResonanceAlphaInit
                  << " fastKernel=" << (cfg.UseFastKernel ? "on" : "off")
                  << " fastSubdiv=" << cfg.FastKernelSubdivisions
                  << "\n";
    }
    else {
        kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);
        std::cout << "CATS3D source mode: Gaussian"
                  << " fastKernel=" << (cfg.UseFastKernel ? "on" : "off")
                  << " fastSubdiv=" << cfg.FastKernelSubdivisions
                  << "\n";
    }

    if(cfg.Identical) kitty.SetPdgId(211, 211);
    else              kitty.SetPdgId(-211, 211);

    kitty.SetQuantumStatistics(cfg.Identical ? 1 : 0);
    kitty.SetQ1Q2(cfg.WithCoulomb ? 1 : 0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(cfg.WithCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
}

TH3D* BuildModelFromKitty(const TH3D& data,
                          CATS3D& kitty,
                          const TString& name)
{
    FitState3D& st = *gFitState3D;
    TH3D* model = dynamic_cast<TH3D*>(data.Clone(name));
    model->SetDirectory(nullptr);
    model->Reset("ICES");

    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                const double cCore = kitty.GetCorrFun(ix-1, iy-1, iz-1);
                const double eCore = kitty.GetCorrFunErr(ix-1, iy-1, iz-1);
                const double kstar = RelativeK3MeV(&data, ix, iy, iz, st.cfg);
                const double lambdaMain = LambdaMainAtKstar(st, kstar);
                const double cTotal = 1.0 + lambdaMain * (cCore - 1.0);

                model->SetBinContent(ix, iy, iz,
                                     cTotal);
                model->SetBinError(ix, iy, iz,
                                   lambdaMain * eCore);
            }
        }
    }

    NormalizeModelCk3D(model, data, st.cfg);
    // PrintModelParityCheck(*model, name, st.cfg);
    return model;
}

bool SameModelCachePoint(const ModelCacheEntry3D& entry,
                         const double Rout,
                         const double Rside,
                         const double Rlong,
                         const double Alpha,
                         const double tolerance)
{
    return std::fabs(entry.Rout  - Rout)  <= tolerance &&
           std::fabs(entry.Rside - Rside) <= tolerance &&
           std::fabs(entry.Rlong - Rlong) <= tolerance &&
           std::fabs(entry.Alpha - Alpha) <= 0.02;
}

TH3D* GetCachedModelForPars(const double* par)
{
    FitState3D& st = *gFitState3D;
    CATS3D& kitty = *st.kitty;
    TH3D& data = *st.data;

    const double Rout = par[0];
    const double Rside = par[1];
    const double Rlong = par[2];
    const double Alpha = st.cfg.UseResonanceSource ? par[3] : 2.0;

    ++st.CacheClock;

    const auto tLookup0 = Clock3D::now();
    if(st.cfg.UseModelCache){
        for(ModelCacheEntry3D& entry : st.ModelCache){
            if(SameModelCachePoint(entry, Rout, Rside, Rlong, Alpha,
                                   st.cfg.ModelCacheToleranceFm)){
                st.Profile.CacheLookupSec += SecondsSince(tLookup0);
                ++st.Profile.CacheLookupCalls;
                entry.LastUsed = st.CacheClock;
                ++st.nCacheHits;
                if(st.nFcnCalls <= 3 || st.nFcnCalls % 10 == 0){
                    std::cout << "CATS3D cache hit " << st.nCacheHits
                              << ": requested=(" << Rout << ", " << Rside << ", " << Rlong << ")"
                              << " alpha=" << Alpha
                              << " cached=(" << entry.Rout << ", " << entry.Rside << ", " << entry.Rlong
                              << ") alpha=" << entry.Alpha << "\n";
                }
                return entry.Model.get();
            }
        }
    }
    st.Profile.CacheLookupSec += SecondsSince(tLookup0);
    ++st.Profile.CacheLookupCalls;

    ++st.nCatsEvals;
    const auto tSet0 = Clock3D::now();
    kitty.SetAnaSource(0, Rout, true);
    kitty.SetAnaSource(1, Rside, true);
    kitty.SetAnaSource(2, Rlong, true);
    if(st.cfg.UseResonanceSource){
        kitty.SetAnaSource(3, Alpha, true);
    }
    st.Profile.SetAnaSourceSec += SecondsSince(tSet0);
    ++st.Profile.SetAnaSourceCalls;

    const bool printTiming = st.nCatsEvals <= 3 || st.nCatsEvals % 10 == 0;
    const auto t0 = Clock3D::now();
    if(printTiming){
        std::cout << "CATS3D build " << st.nCatsEvals
                  << ": Rout=" << Rout
                  << " Rside=" << Rside
                  << " Rlong=" << Rlong
                  << " alpha=" << Alpha
                  << " ... " << std::flush;
    }
    kitty.KillTheCat();
    const double killSec = SecondsSince(t0);
    st.Profile.KillTheCatSec += killSec;
    ++st.Profile.KillTheCatCalls;
    if(printTiming){
        std::cout << killSec << " s";
    }

    const auto tModel0 = Clock3D::now();
    std::unique_ptr<TH3D> newModel(BuildModelFromKitty(data, kitty,
                                                       Form("hCATS3D_cache_%u", st.CacheClock)));
    const double buildModelSec = SecondsSince(tModel0);
    st.Profile.BuildModelSec += buildModelSec;
    ++st.Profile.BuildModelCalls;
    if(printTiming && st.cfg.EnableProfiling){
        std::cout << " (KillTheCat=" << killSec
                  << " s, BuildModel=" << buildModelSec
                  << " s)\n";
    }
    else if(printTiming){
        std::cout << "\n";
    }

    if(!st.cfg.UseModelCache || st.cfg.MaxModelCacheEntries == 0){
        st.ScratchModel = std::move(newModel);
        return st.ScratchModel.get();
    }

    ModelCacheEntry3D newEntry;
    newEntry.Rout = Rout;
    newEntry.Rside = Rside;
    newEntry.Rlong = Rlong;
    newEntry.Alpha = Alpha;
    newEntry.LastUsed = st.CacheClock;
    newEntry.Model = std::move(newModel);

    if(st.ModelCache.size() < st.cfg.MaxModelCacheEntries){
        st.ModelCache.push_back(std::move(newEntry));
        return st.ModelCache.back().Model.get();
    }

    unsigned evict = 0;
    for(unsigned iEntry=1; iEntry<st.ModelCache.size(); ++iEntry){
        if(st.ModelCache[iEntry].LastUsed < st.ModelCache[evict].LastUsed){
            evict = iEntry;
        }
    }
    st.ModelCache[evict] = std::move(newEntry);
    return st.ModelCache[evict].Model.get();
}

double ComputeChi2ForPars(const double* par)
{
    FitState3D& st = *gFitState3D;
    TH3D& data = *st.data;
    ++st.nFcnCalls;
    ++st.Profile.Chi2Calls;

    const auto tModel0 = Clock3D::now();
    TH3D* catsModel = GetCachedModelForPars(par);
    st.Profile.Chi2GetModelSec += SecondsSince(tModel0);

    double chi2 = 0.0;
    int nBins = 0;
    const bool usePml = st.cfg.UsePmlChi2 && st.same && st.mixed;

    const auto tLoop0 = Clock3D::now();
    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                if(!IsInsideFitBoxMeV(&data, ix, iy, iz, st.cfg)) continue;

                double background = par[4];
                if(st.cfg.UseBackground){
                    background = Background3DValue(data, ix, iy, iz, par);
                }
                if(!std::isfinite(background) || background <= 0.0) continue;

                const double model = background * catsModel->GetBinContent(ix, iy, iz);
                if(!std::isfinite(model) || model <= 0.0) continue;

                if(usePml){
                    const double same = st.same->GetBinContent(ix, iy, iz);
                    const double mixed = st.mixed->GetBinContent(ix, iy, iz) *
                                         st.DataNormFactor;
                    if(!std::isfinite(same) || !std::isfinite(mixed)) continue;
                    if(same < 0.0 || mixed < 0.0) continue;
                    const double total = same + mixed;
                    if(total <= 0.0) continue;

                    double logLike = 0.0;
                    if(same > 0.0){
                        const double argSame =
                            model * total / (same * (model + 1.0));
                        if(argSame <= 0.0 || !std::isfinite(argSame)) continue;
                        logLike += same * std::log(argSame);
                    }
                    if(mixed > 0.0){
                        const double argMixed =
                            total / (mixed * (model + 1.0));
                        if(argMixed <= 0.0 || !std::isfinite(argMixed)) continue;
                        logLike += mixed * std::log(argMixed);
                    }
                    chi2 += -2.0 * logLike;
                }
                else {
                    const double y = data.GetBinContent(ix, iy, iz);
                    const double e = data.GetBinError(ix, iy, iz);
                    if(!std::isfinite(y) || !std::isfinite(e) || e <= 0.0) continue;
                    const double pull = (y - model) / e;
                    chi2 += pull*pull;
                }
                ++nBins;
            }
        }
    }
    st.Profile.Chi2LoopSec += SecondsSince(tLoop0);

    st.nFitBins = nBins;
    return chi2;
}

void FitFCN3D(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag)
{
    (void)npar;
    (void)gin;
    (void)iflag;

    if(!gFitState3D || !gFitState3D->data || !gFitState3D->kitty){
        f = 1.e30;
        return;
    }

    if(par[0] <= 0.0 || par[1] <= 0.0 || par[2] <= 0.0 || par[4] <= 0.0){
        f = 1.e30;
        return;
    }
    if(gFitState3D->cfg.UseResonanceSource &&
       (par[3] < gFitState3D->cfg.ResonanceAlphaMin ||
        par[3] > gFitState3D->cfg.ResonanceAlphaMax)){
        f = 1.e30;
        return;
    }

    f = ComputeChi2ForPars(par);
}

TH1D* ProjectFemtoSlice1D(TH3D* h3,
                          const TString& target,
                          const TString& name,
                          const double sliceCutAbs)
{
    if(!h3){
        throw std::runtime_error("ProjectFemtoSlice1D: input histogram is null.");
    }

    TH1D* h1 = nullptr;
    auto centeredBinRange = [](TAxis* axis, const double cutAbs, int& firstBin, int& lastBin) {
        firstBin = 1;
        lastBin = 0;
        for(int iBin=1; iBin<=axis->GetNbins(); ++iBin){
            if(std::fabs(axis->GetBinCenter(iBin)) <= cutAbs + 1e-12){
                if(lastBin == 0) firstBin = iBin;
                lastBin = iBin;
            }
        }
        if(lastBin == 0){
            const int centerBin = axis->FindBin(0.0);
            firstBin = centerBin;
            lastBin = centerBin;
        }
    };
    auto printSliceAxisCheck = [](TAxis* axis,
                                  const TString& axisName,
                                  const TString& histName,
                                  const int firstBin,
                                  const int lastBin) {
        if(!gFitState3D || !gFitState3D->cfg.DebugQSymmetry || !axis) return;
        static int nPrint = 0;
        if(nPrint >= 24) return;
        ++nPrint;

        std::cout << "CATS3D projection slice mirror check "
                  << histName << " axis=" << axisName
                  << " bins " << firstBin << "-" << lastBin
                  << " centers:";
        for(int iBin=firstBin; iBin<=lastBin; ++iBin){
            std::cout << " " << axis->GetBinCenter(iBin);
        }
        std::cout << "\n";

        for(int iBin=firstBin; iBin<=lastBin; ++iBin){
            const double center = axis->GetBinCenter(iBin);
            const int mirrorBin = axis->FindBin(-center);
            const bool mirrorInside = mirrorBin >= firstBin && mirrorBin <= lastBin;
            const double mirrorCenter = axis->GetBinCenter(mirrorBin);
            std::cout << "  bin " << iBin
                      << " q=" << center
                      << " mirrorBin=" << mirrorBin
                      << " mirrorQ=" << mirrorCenter
                      << " q+mirrorQ=" << center + mirrorCenter
                      << " insideSlice=" << (mirrorInside ? "yes" : "no")
                      << "\n";
        }
    };

    if(target.EqualTo("out")){
        int y1, y2, z1, z2;
        centeredBinRange(h3->GetYaxis(), sliceCutAbs, y1, y2);
        centeredBinRange(h3->GetZaxis(), sliceCutAbs, z1, z2);
        // printSliceAxisCheck(h3->GetYaxis(), "q_side", name, y1, y2);
        // printSliceAxisCheck(h3->GetZaxis(), "q_long", name, z1, z2);

        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, y1, y2, z1, z2, "e"));
        const int numBins = (y2-y1+1) * (z2-z1+1);
        if(numBins > 0) h1->Scale(1.0 / double(numBins));
    }
    else if(target.EqualTo("side")){
        int x1, x2, z1, z2;
        centeredBinRange(h3->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(h3->GetZaxis(), sliceCutAbs, z1, z2);
        // printSliceAxisCheck(h3->GetXaxis(), "q_out", name, x1, x2);
        // printSliceAxisCheck(h3->GetZaxis(), "q_long", name, z1, z2);

        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, x1, x2, z1, z2, "e"));
        const int numBins = (x2-x1+1) * (z2-z1+1);
        if(numBins > 0) h1->Scale(1.0 / double(numBins));
    }
    else if(target.EqualTo("long")){
        int x1, x2, y1, y2;
        centeredBinRange(h3->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(h3->GetYaxis(), sliceCutAbs, y1, y2);
        // printSliceAxisCheck(h3->GetXaxis(), "q_out", name, x1, x2);
        // printSliceAxisCheck(h3->GetYaxis(), "q_side", name, y1, y2);

        h1 = dynamic_cast<TH1D*>(h3->ProjectionZ(name, x1, x2, y1, y2, "e"));
        const int numBins = (x2-x1+1) * (y2-y1+1);
        if(numBins > 0) h1->Scale(1.0 / double(numBins));
    }
    else {
        throw std::runtime_error(Form("ProjectFemtoSlice1D: unknown target '%s'", target.Data()));
    }

    if(!h1){
        throw std::runtime_error(Form("ProjectFemtoSlice1D: failed for target '%s'", target.Data()));
    }

    h1->SetDirectory(nullptr);
    return h1;
}

void WriteFemtoSliceProjections(TH3D* h3,
                                const TString& baseName,
                                const double sliceCutAbs)
{
    std::unique_ptr<TH1D> hOut(ProjectFemtoSlice1D(h3, "out",
                                                   baseName + "_proj_out",
                                                   sliceCutAbs));
    std::unique_ptr<TH1D> hSide(ProjectFemtoSlice1D(h3, "side",
                                                    baseName + "_proj_side",
                                                    sliceCutAbs));
    std::unique_ptr<TH1D> hLong(ProjectFemtoSlice1D(h3, "long",
                                                    baseName + "_proj_long",
                                                    sliceCutAbs));

    hOut->Write("", TObject::kOverwrite);
    hSide->Write("", TObject::kOverwrite);
    hLong->Write("", TObject::kOverwrite);
}

TH1D* ProjectMixedWeightedSlice1D(TH3D* hCk,
                                  TH3D* hWeight,
                                  const TString& target,
                                  const TString& name,
                                  const double sliceCutAbs)
{
    if(!hCk || !hWeight){
        throw std::runtime_error("ProjectMixedWeightedSlice1D: null input histogram.");
    }

    std::unique_ptr<TH1D> hWeightedSum;
    std::unique_ptr<TH1D> hWeightSum;

    auto centeredBinRange = [](TAxis* axis, const double cutAbs, int& firstBin, int& lastBin) {
        firstBin = 1;
        lastBin = 0;
        for(int iBin=1; iBin<=axis->GetNbins(); ++iBin){
            if(std::fabs(axis->GetBinCenter(iBin)) <= cutAbs + 1e-12){
                if(lastBin == 0) firstBin = iBin;
                lastBin = iBin;
            }
        }
        if(lastBin == 0){
            const int centerBin = axis->FindBin(0.0);
            firstBin = centerBin;
            lastBin = centerBin;
        }
    };

    TH3D* hWeighted3D = dynamic_cast<TH3D*>(hCk->Clone(name + "_weighted3D"));
    hWeighted3D->SetDirectory(nullptr);
    hWeighted3D->Multiply(hWeight);

    if(target.EqualTo("out")){
        int y1, y2, z1, z2;
        centeredBinRange(hCk->GetYaxis(), sliceCutAbs, y1, y2);
        centeredBinRange(hCk->GetZaxis(), sliceCutAbs, z1, z2);
        hWeightedSum.reset(dynamic_cast<TH1D*>(hWeighted3D->ProjectionX(name + "_weightedSum", y1, y2, z1, z2, "e")));
        hWeightSum.reset(dynamic_cast<TH1D*>(hWeight->ProjectionX(name + "_weightSum", y1, y2, z1, z2, "e")));
    }
    else if(target.EqualTo("side")){
        int x1, x2, z1, z2;
        centeredBinRange(hCk->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(hCk->GetZaxis(), sliceCutAbs, z1, z2);
        hWeightedSum.reset(dynamic_cast<TH1D*>(hWeighted3D->ProjectionY(name + "_weightedSum", x1, x2, z1, z2, "e")));
        hWeightSum.reset(dynamic_cast<TH1D*>(hWeight->ProjectionY(name + "_weightSum", x1, x2, z1, z2, "e")));
    }
    else if(target.EqualTo("long")){
        int x1, x2, y1, y2;
        centeredBinRange(hCk->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(hCk->GetYaxis(), sliceCutAbs, y1, y2);
        hWeightedSum.reset(dynamic_cast<TH1D*>(hWeighted3D->ProjectionZ(name + "_weightedSum", x1, x2, y1, y2, "e")));
        hWeightSum.reset(dynamic_cast<TH1D*>(hWeight->ProjectionZ(name + "_weightSum", x1, x2, y1, y2, "e")));
    }
    else {
        delete hWeighted3D;
        throw std::runtime_error(Form("ProjectMixedWeightedSlice1D: unknown target '%s'", target.Data()));
    }

    delete hWeighted3D;

    if(!hWeightedSum || !hWeightSum){
        throw std::runtime_error(Form("ProjectMixedWeightedSlice1D: failed for target '%s'", target.Data()));
    }

    hWeightedSum->SetDirectory(nullptr);
    hWeightSum->SetDirectory(nullptr);

    TH1D* hOut = dynamic_cast<TH1D*>(hWeightedSum->Clone(name));
    hOut->SetDirectory(nullptr);
    hOut->Divide(hWeightedSum.get(), hWeightSum.get(), 1.0, 1.0);
    return hOut;
}

void WriteMixedWeightedFemtoSliceProjections(TH3D* hCk,
                                             TH3D* hWeight,
                                             const TString& baseName,
                                             const double sliceCutAbs)
{
    if(!hCk || !hWeight) return;

    std::unique_ptr<TH1D> hOut(ProjectMixedWeightedSlice1D(
        hCk, hWeight, "out", baseName + "_mixedWeighted_proj_out", sliceCutAbs));
    std::unique_ptr<TH1D> hSide(ProjectMixedWeightedSlice1D(
        hCk, hWeight, "side", baseName + "_mixedWeighted_proj_side", sliceCutAbs));
    std::unique_ptr<TH1D> hLong(ProjectMixedWeightedSlice1D(
        hCk, hWeight, "long", baseName + "_mixedWeighted_proj_long", sliceCutAbs));

    hOut->Write("", TObject::kOverwrite);
    hSide->Write("", TObject::kOverwrite);
    hLong->Write("", TObject::kOverwrite);
}

TH3D* LoadMixedWeightForProjection(const FitCATS3DConfig& cfg)
{
    if(!cfg.UseSameMixed) return nullptr;

    TFile fin(cfg.InputFileName, "READ");
    if(fin.IsZombie()){
        std::cout << "WARNING LoadMixedWeightForProjection: cannot open "
                  << cfg.InputFileName << "\n";
        return nullptr;
    }

    TH3D* hMixedIn = dynamic_cast<TH3D*>(fin.Get(cfg.MixedHistoName));
    if(!hMixedIn){
        std::cout << "WARNING LoadMixedWeightForProjection: missing "
                  << cfg.MixedHistoName << "\n";
        return nullptr;
    }

    TH3D* hMixed = RebinClone3D(hMixedIn, cfg.RebinFactor,
                                "hMixed3D_projectionWeight");
    return hMixed;
}

TH3D* BuildLambdaMainMap3D(const TH3D& data)
{
    FitState3D& st = *gFitState3D;
    if(!st.cfg.UseLambdaFeedDown || !st.LambdaFeedDown) return nullptr;

    TH3D* lambdaMap = dynamic_cast<TH3D*>(data.Clone("hLambdaMain3D"));
    lambdaMap->SetDirectory(nullptr);
    lambdaMap->Reset("ICES");

    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                const double kstar = RelativeK3MeV(&data, ix, iy, iz, st.cfg);
                lambdaMap->SetBinContent(ix, iy, iz, LambdaMainAtKstar(st, kstar));
            }
        }
    }

    return lambdaMap;
}

TH3D* BuildFeedDownMap3D(const TH3D& data)
{
    FitState3D& st = *gFitState3D;
    if(!st.cfg.UseLambdaFeedDown || !st.LambdaFeedDown) return nullptr;

    TH3D* feedDownMap = dynamic_cast<TH3D*>(data.Clone("hFeedDown3D"));
    feedDownMap->SetDirectory(nullptr);
    feedDownMap->Reset("ICES");

    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                const double kstar = RelativeK3MeV(&data, ix, iy, iz, st.cfg);
                feedDownMap->SetBinContent(ix, iy, iz, FeedDownAtKstar(st, kstar));
            }
        }
    }

    return feedDownMap;
}

void SaveFitOutput(const FitCATS3DConfig& cfg,
                   const FitCATS3DResult& result,
                   TH3D* data,
                   TH3D* rawModel)
{
    if(!rawModel){
        throw std::runtime_error("SaveFitOutput: rawModel is null.");
    }

    std::unique_ptr<TH3D> modelShape(dynamic_cast<TH3D*>(rawModel->Clone("hCATS3D_shape")));
    modelShape->SetDirectory(nullptr);

    std::unique_ptr<TH3D> model(dynamic_cast<TH3D*>(rawModel->Clone("hCATS3D_best")));
    model->SetDirectory(nullptr);

    std::unique_ptr<TH3D> background(dynamic_cast<TH3D*>(rawModel->Clone("hCATS3D_background")));
    background->SetDirectory(nullptr);
    background->Reset("ICES");

    const double pars[14] = {result.Rout, result.Rside, result.Rlong,
                             result.Alpha,
                             result.BkgPar[0], result.BkgPar[1],
                             result.BkgPar[2], result.BkgPar[3],
                             result.BkgPar[4], result.BkgPar[5],
                             result.BkgPar[6], result.BkgPar[7],
                             result.BkgPar[8], result.BkgPar[9]};

    for(int ix=1; ix<=model->GetNbinsX(); ++ix){
        for(int iy=1; iy<=model->GetNbinsY(); ++iy){
            for(int iz=1; iz<=model->GetNbinsZ(); ++iz){
                double bkg = result.Norm;
                if(cfg.UseBackground){
                    bkg = Background3DValue(*model, ix, iy, iz, pars);
                }
                background->SetBinContent(ix, iy, iz, bkg);
                model->SetBinContent(ix, iy, iz,
                                     model->GetBinContent(ix, iy, iz) * bkg);
                model->SetBinError(ix, iy, iz,
                                   model->GetBinError(ix, iy, iz) * std::fabs(bkg));
            }
        }
    }

    std::unique_ptr<TH3D> ratio(dynamic_cast<TH3D*>(data->Clone("hDataOverCATS3D")));
    ratio->SetDirectory(nullptr);
    ratio->Reset("ICES");

    for(int ix=1; ix<=data->GetNbinsX(); ++ix){
        for(int iy=1; iy<=data->GetNbinsY(); ++iy){
            for(int iz=1; iz<=data->GetNbinsZ(); ++iz){
                const double m = model->GetBinContent(ix, iy, iz);
                if(m == 0.0) continue;
                ratio->SetBinContent(ix, iy, iz, data->GetBinContent(ix, iy, iz)/m);
                ratio->SetBinError(ix, iy, iz, data->GetBinError(ix, iy, iz)/std::fabs(m));
            }
        }
    }

    TFile fout(cfg.OutputFileName, "RECREATE");
    data->Write("hDataCk3D");
    if(gFitState3D && gFitState3D->same && gFitState3D->mixed){
        gFitState3D->same->Write("hSame3D_forPML", TObject::kOverwrite);
        gFitState3D->mixed->Write("hMixed3D_forPML_raw", TObject::kOverwrite);
        std::unique_ptr<TH3D> mixedScaled(dynamic_cast<TH3D*>(
            gFitState3D->mixed->Clone("hMixed3D_forPML_scaledByDataNorm")));
        mixedScaled->SetDirectory(nullptr);
        mixedScaled->Scale(gFitState3D->DataNormFactor);
        mixedScaled->Write("", TObject::kOverwrite);
    }
    modelShape->Write("hCATS3D_shape");
    background->Write("hCATS3D_background");
    model->Write("hCATS3D_best");
    ratio->Write("hDataOverCATS3D");
    std::unique_ptr<TH3D> shapeDiffOut(BuildParityDiff3D(*modelShape, "hCATS3D_shape_parityDiff_out", "out"));
    std::unique_ptr<TH3D> shapeDiffSide(BuildParityDiff3D(*modelShape, "hCATS3D_shape_parityDiff_side", "side"));
    std::unique_ptr<TH3D> shapeDiffLong(BuildParityDiff3D(*modelShape, "hCATS3D_shape_parityDiff_long", "long"));
    std::unique_ptr<TH3D> shapeDiffFull(BuildParityDiff3D(*modelShape, "hCATS3D_shape_parityDiff_full", "full"));
    std::unique_ptr<TH3D> bestDiffOut(BuildParityDiff3D(*model, "hCATS3D_best_parityDiff_out", "out"));
    std::unique_ptr<TH3D> bestDiffSide(BuildParityDiff3D(*model, "hCATS3D_best_parityDiff_side", "side"));
    std::unique_ptr<TH3D> bestDiffLong(BuildParityDiff3D(*model, "hCATS3D_best_parityDiff_long", "long"));
    std::unique_ptr<TH3D> bestDiffFull(BuildParityDiff3D(*model, "hCATS3D_best_parityDiff_full", "full"));
    shapeDiffOut->Write("", TObject::kOverwrite);
    shapeDiffSide->Write("", TObject::kOverwrite);
    shapeDiffLong->Write("", TObject::kOverwrite);
    shapeDiffFull->Write("", TObject::kOverwrite);
    bestDiffOut->Write("", TObject::kOverwrite);
    bestDiffSide->Write("", TObject::kOverwrite);
    bestDiffLong->Write("", TObject::kOverwrite);
    bestDiffFull->Write("", TObject::kOverwrite);
    std::unique_ptr<TH3D> lambdaMap(BuildLambdaMainMap3D(*data));
    if(lambdaMap){
        lambdaMap->Write("hLambdaMain3D");
        WriteFemtoSliceProjections(lambdaMap.get(), "hLambdaMain3D", cfg.sliceCutAbs);
    }
    std::unique_ptr<TH3D> feedDownMap(BuildFeedDownMap3D(*data));
    if(feedDownMap){
        feedDownMap->Write("hFeedDown3D");
        WriteFemtoSliceProjections(feedDownMap.get(), "hFeedDown3D", cfg.sliceCutAbs);
    }
    WriteFemtoSliceProjections(data, "hDataCk3D", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(modelShape.get(), "hCATS3D_shape", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(background.get(), "hCATS3D_background", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(model.get(), "hCATS3D_best", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(ratio.get(), "hDataOverCATS3D", cfg.sliceCutAbs);

    std::unique_ptr<TH3D> mixedWeight(LoadMixedWeightForProjection(cfg));
    if(mixedWeight){
        mixedWeight->Write("hMixed3D_projectionWeight", TObject::kOverwrite);
        WriteMixedWeightedFemtoSliceProjections(data, mixedWeight.get(),
                                                "hDataCk3D", cfg.sliceCutAbs);
        WriteMixedWeightedFemtoSliceProjections(modelShape.get(), mixedWeight.get(),
                                                "hCATS3D_shape", cfg.sliceCutAbs);
        WriteMixedWeightedFemtoSliceProjections(model.get(), mixedWeight.get(),
                                                "hCATS3D_best", cfg.sliceCutAbs);
        std::cout << "CATS3D wrote mixed-weighted slice projections for data/model.\n";
    }

    TNtuple summary("fit_summary", "fit_summary",
                    "Rout:Rside:Rlong:Alpha:Norm:Chi2:Ndf:Chi2Ndf:NFitBins:"
                    "FitStatus:MinuitIstat:FitQuality:MinuitAmin:MinuitEdm:"
                    "MinuitErrDef:MinuitNvpar:MinuitNparx:"
                    "BkgP0:BkgOut:BkgSide:BkgLong:"
                    "BkgOut2:BkgSide2:BkgLong2:"
                    "BkgOutSide:BkgOutLong:BkgSideLong");
    Float_t summaryValues[27] = {
        float(result.Rout), float(result.Rside), float(result.Rlong),
        float(result.Alpha), float(result.Norm), float(result.Chi2),
        float(result.Ndf), float(result.Chi2Ndf()), float(result.NFitBins),
        float(result.FitStatus), float(result.MinuitIstat),
        float(result.FitQuality), float(result.MinuitAmin),
        float(result.MinuitEdm), float(result.MinuitErrDef),
        float(result.MinuitNvpar), float(result.MinuitNparx),
        float(result.BkgPar[0]), float(result.BkgPar[1]),
        float(result.BkgPar[2]), float(result.BkgPar[3]),
        float(result.BkgPar[4]), float(result.BkgPar[5]),
        float(result.BkgPar[6]), float(result.BkgPar[7]),
        float(result.BkgPar[8]), float(result.BkgPar[9])
    };
    summary.Fill(summaryValues);
    summary.Write();
}

void PrintFitProfile3D(const FitState3D& state)
{
    if(!state.cfg.EnableProfiling) return;

    const FitProfile3D& p = state.Profile;
    std::cout << "\nCATS3D profiling summary\n"
              << "  setup: data=" << p.BuildDataSec
              << " s, resonanceSource=" << p.BuildResonanceSourceSec
              << " s, configureKitty=" << p.ConfigureKittySec
              << " s, externalKernel=" << p.LoadExternalKernelSec
              << " s, lambda=" << p.LoadLambdaSec << " s\n"
              << "  minuit=" << p.MinuitSec
              << " s, finalChi2=" << p.FinalChi2Sec
              << " s, saveOutput=" << p.SaveOutputSec << " s\n"
              << "  cache lookup: " << p.CacheLookupSec
              << " s / " << p.CacheLookupCalls << " calls\n"
              << "  SetAnaSource: " << p.SetAnaSourceSec
              << " s / " << p.SetAnaSourceCalls << " calls\n"
              << "  KillTheCat: " << p.KillTheCatSec
              << " s / " << p.KillTheCatCalls << " calls";
    if(p.KillTheCatCalls>0){
        std::cout << " avg=" << p.KillTheCatSec/double(p.KillTheCatCalls) << " s";
    }
    std::cout << "\n"
              << "  BuildModelFromKitty: " << p.BuildModelSec
              << " s / " << p.BuildModelCalls << " calls\n"
              << "  chi2 getModel: " << p.Chi2GetModelSec
              << " s / " << p.Chi2Calls << " calls\n"
              << "  chi2 loop: " << p.Chi2LoopSec
              << " s / " << p.Chi2Calls << " calls\n";
}

TH3D* BuildDataCk3DInternal(const FitCATS3DConfig& cfg,
                            std::unique_ptr<TH3D>* sameForPml,
                            std::unique_ptr<TH3D>* mixedForPml,
                            double* normFactor)
{
    if(normFactor) *normFactor = 1.0;

    TFile fin(cfg.InputFileName, "READ");
    if(fin.IsZombie()){
        throw std::runtime_error(Form("BuildDataCk3D: cannot open '%s'", cfg.InputFileName.Data()));
    }

    TH3D* hCk = nullptr;

    if(cfg.UseSameMixed){
        std::unique_ptr<TH3D> hSame(GetTH3DClone(fin, cfg.SameHistoName, "hSame3D_forCk"));
        std::unique_ptr<TH3D> hMixed(GetTH3DClone(fin, cfg.MixedHistoName, "hMixed3D_forCk"));

        std::unique_ptr<TH3D> hSameFit(RebinClone3D(hSame.get(), cfg.RebinFactor,
                                                    "hSame3D_forCk_rebin"));
        std::unique_ptr<TH3D> hMixedFit(RebinClone3D(hMixed.get(), cfg.RebinFactor,
                                                     "hMixed3D_forCk_rebin"));

        if(sameForPml){
            sameForPml->reset(dynamic_cast<TH3D*>(hSameFit->Clone("hSame3D_forPML")));
            (*sameForPml)->SetDirectory(nullptr);
        }
        if(mixedForPml){
            mixedForPml->reset(dynamic_cast<TH3D*>(hMixedFit->Clone("hMixed3D_forPML")));
            (*mixedForPml)->SetDirectory(nullptr);
        }

        hCk = dynamic_cast<TH3D*>(hSameFit->Clone("hDataCk3D"));
        hCk->SetDirectory(nullptr);
        hCk->Divide(hSameFit.get(), hMixedFit.get(), 1.0, 1.0);
    }
    else {
        hCk = GetTH3DClone(fin, cfg.CfHistoName, "hDataCk3D");
        if(cfg.RebinFactor > 1){
            std::cout << "WARNING BuildDataCk3D: RebinFactor is only applied before same/mixed division. "
                      << "Direct CF input is left unrebinned.\n";
        }
    }

    if(cfg.NormalizeDataInNormWindow){
        const double norm = NormalizeCk3D(hCk, cfg);
        if(normFactor) *normFactor = norm;
    }
    return hCk;
}

} // namespace

TH3D* BuildDataCk3D(const FitCATS3DConfig& cfg)
{
    return BuildDataCk3DInternal(cfg, nullptr, nullptr, nullptr);
}

namespace {

FitCATS3DResult FitDataCATS3DImpl(const FitCATS3DConfig& cfg,
                                  DLM_CleverReso3D* resonanceSource)
{
    const auto tData0 = Clock3D::now();
    std::unique_ptr<TH3D> sameForPml;
    std::unique_ptr<TH3D> mixedForPml;
    double dataNormFactor = 1.0;
    std::unique_ptr<TH3D> data(BuildDataCk3DInternal(cfg, &sameForPml,
                                                     &mixedForPml,
                                                     &dataNormFactor));
    CATS3D kitty;

    FitState3D state;
    state.cfg = cfg;
    state.Profile.BuildDataSec = SecondsSince(tData0);
    state.DataNormFactor = dataNormFactor;

    if(cfg.UseResonanceSource && !resonanceSource){
        throw std::runtime_error("FitDataCATS3D: UseResonanceSource=true requires a prepared DLM_CleverReso3D source. Build it in the macro and pass it to FitDataCATS3D(cfg, source).");
    }

    state.ResonanceSource = cfg.UseResonanceSource ? resonanceSource : nullptr;

    const auto tConfig0 = Clock3D::now();
    ConfigureKittyForData(kitty, *data, cfg, state.ResonanceSource);
    state.Profile.ConfigureKittySec = SecondsSince(tConfig0);

    if(cfg.UseExternalFastKernel){
        if(!cfg.UseFastKernel){
            throw std::runtime_error("FitDataCATS3D: UseExternalFastKernel requires UseFastKernel=true.");
        }
        if(!cfg.ExternalFastKernelFileName.Length()){
            throw std::runtime_error("FitDataCATS3D: ExternalFastKernelFileName is empty.");
        }
        const auto tKernel0 = Clock3D::now();
        const bool loaded = cfg.ExternalFastKernelByQ
            ? LoadCATS3DFastKernelByQ(kitty,
                                      cfg.ExternalFastKernelFileName,
                                      cfg.ExternalFastKernelTreeName)
            : LoadCATS3DFastKernel(kitty,
                                   cfg.ExternalFastKernelFileName,
                                   cfg.ExternalFastKernelTreeName);
        state.Profile.LoadExternalKernelSec = SecondsSince(tKernel0);
        if(!loaded){
            throw std::runtime_error(Form("FitDataCATS3D: failed to load external fast kernel '%s'.",
                                          cfg.ExternalFastKernelFileName.Data()));
        }
        std::cout << "CATS3D external fast kernel loaded in "
                  << state.Profile.LoadExternalKernelSec << " s\n";
    }

    state.data = data.get();
    state.same = sameForPml.get();
    state.mixed = mixedForPml.get();
    state.kitty = &kitty;
    state.BkgNormInit = cfg.NormInit;
    state.BkgParInit[0] = cfg.NormInit;
    PrefitBackground3D(state);

    const auto tLambda0 = Clock3D::now();
    state.LambdaFeedDown.reset(LoadDlmFeedDownFraction(cfg));
    state.Profile.LoadLambdaSec = SecondsSince(tLambda0);
    gFitState3D = &state;

    if(state.LambdaFeedDown){
        std::cout << "Loaded lambda feeddown histogram: "
                  << cfg.LambdaFileName << " :: " << cfg.LambdaHistoName
                  << "\n";
        PrintLambdaFeedDownSummary(state);
    }
    std::cout << "CATS3D objective: "
              << ((cfg.UsePmlChi2 && state.same && state.mixed)
                  ? "Poisson maximum-likelihood chi2_PML"
                  : "Gaussian chi2 on normalized Ck")
              << " dataNormFactor=" << state.DataNormFactor << "\n";

    TMinuit minuit(14);
    minuit.SetFCN(FitFCN3D);
    minuit.SetPrintLevel(1);

    Double_t arglist[2] = {1.0, 0.0};
    Int_t ierflg = 0;
    minuit.mnexcm("SET ERR", arglist, 1, ierflg);

    minuit.DefineParameter(0, "Rout",  cfg.RoutInit,  cfg.MinuitRadiusStep,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(1, "Rside", cfg.RsideInit, cfg.MinuitRadiusStep,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(2, "Rlong", cfg.RlongInit, cfg.MinuitRadiusStep,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(3, "Alpha", cfg.ResonanceAlphaInit, 0.02,
                           cfg.ResonanceAlphaMin, cfg.ResonanceAlphaMax);
    const char* bkgParNames[10] = {
        "BkgP0", "BkgOut", "BkgSide", "BkgLong",
        "BkgOut2", "BkgSide2", "BkgLong2",
        "BkgOutSide", "BkgOutLong", "BkgSideLong"
    };
    const double bkgStep[10] = {
        0.01, 0.05, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2
    };
    for(int iPar=0; iPar<10; ++iPar){
        const int minuitIndex = 4 + iPar;
        const double minLimit = iPar == 0 ? 0.5 : 0.0;
        const double maxLimit = iPar == 0 ? 1.5 : 0.0;
        minuit.DefineParameter(minuitIndex, bkgParNames[iPar],
                               state.BkgParInit[iPar], bkgStep[iPar],
                               minLimit, maxLimit);
    }

    if(!cfg.UseResonanceSource || cfg.FixResonanceAlpha){
        minuit.FixParameter(3);
    }

    if(cfg.UseBackground){
        if(cfg.FixBackground){
            for(int iPar=4; iPar<14; ++iPar){
                minuit.FixParameter(iPar);
            }
        }
    }
    else {
        for(int iPar=5; iPar<14; ++iPar){
            minuit.FixParameter(iPar);
        }
        if(cfg.FixNorm){
            minuit.FixParameter(4);
        }
    }

    const auto tMinuit0 = Clock3D::now();
    if(cfg.RunMinuitSimplexFirst){
        arglist[0] = cfg.MinuitSimplexMaxCalls;
        arglist[1] = cfg.MinuitTolerance;
        std::cout << "CATS3D Minuit preconditioning: SIMPLEX calls="
                  << cfg.MinuitSimplexMaxCalls
                  << " tolerance=" << cfg.MinuitTolerance
                  << " radiusStep=" << cfg.MinuitRadiusStep << "\n";
        minuit.mnexcm("SIMPLEX", arglist, 2, ierflg);
    }
    arglist[0] = cfg.MinuitMaxCalls;
    arglist[1] = cfg.MinuitTolerance;
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
    if(cfg.RunMinuitHesse){
        Double_t hesseArgs[1] = {double(cfg.MinuitMaxCalls)};
        Int_t hesseFlg = 0;
        std::cout << "CATS3D Minuit covariance: HESSE calls="
                  << cfg.MinuitMaxCalls << "\n";
        minuit.mnexcm("HESSE", hesseArgs, 1, hesseFlg);
        std::cout << "CATS3D Minuit HESSE ierflg=" << hesseFlg << "\n";
        if(ierflg == 0) ierflg = hesseFlg;
    }
    state.Profile.MinuitSec = SecondsSince(tMinuit0);

    FitCATS3DResult result;
    double err = 0.0;
    minuit.GetParameter(0, result.Rout, err);
    minuit.GetParameter(1, result.Rside, err);
    minuit.GetParameter(2, result.Rlong, err);
    minuit.GetParameter(3, result.Alpha, err);
    for(int iPar=0; iPar<10; ++iPar){
        minuit.GetParameter(4+iPar, result.BkgPar[iPar], err);
    }
    result.Norm = result.BkgPar[0];

    double pars[14] = {result.Rout, result.Rside, result.Rlong,
                       result.Alpha,
                       result.BkgPar[0], result.BkgPar[1],
                       result.BkgPar[2], result.BkgPar[3],
                       result.BkgPar[4], result.BkgPar[5],
                       result.BkgPar[6], result.BkgPar[7],
                       result.BkgPar[8], result.BkgPar[9]};
    const auto tFinalChi20 = Clock3D::now();
    result.Chi2 = ComputeChi2ForPars(pars);
    state.Profile.FinalChi2Sec = SecondsSince(tFinalChi20);
    result.NFitBins = state.nFitBins;
    int nFreePars = 3;
    if(cfg.UseResonanceSource && !cfg.FixResonanceAlpha){
        ++nFreePars;
    }
    if(cfg.UseBackground){
        nFreePars += cfg.FixBackground ? 0 : 10;
    }
    else {
        nFreePars += cfg.FixNorm ? 0 : 1;
    }
    result.Ndf = result.NFitBins - nFreePars;
    result.FitStatus = ierflg;
    minuit.mnstat(result.MinuitAmin, result.MinuitEdm, result.MinuitErrDef,
                  result.MinuitNvpar, result.MinuitNparx,
                  result.MinuitIstat);
    if(result.FitStatus == 0 && result.MinuitIstat == 3){
        result.FitQuality = 1;
    }
    else if(result.FitStatus == 0 && result.MinuitIstat == 2){
        result.FitQuality = 2;
    }
    else {
        result.FitQuality = 3;
    }

    TH3D* finalRawModel = GetCachedModelForPars(pars);
    const auto tSave0 = Clock3D::now();
    SaveFitOutput(cfg, result, data.get(), finalRawModel);
    state.Profile.SaveOutputSec = SecondsSince(tSave0);
    gFitState3D = nullptr;

    std::cout << "CATS3D fit result: Rout=" << result.Rout
              << " Rside=" << result.Rside
              << " Rlong=" << result.Rlong
              << " Alpha=" << result.Alpha
              << " Norm=" << result.Norm
              << " "
              << ((cfg.UsePmlChi2 && state.same && state.mixed)
                  ? "chi2_PML/ndf=" : "chi2/ndf=")
              << result.Chi2 << "/" << result.Ndf
              << " = " << result.Chi2Ndf() << "\n"
              << "Minuit status: ierflg=" << result.FitStatus
              << " istat=" << result.MinuitIstat
              << " fitQuality=" << result.FitQuality
              << " amin=" << result.MinuitAmin
              << " edm=" << result.MinuitEdm
              << " errdef=" << result.MinuitErrDef
              << " nvpar=" << result.MinuitNvpar
              << " nparx=" << result.MinuitNparx << "\n"
              << "CATS3D builds/cache hits: " << state.nCatsEvals
              << "/" << state.nCacheHits << "\n";
    PrintFitProfile3D(state);

    return result;
}

} // namespace

FitCATS3DResult FitDataCATS3D(const FitCATS3DConfig& cfg)
{
    return FitDataCATS3DImpl(cfg, nullptr);
}

FitCATS3DResult FitDataCATS3D(const FitCATS3DConfig& cfg,
                              DLM_CleverReso3D* resonanceSource)
{
    return FitDataCATS3DImpl(cfg, resonanceSource);
}

void fitDataCATS3D()
{
    FitCATS3DConfig cfg;
    FitDataCATS3D(cfg);
}
