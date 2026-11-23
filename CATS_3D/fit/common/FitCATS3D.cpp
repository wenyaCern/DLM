#include "FitCATS3D.h"
#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"
#include "AdaptivePipiKernel.h"

#include "CATS3D.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"

#include "TFile.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TString.h"

#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
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
    std::shared_ptr<AdaptivePipiKernel> AdaptiveKernel;
    std::vector<double> QSmearing[3];
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
    double BkgParErrInit[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0, 0.0};
    double DataNormFactor = 1.0;
    double DataNormRelError = 0.0;
};

FitState3D* gFitState3D = nullptr;
FitState3D* gBkgPrefitState3D = nullptr;
std::shared_ptr<AdaptivePipiKernel> gSharedAdaptiveKernel;
std::string gSharedAdaptiveKernelPath;
double gBkgRejectQminMeV = 100.0;
double gBkgRejectQmaxMeV = 200.0;
double gBkgQAxisScaleToMeV = 1000.0;
bool gBkgUseSphericalRegion = false;
bool gBkgUseRadialQuadratic = false;
bool gBkgUseDiagonalQuartic = false;
double gBkgBasisMean[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0};

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
    const bool inside = std::fabs(qOut) < cfg.fitQmaxMeV &&
                        std::fabs(qSide) < cfg.fitQmaxMeV &&
                        std::fabs(qLong) < cfg.fitQmaxMeV;
    if(!inside || !cfg.FitAxisSlicesOnly) return inside;
    const int nNearAxis = (std::fabs(qOut) < cfg.FitAxisSliceCutMeV ? 1 : 0) +
                          (std::fabs(qSide) < cfg.FitAxisSliceCutMeV ? 1 : 0) +
                          (std::fabs(qLong) < cfg.FitAxisSliceCutMeV ? 1 : 0);
    return nNearAxis >= 2;
}

bool IsInsideNormalizationRegionMeV(const TH3D* h,
                                    const int ix,
                                    const int iy,
                                    const int iz,
                                    const FitCATS3DConfig& cfg)
{
    const double qOut = std::fabs(h->GetXaxis()->GetBinCenter(ix) * cfg.qAxisScaleToMeV);
    const double qSide = std::fabs(h->GetYaxis()->GetBinCenter(iy) * cfg.qAxisScaleToMeV);
    const double qLong = std::fabs(h->GetZaxis()->GetBinCenter(iz) * cfg.qAxisScaleToMeV);

    if(cfg.UseWedgeNormalizationRegion){
        const double qMag = std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);
        return qMag >= cfg.normQminMeV && qMag < cfg.normQmaxMeV &&
               qSide > qOut + cfg.NormWedgeSideOffsetMeV;
    }
    if(cfg.UseSphericalNormalizationRegion){
        const double qMag = std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);
        return qMag >= cfg.normQminMeV && qMag <= cfg.normQmaxMeV;
    }

    // Use the full Cartesian box after removing the central femtoscopic cube.
    // Thus at least one component must lie outside |q_i| < normQminMeV.
    const bool insideOuterBox = qOut <= cfg.normQmaxMeV &&
                                qSide <= cfg.normQmaxMeV &&
                                qLong <= cfg.normQmaxMeV;
    const bool outsideFemtoCube = qOut >= cfg.normQminMeV ||
                                  qSide >= cfg.normQminMeV ||
                                  qLong >= cfg.normQminMeV;
    return insideOuterBox && outsideFemtoCube;
}

bool IsInsideBackgroundRegionMeV(const TH3D* h,
                                 const int ix,
                                 const int iy,
                                 const int iz,
                                 const FitCATS3DConfig& cfg)
{
    const double qo = std::fabs(h->GetXaxis()->GetBinCenter(ix)*cfg.qAxisScaleToMeV);
    const double qs = std::fabs(h->GetYaxis()->GetBinCenter(iy)*cfg.qAxisScaleToMeV);
    const double ql = std::fabs(h->GetZaxis()->GetBinCenter(iz)*cfg.qAxisScaleToMeV);
    if(cfg.UseHybridBackgroundRegion){
        const double slice = cfg.sliceCutAbs*cfg.qAxisScaleToMeV;
        const double qmin = cfg.BackgroundNonFemtoQminMeV;
        const double qmax = cfg.BackgroundNonFemtoQmaxMeV;
        const bool corner = qo >= qmin && qs >= qmin && ql >= qmin &&
                            qo <= qmax && qs <= qmax && ql <= qmax;
        const bool outTube = qo >= qmin && qo <= qmax && qs <= slice && ql <= slice;
        const bool sideTube = qs >= qmin && qs <= qmax && qo <= slice && ql <= slice;
        const bool longTube = ql >= qmin && ql <= qmax && qo <= slice && qs <= slice;
        return corner || outTube || sideTube || longTube;
    }
    if(cfg.UseSphericalBackgroundRegion){
        const double qmag = std::sqrt(qo*qo + qs*qs + ql*ql);
        return qmag >= cfg.BackgroundNonFemtoQminMeV &&
               qmag <= cfg.BackgroundNonFemtoQmaxMeV;
    }
    const bool insideOuter = qo <= cfg.BackgroundNonFemtoQmaxMeV &&
                             qs <= cfg.BackgroundNonFemtoQmaxMeV &&
                             ql <= cfg.BackgroundNonFemtoQmaxMeV;
    const bool outsideInner = qo >= cfg.BackgroundNonFemtoQminMeV ||
                              qs >= cfg.BackgroundNonFemtoQminMeV ||
                              ql >= cfg.BackgroundNonFemtoQminMeV;
    return insideOuter && outsideInner;
}

double Background3DValue(const double qOut,
                         const double qSide,
                         const double qLong,
                         const double* par)
{
    if(gBkgUseDiagonalQuartic){
        const double basis[6] = {
            qOut*qOut, qSide*qSide, qLong*qLong,
            qOut*qOut*qOut*qOut,
            qSide*qSide*qSide*qSide,
            qLong*qLong*qLong*qLong
        };
        double value = par[0];
        for(int i=0; i<3; ++i) value += par[i+4]*(basis[i]-gBkgBasisMean[i+3]);
        for(int i=0; i<3; ++i) value += par[i+7]*(basis[i+3]-gBkgBasisMean[i+6]);
        return value;
    }
    if(gBkgUseRadialQuadratic){
        const double qMag2 = qOut*qOut + qSide*qSide + qLong*qLong;
        const double meanQMag2 = gBkgBasisMean[3] + gBkgBasisMean[4] +
                                 gBkgBasisMean[5];
        return par[0] + par[4] * (qMag2 - meanQMag2);
    }
    const double basis[9] = {
        qOut, qSide, qLong,
        qOut*qOut, qSide*qSide, qLong*qLong,
        qOut*qSide, qOut*qLong, qSide*qLong
    };
    double value = par[0];
    for(int i=0; i<9; ++i){
        value += par[i+1] * (basis[i] - gBkgBasisMean[i]);
    }
    return value;
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

double FullBackground3DValue(const TH3D& h,
                             const int ix, const int iy, const int iz,
                             const double* par, const FitCATS3DConfig& cfg)
{
    const double smooth = Background3DValue(h,ix,iy,iz,par);
    if(!cfg.UseDiagonalBandBackground) return smooth;
    const double qo=h.GetXaxis()->GetBinCenter(ix);
    const double qs=h.GetYaxis()->GetBinCenter(iy);
    const double ql=h.GetZaxis()->GetBinCenter(iz);
    const double perpendicular=(qs-cfg.BandSlopeSidePerOut*qo)/cfg.BandWidthPerpGeV;
    const double longitudinal=ql/cfg.BandWidthLongGeV;
    const double depletion=par[14]*std::exp(-0.5*(perpendicular*perpendicular+
                                                  longitudinal*longitudinal));
    return smooth*(1.0-depletion);
}

double Background3DRejectFemtoTF3(double* x, double* par)
{
    const double qOut = std::fabs(x[0]) * gBkgQAxisScaleToMeV;
    const double qSide = std::fabs(x[1]) * gBkgQAxisScaleToMeV;
    const double qLong = std::fabs(x[2]) * gBkgQAxisScaleToMeV;

    if(gBkgUseSphericalRegion){
        const double qMag = std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);
        if(qMag < gBkgRejectQminMeV || qMag > gBkgRejectQmaxMeV){
            TF1::RejectPoint();
        }
        return Background3DValue(x[0], x[1], x[2], par);
    }

    const bool insideOuterBox = qOut <= gBkgRejectQmaxMeV &&
                                qSide <= gBkgRejectQmaxMeV &&
                                qLong <= gBkgRejectQmaxMeV;
    const bool outsideFemtoCube = qOut >= gBkgRejectQminMeV ||
                                  qSide >= gBkgRejectQminMeV ||
                                  qLong >= gBkgRejectQminMeV;

    if(!insideOuterBox || !outsideFemtoCube){
        TF1::RejectPoint();
    }
    return Background3DValue(x[0], x[1], x[2], par);
}

double ComputeBackgroundPml(const double* par)
{
    FitState3D& st = *gBkgPrefitState3D;
    double pml = 0.0;
    for(int ix=1; ix<=st.same->GetNbinsX(); ++ix){
        for(int iy=1; iy<=st.same->GetNbinsY(); ++iy){
            for(int iz=1; iz<=st.same->GetNbinsZ(); ++iz){
                if(!IsInsideBackgroundRegionMeV(st.same, ix, iy, iz, st.cfg)) continue;
                const double same = st.same->GetBinContent(ix, iy, iz);
                const double mixed = st.mixed->GetBinContent(ix, iy, iz);
                const double total = same + mixed;
                if(!std::isfinite(same) || !std::isfinite(mixed) ||
                   same < 0.0 || mixed < 0.0 || total <= 0.0) continue;
                const double qo = st.same->GetXaxis()->GetBinCenter(ix);
                const double qs = st.same->GetYaxis()->GetBinCenter(iy);
                const double ql = st.same->GetZaxis()->GetBinCenter(iz);
                const double background = Background3DValue(qo, qs, ql, par);
                const double rawRatioModel = st.DataNormFactor * background;
                if(!(rawRatioModel > 0.0) || !std::isfinite(rawRatioModel)) continue;
                double logLike = 0.0;
                if(same > 0.0){
                    const double arg = rawRatioModel*total /
                                       (same*(rawRatioModel+1.0));
                    if(!(arg > 0.0) || !std::isfinite(arg)) continue;
                    logLike += same*std::log(arg);
                }
                if(mixed > 0.0){
                    const double arg = total /
                                       (mixed*(rawRatioModel+1.0));
                    if(!(arg > 0.0) || !std::isfinite(arg)) continue;
                    logLike += mixed*std::log(arg);
                }
                pml += -2.0*logLike;
            }
        }
    }
    return pml;
}

void BackgroundPmlFCN(Int_t&, Double_t*, Double_t& f, Double_t* par, Int_t)
{
    f = ComputeBackgroundPml(par);
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

    if(cfg.LambdaMtBinIndex > 0){
        TH2* h2 = dynamic_cast<TH2*>(hIn);
        if(!h2 || cfg.LambdaMtBinIndex > h2->GetNbinsY()){
            throw std::runtime_error(Form(
                "LoadLambdaHistogram: '%s' is not a TH2 with mT row %d",
                cfg.LambdaHistoName.Data(), cfg.LambdaMtBinIndex));
        }
        std::unique_ptr<TH1D> projection(h2->ProjectionX(
            "hLambdaSelectedMt", cfg.LambdaMtBinIndex,
            cfg.LambdaMtBinIndex, "e"));
        projection->SetDirectory(nullptr);
        return MakeDlmFeedDownFractionFromLambdaTH1(projection.get(), cfg);
    }
    if(dynamic_cast<TH2*>(hIn)){
        throw std::runtime_error("LoadLambdaHistogram: TH2 requires LambdaMtBinIndex");
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
    double sumC = 0.0;
    unsigned long nBins = 0;

    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                if(!IsInsideNormalizationRegionMeV(h, ix, iy, iz, cfg)) continue;

                const double c = h->GetBinContent(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;

                sumC += c;
                ++nBins;
            }
        }
    }

    if(nBins == 0){
        std::cout << "WARNING NormalizeCk3D: no bins in normalization window; data is left unscaled.\n";
        return 1.0;
    }

    const double norm = sumC / static_cast<double>(nBins);
    if(norm > 0.0 && std::isfinite(norm)){
        h->Scale(1.0 / norm);
        return norm;
    }
    return 1.0;
}

double SameMixedExposureNorm3D(const TH3D& same,
                               const TH3D& mixed,
                               const FitCATS3DConfig& cfg,
                               double* relativeError)
{
    if(relativeError) *relativeError = 0.0;
    double sumSame = 0.0;
    double sumMixed = 0.0;
    for(int ix=1; ix<=same.GetNbinsX(); ++ix){
        for(int iy=1; iy<=same.GetNbinsY(); ++iy){
            for(int iz=1; iz<=same.GetNbinsZ(); ++iz){
                if(!IsInsideNormalizationRegionMeV(&same, ix, iy, iz, cfg)) continue;
                const double s = same.GetBinContent(ix, iy, iz);
                const double m = mixed.GetBinContent(ix, iy, iz);
                if(!std::isfinite(s) || !std::isfinite(m) || s < 0.0 || m < 0.0) continue;
                sumSame += s;
                sumMixed += m;
            }
        }
    }
    if(sumSame <= 0.0 || sumMixed <= 0.0){
        std::cout << "WARNING SameMixedExposureNorm3D: empty normalization region; using 1.\n";
        return 1.0;
    }
    if(relativeError){
        *relativeError = std::sqrt(1.0/sumSame + 1.0/sumMixed);
    }
    return sumSame / sumMixed;
}

double MeanInNormWindow3D(const TH3D& h,
                          const TH3D* weightRef,
                          const FitCATS3DConfig& cfg)
{
    (void)weightRef;
    double sumC = 0.0;
    unsigned long nBins = 0;

    for(int ix=1; ix<=h.GetNbinsX(); ++ix){
        for(int iy=1; iy<=h.GetNbinsY(); ++iy){
            for(int iz=1; iz<=h.GetNbinsZ(); ++iz){
                if(!IsInsideNormalizationRegionMeV(&h, ix, iy, iz, cfg)) continue;

                const double c = h.GetBinContent(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;

                sumC += c;
                ++nBins;
            }
        }
    }

    if(nBins == 0) return 0.0;
    return sumC / static_cast<double>(nBins);
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

void PrefitBackground3D(FitState3D& state)
{
    const FitCATS3DConfig& cfg = state.cfg;
    if(!cfg.UseBackground || !cfg.PrefitBackground || !state.data) return;

    TH3D* data = state.data;
    gBkgRejectQminMeV = cfg.BackgroundNonFemtoQminMeV;
    gBkgRejectQmaxMeV = cfg.BackgroundNonFemtoQmaxMeV;
    gBkgQAxisScaleToMeV = cfg.qAxisScaleToMeV;
    gBkgUseSphericalRegion = cfg.UseSphericalBackgroundRegion;
    gBkgUseRadialQuadratic = cfg.UseRadialQuadraticBackground;
    gBkgUseDiagonalQuartic = cfg.UseDiagonalQuarticBackground;

    for(double& mean : gBkgBasisMean) mean = 0.0;
    unsigned long basisBins = 0;
    for(int ix=1; ix<=data->GetNbinsX(); ++ix){
        for(int iy=1; iy<=data->GetNbinsY(); ++iy){
            for(int iz=1; iz<=data->GetNbinsZ(); ++iz){
                const double qoMeV = std::fabs(data->GetXaxis()->GetBinCenter(ix) *
                                               cfg.qAxisScaleToMeV);
                const double qsMeV = std::fabs(data->GetYaxis()->GetBinCenter(iy) *
                                               cfg.qAxisScaleToMeV);
                const double qlMeV = std::fabs(data->GetZaxis()->GetBinCenter(iz) *
                                               cfg.qAxisScaleToMeV);
                if(cfg.UseHybridBackgroundRegion){
                    if(!IsInsideBackgroundRegionMeV(data, ix, iy, iz, cfg)) continue;
                }
                else if(cfg.UseSphericalBackgroundRegion){
                    const double qMagMeV = std::sqrt(
                        qoMeV*qoMeV + qsMeV*qsMeV + qlMeV*qlMeV);
                    if(qMagMeV < cfg.BackgroundNonFemtoQminMeV ||
                       qMagMeV > cfg.BackgroundNonFemtoQmaxMeV) continue;
                }
                else {
                    const bool insideBaselineOuterBox =
                        qoMeV <= cfg.BackgroundNonFemtoQmaxMeV &&
                        qsMeV <= cfg.BackgroundNonFemtoQmaxMeV &&
                        qlMeV <= cfg.BackgroundNonFemtoQmaxMeV;
                    const bool outsideBaselineInnerCube =
                        qoMeV >= cfg.BackgroundNonFemtoQminMeV ||
                        qsMeV >= cfg.BackgroundNonFemtoQminMeV ||
                        qlMeV >= cfg.BackgroundNonFemtoQminMeV;
                    if(!insideBaselineOuterBox || !outsideBaselineInnerCube) continue;
                }
                const double c = data->GetBinContent(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;
                const double qo = data->GetXaxis()->GetBinCenter(ix);
                const double qs = data->GetYaxis()->GetBinCenter(iy);
                const double ql = data->GetZaxis()->GetBinCenter(iz);
                const double basis[9] = {
                    qo, qs, ql, qo*qo, qs*qs, ql*ql,
                    cfg.UseDiagonalQuarticBackground ? qo*qo*qo*qo : qo*qs,
                    cfg.UseDiagonalQuarticBackground ? qs*qs*qs*qs : qo*ql,
                    cfg.UseDiagonalQuarticBackground ? ql*ql*ql*ql : qs*ql
                };
                for(int i=0; i<9; ++i) gBkgBasisMean[i] += basis[i];
                ++basisBins;
            }
        }
    }
    if(basisBins > 0){
        for(double& mean : gBkgBasisMean){
            mean /= static_cast<double>(basisBins);
        }
    }

    if(cfg.LoadBackgroundShape){
        TFile shapeFile(cfg.BackgroundShapeFileName, "READ");
        TNtuple* summary = shapeFile.IsZombie() ? nullptr
                           : dynamic_cast<TNtuple*>(shapeFile.Get("fit_summary"));
        if(!summary || summary->GetEntries() < 1){
            throw std::runtime_error(Form(
                "PrefitBackground3D: cannot load baseline summary from '%s'.",
                cfg.BackgroundShapeFileName.Data()));
        }
        const char* branchNames[10] = {
            "BkgP0", "BkgOut", "BkgSide", "BkgLong",
            "BkgOut2", "BkgSide2", "BkgLong2",
            "BkgOutSide", "BkgOutLong", "BkgSideLong"
        };
        Float_t values[10] = {};
        for(int i=0; i<10; ++i){
            if(!summary->GetBranch(branchNames[i])){
                throw std::runtime_error(Form(
                    "PrefitBackground3D: missing branch '%s' in '%s'; rerun default fit.",
                    branchNames[i], cfg.BackgroundShapeFileName.Data()));
            }
            summary->SetBranchAddress(branchNames[i], &values[i]);
        }
        summary->GetEntry(0);
        for(int i=0; i<10; ++i) state.BkgParInit[i] = values[i];
        state.BkgNormInit = state.BkgParInit[0];
        std::cout << "CATS3D loaded default baseline shape from "
                  << cfg.BackgroundShapeFileName
                  << "; p0 initial=" << state.BkgParInit[0]
                  << " and shape parameters fixed in subSamp fit\n";
        return;
    }

    if(cfg.UsePmlBackgroundPrefit && state.same && state.mixed){
        TMinuit baselineMinuit(10);
        baselineMinuit.SetFCN(BackgroundPmlFCN);
        baselineMinuit.SetPrintLevel(1);
        gBkgPrefitState3D = &state;

        const char* parNames[10] = {
            "BkgP0", "BkgOut", "BkgSide", "BkgLong",
            "BkgOut2", "BkgSide2", "BkgLong2",
            "BkgOutSide", "BkgOutLong", "BkgSideLong"
        };
        const double steps[10] = {
            0.001, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05
        };
        for(int iPar=0; iPar<10; ++iPar){
            baselineMinuit.DefineParameter(iPar, parNames[iPar],
                                           state.BkgParInit[iPar], steps[iPar],
                                           iPar == 0 ? 0.5 : 0.0,
                                           iPar == 0 ? 1.5 : 0.0);
        }
        if(cfg.UseRadialQuadraticBackground){
            for(int iPar=1; iPar<10; ++iPar){
                if(iPar != 4) baselineMinuit.FixParameter(iPar);
            }
        }
        else if(cfg.UseDiagonalQuarticBackground){
            baselineMinuit.FixParameter(1);
            baselineMinuit.FixParameter(2);
            baselineMinuit.FixParameter(3);
        }
        else if(cfg.UseDiagonalQuadraticBackground){
            baselineMinuit.FixParameter(1);
            baselineMinuit.FixParameter(2);
            baselineMinuit.FixParameter(3);
            baselineMinuit.FixParameter(7);
            baselineMinuit.FixParameter(8);
            baselineMinuit.FixParameter(9);
        }

        Double_t args[2] = {50000.0, 0.01};
        Int_t status = 0;
        baselineMinuit.mnexcm("MIGRAD", args, 2, status);
        for(int iPar=0; iPar<10; ++iPar){
            baselineMinuit.GetParameter(iPar, state.BkgParInit[iPar],
                                        state.BkgParErrInit[iPar]);
        }
        state.BkgNormInit = state.BkgParInit[0];
        const double pml = ComputeBackgroundPml(state.BkgParInit);
        gBkgPrefitState3D = nullptr;
        std::cout << "CATS3D background prefit: raw SAME/MIXED PML"
                  << " p0=" << state.BkgParInit[0]
                  << " out2=" << state.BkgParInit[4]
                  << " side2=" << state.BkgParInit[5]
                  << " long2=" << state.BkgParInit[6]
                  << (cfg.UseDiagonalQuarticBackground ?
                      Form(" out4=%g side4=%g long4=%g",
                           state.BkgParInit[7], state.BkgParInit[8],
                           state.BkgParInit[9]) : "")
                  << " PML=" << pml
                  << " status=" << status
                  << " baseline region: "
                  << (cfg.UseHybridBackgroundRegion ? "hybrid sideband " :
                      (cfg.UseSphericalBackgroundRegion ? "spherical shell " :
                       "Cartesian outer region "))
                  << cfg.BackgroundNonFemtoQminMeV << "-"
                  << cfg.BackgroundNonFemtoQmaxMeV << " MeV\n";
        return;
    }

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

    if(cfg.UseRadialQuadraticBackground){
        // Keep only p0 and one common coefficient multiplying q_mag^2.
        // Parameter 4 is reused as the radial quadratic coefficient; all
        // direction-dependent and cross terms are disabled.
        for(int iPar=1; iPar<10; ++iPar){
            if(iPar != 4) baseline.FixParameter(iPar, 0.0);
        }
        baseline.SetParName(4, "BkgQmag2");
    }
    else if(cfg.UseDiagonalQuarticBackground){
        baseline.FixParameter(1, 0.0);
        baseline.FixParameter(2, 0.0);
        baseline.FixParameter(3, 0.0);
        baseline.SetParName(7, "BkgOut4");
        baseline.SetParName(8, "BkgSide4");
        baseline.SetParName(9, "BkgLong4");
    }
    else if(cfg.UseDiagonalQuadraticBackground){
        // Direction-dependent even baseline: p0 + p_out*q_out^2
        // + p_side*q_side^2 + p_long*q_long^2.  Linear and mixed
        // terms are excluded so the baseline respects q_i -> -q_i.
        baseline.FixParameter(1, 0.0);
        baseline.FixParameter(2, 0.0);
        baseline.FixParameter(3, 0.0);
        baseline.FixParameter(7, 0.0);
        baseline.FixParameter(8, 0.0);
        baseline.FixParameter(9, 0.0);
    }

    // The centered basis guarantees <B>_norm = p0.  Full data determine all
    // ten baseline parameters in the sideband.
    const int fitStatus = data->Fit(&baseline, "RQ0");
    for(int iPar=0; iPar<10; ++iPar){
        state.BkgParInit[iPar] = baseline.GetParameter(iPar);
        state.BkgParErrInit[iPar] = baseline.GetParError(iPar);
    }
    state.BkgNormInit = state.BkgParInit[0];

    std::cout << "CATS3D background prefit: "
              << (cfg.UseRadialQuadraticBackground ? "radial-pol2 p0=" :
                  (cfg.UseDiagonalQuadraticBackground ? "directional-even-pol2 p0=" :
                   "pol2 p0="))
              << state.BkgParInit[0]
              << " p1=" << state.BkgParInit[1]
              << " p2=" << state.BkgParInit[2]
              << " p3=" << state.BkgParInit[3]
              << " arithmetic mean=p0"
              << " status=" << fitStatus
              << " baseline region: "
              << (cfg.UseSphericalBackgroundRegion ? "spherical shell " : "Cartesian outer region ")
              << cfg.BackgroundNonFemtoQminMeV << "-"
              << cfg.BackgroundNonFemtoQmaxMeV << " MeV\n";
}

void ConfigureKittyForData(CATS3D& kitty,
                           const TH3D& data,
                           const FitCATS3DConfig& cfg,
                           DLM_CleverReso3D* resonanceSource)
{
    if(cfg.UseFineModelGrid){
        // The fine kernel has the same number of bins as the data but spans
        // only the central half of each momentum axis.
        kitty.SetMomBinsOut(data.GetNbinsX(), -cfg.fitQmaxMeV, cfg.fitQmaxMeV);
        kitty.SetMomBinsSide(data.GetNbinsY(), -cfg.fitQmaxMeV, cfg.fitQmaxMeV);
        kitty.SetMomBinsLong(data.GetNbinsZ(), -cfg.fitQmaxMeV, cfg.fitQmaxMeV);
    }
    else {
        kitty.SetMomBinsOut(data.GetNbinsX(),
                            data.GetXaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                            data.GetXaxis()->GetXmax()*cfg.qAxisScaleToMeV);
        kitty.SetMomBinsSide(data.GetNbinsY(),
                             data.GetYaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                             data.GetYaxis()->GetXmax()*cfg.qAxisScaleToMeV);
        kitty.SetMomBinsLong(data.GetNbinsZ(),
                             data.GetZaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                             data.GetZaxis()->GetXmax()*cfg.qAxisScaleToMeV);
    }

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

    if(st.cfg.UseFineModelGrid){
        const double qMax = st.cfg.fitQmaxMeV;
        const int nFine = int(kitty.GetNumMomBinsOut());
        const double fineStep = 2.0*qMax/double(nFine);
        const TAxis* axes[3] = {data.GetXaxis(), data.GetYaxis(), data.GetZaxis()};
        for(const TAxis* axis : axes){
            if(std::fabs(axis->GetBinWidth(1)*st.cfg.qAxisScaleToMeV - 2.0*fineStep) > 1e-6)
                throw std::runtime_error("Fine model requires data bins twice the fine q width");
        }
        for(int ix=1; ix<=data.GetNbinsX(); ++ix){
            const double loX=data.GetXaxis()->GetBinLowEdge(ix)*st.cfg.qAxisScaleToMeV;
            for(int iy=1; iy<=data.GetNbinsY(); ++iy){
                const double loY=data.GetYaxis()->GetBinLowEdge(iy)*st.cfg.qAxisScaleToMeV;
                for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                    const double loZ=data.GetZaxis()->GetBinLowEdge(iz)*st.cfg.qAxisScaleToMeV;
                    const double lo[3]={loX,loY,loZ};
                    int start[3];
                    bool inFine=true;
                    for(int a=0;a<3;++a){
                        const double index=(lo[a]+qMax)/fineStep;
                        start[a]=int(std::lround(index));
                        if(std::fabs(index-double(start[a]))>1e-5 || start[a]<0 || start[a]+1>=nFine)
                            inFine=false;
                    }
                    if(!inFine){
                        model->SetBinContent(ix,iy,iz,1.0);
                        continue;
                    }
                    double total=0.0;
                    for(int dx=0;dx<2;++dx)
                        for(int dy=0;dy<2;++dy)
                            for(int dz=0;dz<2;++dz){
                                const int fx=start[0]+dx,fy=start[1]+dy,fz=start[2]+dz;
                                const double qx=-qMax+(double(fx)+0.5)*fineStep;
                                const double qy=-qMax+(double(fy)+0.5)*fineStep;
                                const double qz=-qMax+(double(fz)+0.5)*fineStep;
                                const double kstar=st.cfg.relativeMomentumFactor*
                                    std::sqrt(qx*qx+qy*qy+qz*qz);
                                const double lambda=LambdaMainAtKstar(st,kstar);
                                total+=1.0+lambda*(kitty.GetCorrFun(fx,fy,fz)-1.0);
                            }
                    model->SetBinContent(ix,iy,iz,total/8.0);
                }
            }
        }
        NormalizeModelCk3D(model,data,st.cfg);
        return model;
    }

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
    return model;
}

double AxisOverlapFraction(const TAxis& fine, int fineBin,
                           const TAxis& coarse, int coarseBin)
{
    const double low=std::max(fine.GetBinLowEdge(fineBin),coarse.GetBinLowEdge(coarseBin));
    const double high=std::min(fine.GetBinUpEdge(fineBin),coarse.GetBinUpEdge(coarseBin));
    return std::max(0.0,high-low)/fine.GetBinWidth(fineBin);
}

void LoadFactorizedQSmearing(FitState3D& st, const TH3D& data)
{
    if(!st.cfg.UseFactorizedQSmearing) return;
    if(!st.cfg.UseAdaptiveKernel)
        throw std::runtime_error("Factorized q smearing currently requires the adaptive kernel");
    TFile file(st.cfg.QSmearingFileName,"READ");
    if(file.IsZombie()) throw std::runtime_error("Cannot open q smearing ROOT file");
    const char* names[3]={"out","side","long"};
    const TAxis* dataAxes[3]={data.GetXaxis(),data.GetYaxis(),data.GetZaxis()};
    for(int axis=0;axis<3;++axis){
        const TString path=st.cfg.QSmearingHistogramBase+
            Form("hTrueQ%sVsQ%s",names[axis],names[axis]);
        auto* h=dynamic_cast<TH2*>(file.Get(path));
        if(!h) throw std::runtime_error(Form("Missing q response: %s",path.Data()));
        const int n=dataAxes[axis]->GetNbins();
        const int nf=h->GetNbinsX();
        if(n!=data.GetNbinsX() || n!=data.GetNbinsY() || n!=data.GetNbinsZ() ||
           nf!=h->GetNbinsY())
            throw std::runtime_error("Unexpected q response or data binning");
        const TAxis& truth=*h->GetXaxis();
        const TAxis& reco=*h->GetYaxis();
        const TAxis& target=*dataAxes[axis];
        if(std::fabs(truth.GetXmin()-target.GetXmin())>1e-8 ||
           std::fabs(truth.GetXmax()-target.GetXmax())>1e-8 ||
           std::fabs(reco.GetXmin()-target.GetXmin())>1e-8 ||
           std::fabs(reco.GetXmax()-target.GetXmax())>1e-8)
            throw std::runtime_error("Q response range must match fit-data range");
        std::vector<double> joint(size_t(n)*n,0.0), outside(n,0.0);
        for(int j=1;j<=nf;++j){
            for(int b=1;b<=n;++b){
                const double wy=AxisOverlapFraction(reco,j,target,b);
                if(wy==0) continue;
                outside[b-1]+=(h->GetBinContent(0,j)+h->GetBinContent(nf+1,j))*wy;
                for(int i=1;i<=nf;++i){
                    const double count=h->GetBinContent(i,j)*wy;
                    if(count==0) continue;
                    for(int a=1;a<=n;++a){
                        const double wx=AxisOverlapFraction(truth,i,target,a);
                        if(wx>0) joint[size_t(a-1)*n+b-1]+=count*wx;
                    }
                }
            }
        }
        auto& response=st.QSmearing[axis];
        response.resize(size_t(n)*n);
        for(int b=0;b<n;++b){
            double denominator=outside[b];
            for(int a=0;a<n;++a) denominator+=joint[size_t(a)*n+b];
            if(denominator<=0) throw std::runtime_error("Empty reconstructed q response bin");
            for(int a=0;a<n;++a)
                response[size_t(a)*n+b]=joint[size_t(a)*n+b]/denominator;
        }
        std::cout << "Loaded q_" << names[axis] << " response: " << nf
                  << " MC bins -> " << n << " fit bins; outside-truth fraction at q=0: "
                  << outside[n/2]/(outside[n/2]+[&](){double sum=0;for(int a=0;a<n;++a)
                      sum+=joint[size_t(a)*n+n/2];return sum;}()) << "\n";
    }
}

void SmearAdaptiveModel(TH3D& model, const FitState3D& st)
{
    if(!st.cfg.UseFactorizedQSmearing) return;
    const int nx=model.GetNbinsX(), ny=model.GetNbinsY(), nz=model.GetNbinsZ();
    const size_t total=size_t(nx)*ny*nz;
    std::vector<double> delta(total), xStage(total), yStage(total);
    const auto index=[ny,nz](int x,int y,int z){return (size_t(x)*ny+y)*nz+z;};
    for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) for(int z=0;z<nz;++z)
        delta[index(x,y,z)]=model.GetBinContent(x+1,y+1,z+1)-1.0;
    for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) for(int z=0;z<nz;++z)
        for(int a=0;a<nx;++a)
            xStage[index(x,y,z)]+=st.QSmearing[0][size_t(a)*nx+x]*delta[index(a,y,z)];
    for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) for(int z=0;z<nz;++z)
        for(int a=0;a<ny;++a)
            yStage[index(x,y,z)]+=st.QSmearing[1][size_t(a)*ny+y]*xStage[index(x,a,z)];
    for(int x=0;x<nx;++x) for(int y=0;y<ny;++y) for(int z=0;z<nz;++z){
        double correction=0;
        for(int a=0;a<nz;++a)
            correction+=st.QSmearing[2][size_t(a)*nz+z]*yStage[index(x,y,a)];
        model.SetBinContent(x+1,y+1,z+1,1.0+correction);
    }
}

TH3D* BuildModelFromAdaptive(const TH3D& data, const TString& name,
                             double rout, double rside, double rlong)
{
    FitState3D& st=*gFitState3D;
    const AdaptivePipiKernel& kernel=*st.AdaptiveKernel;
    const std::vector<double> core=kernel.Fold(rout,rside,rlong);
    TH3D* model=dynamic_cast<TH3D*>(data.Clone(name));
    model->SetDirectory(nullptr);
    model->Reset("ICES");
    const TAxis* axes[3]={data.GetXaxis(),data.GetYaxis(),data.GetZaxis()};
    std::vector<std::vector<std::pair<unsigned,double>>> overlaps[3];
    for(unsigned a=0;a<3;++a){
        overlaps[a].resize(axes[a]->GetNbins());
        for(int i=1;i<=axes[a]->GetNbins();++i)
            overlaps[a][i-1]=kernel.AxisOverlaps(a,
                axes[a]->GetBinLowEdge(i)*st.cfg.qAxisScaleToMeV,
                axes[a]->GetBinUpEdge(i)*st.cfg.qAxisScaleToMeV);
    }
    for(int ix=1;ix<=data.GetNbinsX();++ix)
      for(int iy=1;iy<=data.GetNbinsY();++iy)
        for(int iz=1;iz<=data.GetNbinsZ();++iz){
            double correction=0;
            const double volume=axes[0]->GetBinWidth(ix)*axes[1]->GetBinWidth(iy)*
                                axes[2]->GetBinWidth(iz)*std::pow(st.cfg.qAxisScaleToMeV,3);
            for(const auto& x:overlaps[0][ix-1])
              for(const auto& y:overlaps[1][iy-1])
                for(const auto& z:overlaps[2][iz-1]){
                    const double qx=kernel.QCenter(0,x.first);
                    const double qy=kernel.QCenter(1,y.first);
                    const double qz=kernel.QCenter(2,z.first);
                    const double kstar=st.cfg.relativeMomentumFactor*std::sqrt(qx*qx+qy*qy+qz*qz);
                    correction+=x.second*y.second*z.second*
                        LambdaMainAtKstar(st,kstar)*(core[kernel.QIndex(x.first,y.first,z.first)]-1.0);
                }
            model->SetBinContent(ix,iy,iz,1.0+correction/volume);
        }
    SmearAdaptiveModel(*model,st);
    NormalizeModelCk3D(model,data,st.cfg);
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
    if(!st.AdaptiveKernel){
        const auto tSet0 = Clock3D::now();
        kitty.SetAnaSource(0, Rout, true);
        kitty.SetAnaSource(1, Rside, true);
        kitty.SetAnaSource(2, Rlong, true);
        if(st.cfg.UseResonanceSource) kitty.SetAnaSource(3, Alpha, true);
        st.Profile.SetAnaSourceSec += SecondsSince(tSet0);
        ++st.Profile.SetAnaSourceCalls;
    }

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
    if(!st.AdaptiveKernel) kitty.KillTheCat();
    const double killSec = SecondsSince(t0);
    st.Profile.KillTheCatSec += killSec;
    ++st.Profile.KillTheCatCalls;
    if(printTiming){
        std::cout << killSec << " s";
    }

    const auto tModel0 = Clock3D::now();
    std::unique_ptr<TH3D> newModel(st.AdaptiveKernel
        ? BuildModelFromAdaptive(data,Form("hCATS3D_cache_%u", st.CacheClock),Rout,Rside,Rlong)
        : BuildModelFromKitty(data, kitty,Form("hCATS3D_cache_%u", st.CacheClock)));
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
                    background = FullBackground3DValue(data, ix, iy, iz, par, st.cfg);
                }
                if(!std::isfinite(background) || background <= 0.0) continue;

                const double model = background * catsModel->GetBinContent(ix, iy, iz);
                if(!std::isfinite(model) || model <= 0.0) continue;

                if(usePml){
                    const double same = st.same->GetBinContent(ix, iy, iz);
                    const double mixed = st.mixed->GetBinContent(ix, iy, iz);
                    if(!std::isfinite(same) || !std::isfinite(mixed) ||
                       same < 0.0 || mixed < 0.0) continue;
                    const double total = same + mixed;
                    if(total <= 0.0) continue;

                    // SAME and MIXED remain the observed raw Poisson counts.
                    // Their relative exposure belongs in the expected raw
                    // SAME/MIXED ratio, not in either observed count.
                    const double rawRatioModel = st.DataNormFactor * model;
                    if(!(rawRatioModel > 0.0) ||
                       !std::isfinite(rawRatioModel)) continue;

                    double logLike = 0.0;
                    if(same > 0.0){
                        const double arg = rawRatioModel * total /
                                           (same * (rawRatioModel + 1.0));
                        if(!(arg > 0.0) || !std::isfinite(arg)) continue;
                        logLike += same * std::log(arg);
                    }
                    if(mixed > 0.0){
                        const double arg = total /
                                           (mixed * (rawRatioModel + 1.0));
                        if(!(arg > 0.0) || !std::isfinite(arg)) continue;
                        logLike += mixed * std::log(arg);
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

    if(target.EqualTo("out")){
        int y1, y2, z1, z2;
        centeredBinRange(h3->GetYaxis(), sliceCutAbs, y1, y2);
        centeredBinRange(h3->GetZaxis(), sliceCutAbs, z1, z2);

        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, y1, y2, z1, z2, "e"));
        const int numBins = (y2-y1+1) * (z2-z1+1);
        if(numBins > 0) h1->Scale(1.0 / double(numBins));
    }
    else if(target.EqualTo("side")){
        int x1, x2, z1, z2;
        centeredBinRange(h3->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(h3->GetZaxis(), sliceCutAbs, z1, z2);

        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, x1, x2, z1, z2, "e"));
        const int numBins = (x2-x1+1) * (z2-z1+1);
        if(numBins > 0) h1->Scale(1.0 / double(numBins));
    }
    else if(target.EqualTo("long")){
        int x1, x2, y1, y2;
        centeredBinRange(h3->GetXaxis(), sliceCutAbs, x1, x2);
        centeredBinRange(h3->GetYaxis(), sliceCutAbs, y1, y2);

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

void WriteMixedWeightedFemtoSliceProjections(TH3D* h3,
                                               const TH3D* mixed,
                                               const TString& baseName,
                                               const double sliceCutAbs)
{
    if(!h3 || !mixed) return;

    std::unique_ptr<TH3D> weighted(dynamic_cast<TH3D*>(
        h3->Clone(baseName + "_weighted3D_tmp")));
    weighted->SetDirectory(nullptr);
    for(int ix=1; ix<=h3->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h3->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h3->GetNbinsZ(); ++iz){
                const double w = mixed->GetBinContent(ix, iy, iz);
                weighted->SetBinContent(ix, iy, iz,
                                        w*h3->GetBinContent(ix, iy, iz));
                weighted->SetBinError(ix, iy, iz,
                                      w*h3->GetBinError(ix, iy, iz));
            }
        }
    }

    const TString directions[3] = {"out", "side", "long"};
    for(const TString& direction : directions){
        std::unique_ptr<TH1D> numerator(ProjectFemtoSlice1D(
            weighted.get(), direction, baseName + "_weighted_num_" + direction,
            sliceCutAbs));
        std::unique_ptr<TH1D> denominator(ProjectFemtoSlice1D(
            const_cast<TH3D*>(mixed), direction,
            baseName + "_weighted_den_" + direction, sliceCutAbs));
        numerator->Divide(denominator.get());
        numerator->SetName(baseName + "_mixedWeighted_proj_" + direction);
        numerator->Write("", TObject::kOverwrite);
    }
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

    const double pars[15] = {result.Rout, result.Rside, result.Rlong,
                             result.Alpha,
                             result.BkgPar[0], result.BkgPar[1],
                             result.BkgPar[2], result.BkgPar[3],
                             result.BkgPar[4], result.BkgPar[5],
                             result.BkgPar[6], result.BkgPar[7],
                             result.BkgPar[8], result.BkgPar[9],
                             result.BandAmplitude};

    for(int ix=1; ix<=model->GetNbinsX(); ++ix){
        for(int iy=1; iy<=model->GetNbinsY(); ++iy){
            for(int iz=1; iz<=model->GetNbinsZ(); ++iz){
                double bkg = result.Norm;
                if(cfg.UseBackground){
                    bkg = FullBackground3DValue(*model, ix, iy, iz, pars, cfg);
                }
                background->SetBinContent(ix, iy, iz, bkg);
                if(cfg.OutputBaselineOutsideFitBox &&
                   !IsInsideFitBoxMeV(model.get(),ix,iy,iz,cfg)){
                    modelShape->SetBinContent(ix,iy,iz,1.0);
                    modelShape->SetBinError(ix,iy,iz,0.0);
                    model->SetBinContent(ix,iy,iz,bkg);
                    model->SetBinError(ix,iy,iz,0.0);
                    continue;
                }
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
    if(cfg.UseAdaptiveKernel)
        TNamed("AdaptiveKernelFile",cfg.AdaptiveKernelFileName.Data()).Write();
    data->Write("hDataCk3D");
    modelShape->Write("hCATS3D_shape");
    background->Write("hCATS3D_background");
    model->Write("hCATS3D_best");
    ratio->Write("hDataOverCATS3D");
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
    if(result.BkgParErr[0] > 0.0 && std::isfinite(result.BkgParErr[0])){
        std::unique_ptr<TH3D> bkgLow(dynamic_cast<TH3D*>(
            background->Clone("hCATS3D_background_1sigma_low")));
        std::unique_ptr<TH3D> bkgHigh(dynamic_cast<TH3D*>(
            background->Clone("hCATS3D_background_1sigma_high")));
        bkgLow->SetDirectory(nullptr);
        bkgHigh->SetDirectory(nullptr);
        for(int ix=1; ix<=background->GetNbinsX(); ++ix){
            for(int iy=1; iy<=background->GetNbinsY(); ++iy){
                for(int iz=1; iz<=background->GetNbinsZ(); ++iz){
                    bkgLow->SetBinContent(ix, iy, iz,
                        std::max(0.0, result.Norm-result.BkgParErr[0]));
                    bkgHigh->SetBinContent(ix, iy, iz,
                        result.Norm+result.BkgParErr[0]);
                }
            }
        }
        bkgLow->Write();
        bkgHigh->Write();
        WriteFemtoSliceProjections(bkgLow.get(),
            "hCATS3D_background_1sigma_low", cfg.sliceCutAbs);
        WriteFemtoSliceProjections(bkgHigh.get(),
            "hCATS3D_background_1sigma_high", cfg.sliceCutAbs);
    }
    WriteFemtoSliceProjections(model.get(), "hCATS3D_best", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(ratio.get(), "hDataOverCATS3D", cfg.sliceCutAbs);
    if(gFitState3D && gFitState3D->mixed){
        WriteMixedWeightedFemtoSliceProjections(
            data, gFitState3D->mixed, "hDataCk3D", cfg.sliceCutAbs);
        WriteMixedWeightedFemtoSliceProjections(
            modelShape.get(), gFitState3D->mixed, "hCATS3D_shape", cfg.sliceCutAbs);
        WriteMixedWeightedFemtoSliceProjections(
            background.get(), gFitState3D->mixed, "hCATS3D_background", cfg.sliceCutAbs);
        WriteMixedWeightedFemtoSliceProjections(
            model.get(), gFitState3D->mixed, "hCATS3D_best", cfg.sliceCutAbs);
    }

    TNtuple summary("fit_summary", "fit_summary",
                    "Rout:Rside:Rlong:Alpha:Norm:Chi2:Ndf:Chi2Ndf:NFitBins:FitStatus:"
                    "BkgP0:BkgOut:BkgSide:BkgLong:BkgOut2:BkgSide2:BkgLong2:"
                    "BkgOutSide:BkgOutLong:BkgSideLong:BkgP0Err:BandAmp:BandAmpErr");
    Float_t summaryValues[23] = {
        float(result.Rout), float(result.Rside), float(result.Rlong),
        float(result.Alpha), float(result.Norm), float(result.Chi2),
        float(result.Ndf), float(result.Chi2Ndf()), float(result.NFitBins),
        float(result.FitStatus),
        float(result.BkgPar[0]), float(result.BkgPar[1]),
        float(result.BkgPar[2]), float(result.BkgPar[3]),
        float(result.BkgPar[4]), float(result.BkgPar[5]),
        float(result.BkgPar[6]), float(result.BkgPar[7]),
        float(result.BkgPar[8]), float(result.BkgPar[9]),
        float(result.BkgParErr[0]),
        float(result.BandAmplitude), float(result.BandAmplitudeErr)
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

} // namespace

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
    bool normalizedFromRawCounts = false;

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
        if(cfg.NormalizeDataInNormWindow){
            double normRelError = 0.0;
            const double norm = SameMixedExposureNorm3D(*hSameFit, *hMixedFit,
                                                        cfg, &normRelError);
            hCk->Scale(1.0 / norm);
            if(normFactor) *normFactor = norm;
            normalizedFromRawCounts = true;
            std::cout << "CATS3D raw-count exposure normalization: sumSame/sumMixed="
                      << norm << " relativeError=" << normRelError << "\n";
        }
    }
    else {
        hCk = GetTH3DClone(fin, cfg.CfHistoName, "hDataCk3D");
        if(cfg.RebinFactor > 1){
            std::cout << "WARNING BuildDataCk3D: RebinFactor is only applied before same/mixed division. "
                      << "Direct CF input is left unrebinned.\n";
        }
    }

    if(cfg.NormalizeDataInNormWindow && !normalizedFromRawCounts){
        const double norm = NormalizeCk3D(hCk, cfg);
        if(normFactor) *normFactor = norm;
    }
    return hCk;
}

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
    if(sameForPml && mixedForPml){
        double ignoredNorm = SameMixedExposureNorm3D(*sameForPml, *mixedForPml,
                                                     cfg, &state.DataNormRelError);
        (void)ignoredNorm;
    }

    if(cfg.UseResonanceSource && !resonanceSource){
        throw std::runtime_error("FitDataCATS3D: UseResonanceSource=true requires a prepared DLM_CleverReso3D source. Build it in the macro and pass it to FitDataCATS3D(cfg, source).");
    }

    state.ResonanceSource = cfg.UseResonanceSource ? resonanceSource : nullptr;

    const auto tConfig0 = Clock3D::now();
    ConfigureKittyForData(kitty, *data, cfg, state.ResonanceSource);
    state.Profile.ConfigureKittySec = SecondsSince(tConfig0);

    if(cfg.UseAdaptiveKernel){
        if(cfg.UseExternalFastKernel || cfg.UseResonanceSource || !cfg.Identical || !cfg.WithCoulomb)
            throw std::runtime_error("Adaptive kernel requires pi+pi+ Coulomb+QS Gaussian source and no legacy external kernel");
        if(!cfg.AdaptiveKernelFileName.Length())
            throw std::runtime_error("AdaptiveKernelFileName is empty");
        const auto tKernel0=Clock3D::now();
        const std::string kernelPath=cfg.AdaptiveKernelFileName.Data();
        if(!gSharedAdaptiveKernel || gSharedAdaptiveKernelPath!=kernelPath){
            gSharedAdaptiveKernel=std::make_shared<AdaptivePipiKernel>(kernelPath);
            gSharedAdaptiveKernelPath=kernelPath;
        } else {
            std::cout << "Reusing adaptive kernel already loaded in this process\n";
        }
        state.AdaptiveKernel=gSharedAdaptiveKernel;
        state.Profile.LoadExternalKernelSec=SecondsSince(tKernel0);
        std::cout << "Adaptive kernel ready in " << state.Profile.LoadExternalKernelSec << " s\n";
    }

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
            throw std::runtime_error(Form(
                "FitDataCATS3D: failed to load external fast kernel '%s'.",
                cfg.ExternalFastKernelFileName.Data()));
        }
        std::cout << "CATS3D external fast kernel loaded in "
                  << state.Profile.LoadExternalKernelSec << " s\n";
    }

    state.data = data.get();
    state.same = sameForPml.get();
    state.mixed = mixedForPml.get();
    state.kitty = &kitty;
    LoadFactorizedQSmearing(state,*data);
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

    const bool usePmlObjective = cfg.UsePmlChi2 && state.same && state.mixed;
    std::cout << "CATS3D objective: "
              << (usePmlObjective
                  ? "Poisson maximum-likelihood chi2_PML"
                  : "Gaussian chi2 on normalized Ck")
              << " dataNormFactor=" << state.DataNormFactor
              << " p0=" << (cfg.FixNorm ? "fixed" : "free") << "\n";

    TMinuit minuit(15);
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
    minuit.DefineParameter(14,"BandAmp",cfg.BandAmplitudeInit,0.001,0.0,0.08);
    if(!cfg.UseDiagonalBandBackground || cfg.FixBandAmplitude) minuit.FixParameter(14);

    if(!cfg.UseResonanceSource || cfg.FixResonanceAlpha){
        minuit.FixParameter(3);
    }

    if(cfg.UseBackground){
        if(cfg.FixBackground){
            // Keep the nine polynomial shape parameters fixed.  The centered
            // constant p0 is controlled separately by FixNorm.
            for(int iPar=5; iPar<14; ++iPar){
                minuit.FixParameter(iPar);
            }
            if(cfg.FixNorm) minuit.FixParameter(4);
        }
        else if(cfg.UseDiagonalQuadraticBackground &&
                cfg.FloatDiagonalQuadraticBackgroundInFinalFit){
            // Float only the physically even diagonal terms in the final fit.
            // Linear and mixed terms remain disabled.
            minuit.FixParameter(5);
            minuit.FixParameter(6);
            minuit.FixParameter(7);
            minuit.FixParameter(11);
            minuit.FixParameter(12);
            minuit.FixParameter(13);
            if(cfg.FixNorm) minuit.FixParameter(4);
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
    const auto tMinuit0 = Clock3D::now();
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
    minuit.GetParameter(0, result.Rout, result.RoutErr);
    minuit.GetParameter(1, result.Rside, result.RsideErr);
    minuit.GetParameter(2, result.Rlong, result.RlongErr);
    minuit.GetParameter(3, result.Alpha, result.AlphaErr);
    minuit.GetParameter(14, result.BandAmplitude, result.BandAmplitudeErr);
    for(int iPar=0; iPar<10; ++iPar){
        minuit.GetParameter(4+iPar, result.BkgPar[iPar], err);
        result.BkgParErr[iPar] = (iPar == 0 && !cfg.FixNorm)
                               ? err : state.BkgParErrInit[iPar];
    }
    result.Norm = result.BkgPar[0];

    double pars[15] = {result.Rout, result.Rside, result.Rlong,
                       result.Alpha,
                       result.BkgPar[0], result.BkgPar[1],
                       result.BkgPar[2], result.BkgPar[3],
                       result.BkgPar[4], result.BkgPar[5],
                       result.BkgPar[6], result.BkgPar[7],
                       result.BkgPar[8], result.BkgPar[9], result.BandAmplitude};
    const auto tFinalChi20 = Clock3D::now();
    result.Chi2 = ComputeChi2ForPars(pars);
    state.Profile.FinalChi2Sec = SecondsSince(tFinalChi20);
    result.NFitBins = state.nFitBins;
    int nFreePars = 3;
    if(cfg.UseResonanceSource && !cfg.FixResonanceAlpha){
        ++nFreePars;
    }
    if(cfg.UseDiagonalBandBackground && !cfg.FixBandAmplitude) ++nFreePars;
    if(cfg.UseBackground){
        nFreePars += cfg.FixBackground ? (cfg.FixNorm ? 0 : 1) : 10;
    }
    else {
        nFreePars += cfg.FixNorm ? 0 : 1;
    }
    result.Ndf = result.NFitBins - nFreePars;
    result.FitStatus = ierflg;

    // Only Rout, Rside, Rlong float in this diagnostic fit.  Save their full
    // covariance so that the displayed model band includes correlations.
    double radiusCovariance[9] = {};
    const bool saveRadiusBand = cfg.UseAdaptiveKernel &&
        cfg.UseFactorizedQSmearing && cfg.FixBackground && cfg.FixNorm &&
        !cfg.UseResonanceSource && !cfg.UseDiagonalBandBackground;
    if(saveRadiusBand) minuit.mnemat(radiusCovariance,3);

    TH3D* finalRawModel = GetCachedModelForPars(pars);
    const auto tSave0 = Clock3D::now();
    SaveFitOutput(cfg, result, data.get(), finalRawModel);
    if(saveRadiusBand){
        TFile output(cfg.OutputFileName,"UPDATE");
        if(output.IsZombie()) throw std::runtime_error("Cannot append radius uncertainty to fit output");
        TH2D covariance("hRadiusCovariance",";R parameter;R parameter",3,0,3,3,0,3);
        const char* labels[3]={"Rout","Rside","Rlong"};
        for(int a=0;a<3;++a){
            covariance.GetXaxis()->SetBinLabel(a+1,labels[a]);
            covariance.GetYaxis()->SetBinLabel(a+1,labels[a]);
            for(int b=0;b<3;++b)
                covariance.SetBinContent(a+1,b+1,radiusCovariance[a*3+b]);
        }
        covariance.Write();
        for(int a=0;a<3;++a){
            const double step=std::max(0.05,2.0*std::sqrt(std::max(0.0,radiusCovariance[a*3+a])));
            double shifted[15];
            std::copy(pars,pars+15,shifted);
            shifted[a]=pars[a]+step;
            std::unique_ptr<TH3D> plus(dynamic_cast<TH3D*>(
                GetCachedModelForPars(shifted)->Clone(Form("hRadiusPlus_%d",a))));
            plus->SetDirectory(nullptr);
            shifted[a]=pars[a]-step;
            std::unique_ptr<TH3D> minus(dynamic_cast<TH3D*>(
                GetCachedModelForPars(shifted)->Clone(Form("hRadiusMinus_%d",a))));
            minus->SetDirectory(nullptr);
            std::unique_ptr<TH3D> derivative(dynamic_cast<TH3D*>(
                plus->Clone(Form("hCATS3D_d%s",labels[a]))));
            derivative->SetDirectory(nullptr);
            for(int ix=1;ix<=derivative->GetNbinsX();++ix)
              for(int iy=1;iy<=derivative->GetNbinsY();++iy)
                for(int iz=1;iz<=derivative->GetNbinsZ();++iz){
                    const double bkg=cfg.UseBackground ?
                        FullBackground3DValue(*derivative,ix,iy,iz,pars,cfg) : result.Norm;
                    derivative->SetBinContent(ix,iy,iz,bkg*
                        (plus->GetBinContent(ix,iy,iz)-minus->GetBinContent(ix,iy,iz))/(2.0*step));
                }
            derivative->Write();
            if(state.mixed) WriteMixedWeightedFemtoSliceProjections(
                derivative.get(),state.mixed,derivative->GetName(),cfg.sliceCutAbs);
        }
        output.Close();
    }
    state.Profile.SaveOutputSec = SecondsSince(tSave0);
    gFitState3D = nullptr;

    std::cout << "CATS3D fit result: Rout=" << result.Rout
              << " Rside=" << result.Rside
              << " Rlong=" << result.Rlong
              << " Alpha=" << result.Alpha
              << " BandAmp=" << result.BandAmplitude
              << " Norm=" << result.Norm
              << " chi2/ndf=" << result.Chi2 << "/" << result.Ndf
              << " = " << result.Chi2Ndf() << "\n"
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
