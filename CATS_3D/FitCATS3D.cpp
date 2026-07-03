#include "FitCATS3D.h"
#include "Basics3D.h"

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
};

FitState3D* gFitState3D = nullptr;
double gBkgRejectQminMeV = 100.0;
double gBkgQAxisScaleToMeV = 1000.0;

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
                         const double norm)
{
    (void)qOut;
    (void)qSide;
    (void)qLong;
    return norm;
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
    return Background3DValue(qOut, qSide, qLong, par[4]);
}

double Background3DRejectFemtoTF3(double* x, double* par)
{
    if(std::fabs(x[0]) * gBkgQAxisScaleToMeV < gBkgRejectQminMeV ||
       std::fabs(x[1]) * gBkgQAxisScaleToMeV < gBkgRejectQminMeV ||
       std::fabs(x[2]) * gBkgQAxisScaleToMeV < gBkgRejectQminMeV){
        TF1::RejectPoint();
    }
    return Background3DValue(x[0], x[1], x[2],
                             par[0]);
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

void NormalizeCk3D(TH3D* h, const FitCATS3DConfig& cfg)
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
        return;
    }

    const double norm = sumWC / sumW;
    if(norm > 0.0 && std::isfinite(norm)){
        h->Scale(1.0 / norm);
    }
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

void PrefitBackground3D(FitState3D& state)
{
    const FitCATS3DConfig& cfg = state.cfg;
    if(!cfg.UseBackground || !cfg.PrefitBackground || !state.data) return;

    TH3D* data = state.data;
    gBkgRejectQminMeV = cfg.BackgroundNonFemtoQminMeV;
    gBkgQAxisScaleToMeV = cfg.qAxisScaleToMeV;

    TF3 baseline("fCATS3D_background_prefit",
                 Background3DRejectFemtoTF3,
                 data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(),
                 data->GetYaxis()->GetXmin(), data->GetYaxis()->GetXmax(),
                 data->GetZaxis()->GetXmin(), data->GetZaxis()->GetXmax(),
                 1);
    baseline.SetParName(0, "BkgNorm");
    baseline.SetParameter(0, state.BkgNormInit);
    baseline.SetParLimits(0, 0.5, 1.5);

    const int fitStatus = data->Fit(&baseline, "RQ0");
    state.BkgNormInit = baseline.GetParameter(0);

    std::cout << "CATS3D background prefit: "
              << "Norm=" << state.BkgNormInit
              << " status=" << fitStatus
              << " region |q_out,side,long| > "
              << cfg.BackgroundNonFemtoQminMeV << " MeV\n";
}

void ConfigureKittyForData(CATS3D& kitty,
                           const TH3D& data,
                           const FitCATS3DConfig& cfg,
                           DLM_CleverReso3D* resonanceSource)
{
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

    const auto tLoop0 = Clock3D::now();
    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                if(!IsInsideFitBoxMeV(&data, ix, iy, iz, st.cfg)) continue;

                const double y = data.GetBinContent(ix, iy, iz);
                const double e = data.GetBinError(ix, iy, iz);
                if(!std::isfinite(y) || !std::isfinite(e) || e <= 0.0) continue;

                double background = par[4];
                if(st.cfg.UseBackground){
                    background = Background3DValue(data, ix, iy, iz, par);
                }
                if(!std::isfinite(background) || background <= 0.0) continue;

                const double model = background * catsModel->GetBinContent(ix, iy, iz);
                const double pull = (y - model) / e;
                chi2 += pull*pull;
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

    const double pars[5] = {result.Rout, result.Rside, result.Rlong,
                            result.Alpha, result.Norm};

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
    WriteFemtoSliceProjections(model.get(), "hCATS3D_best", cfg.sliceCutAbs);
    WriteFemtoSliceProjections(ratio.get(), "hDataOverCATS3D", cfg.sliceCutAbs);

    TNtuple summary("fit_summary", "fit_summary",
                    "Rout:Rside:Rlong:Alpha:Norm:Chi2:Ndf:Chi2Ndf:NFitBins:FitStatus");
    summary.Fill(result.Rout, result.Rside, result.Rlong, result.Alpha, result.Norm,
                 result.Chi2, result.Ndf, result.Chi2Ndf(),
                 result.NFitBins, result.FitStatus);
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

TH3D* BuildDataCk3D(const FitCATS3DConfig& cfg)
{
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
        NormalizeCk3D(hCk, cfg);
    }
    return hCk;
}

namespace {

FitCATS3DResult FitDataCATS3DImpl(const FitCATS3DConfig& cfg,
                                  DLM_CleverReso3D* resonanceSource)
{
    const auto tData0 = Clock3D::now();
    std::unique_ptr<TH3D> data(BuildDataCk3D(cfg));
    CATS3D kitty;

    FitState3D state;
    state.cfg = cfg;
    state.Profile.BuildDataSec = SecondsSince(tData0);

    if(cfg.UseResonanceSource && !resonanceSource){
        throw std::runtime_error("FitDataCATS3D: UseResonanceSource=true requires a prepared DLM_CleverReso3D source. Build it in the macro and pass it to FitDataCATS3D(cfg, source).");
    }

    state.ResonanceSource = cfg.UseResonanceSource ? resonanceSource : nullptr;

    const auto tConfig0 = Clock3D::now();
    ConfigureKittyForData(kitty, *data, cfg, state.ResonanceSource);
    state.Profile.ConfigureKittySec = SecondsSince(tConfig0);

    state.data = data.get();
    state.kitty = &kitty;
    state.BkgNormInit = cfg.NormInit;
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

    TMinuit minuit(5);
    minuit.SetFCN(FitFCN3D);
    minuit.SetPrintLevel(1);

    Double_t arglist[2] = {1.0, 0.0};
    Int_t ierflg = 0;
    minuit.mnexcm("SET ERR", arglist, 1, ierflg);

    minuit.DefineParameter(0, "Rout",  cfg.RoutInit,  0.1,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(1, "Rside", cfg.RsideInit, 0.1,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(2, "Rlong", cfg.RlongInit, 0.1,
                           cfg.RadiusMinFm, cfg.RadiusMaxFm);
    minuit.DefineParameter(3, "Alpha", cfg.ResonanceAlphaInit, 0.02,
                           cfg.ResonanceAlphaMin, cfg.ResonanceAlphaMax);
    minuit.DefineParameter(4, "Norm",  state.BkgNormInit,  0.01, 0.5, 1.5);

    if(!cfg.UseResonanceSource || cfg.FixResonanceAlpha){
        minuit.FixParameter(3);
    }

    if(cfg.UseBackground){
        if(cfg.FixBackground){
            minuit.FixParameter(4);
        }
    }
    if(!cfg.UseBackground && cfg.FixNorm){
        minuit.FixParameter(4);
    }

    arglist[0] = cfg.MinuitMaxCalls;
    arglist[1] = 1.0;
    const auto tMinuit0 = Clock3D::now();
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
    state.Profile.MinuitSec = SecondsSince(tMinuit0);

    FitCATS3DResult result;
    double err = 0.0;
    minuit.GetParameter(0, result.Rout, err);
    minuit.GetParameter(1, result.Rside, err);
    minuit.GetParameter(2, result.Rlong, err);
    minuit.GetParameter(3, result.Alpha, err);
    minuit.GetParameter(4, result.Norm, err);

    double pars[5] = {result.Rout, result.Rside, result.Rlong,
                      result.Alpha, result.Norm};
    const auto tFinalChi20 = Clock3D::now();
    result.Chi2 = ComputeChi2ForPars(pars);
    state.Profile.FinalChi2Sec = SecondsSince(tFinalChi20);
    result.NFitBins = state.nFitBins;
    int nFreePars = 3;
    if(cfg.UseResonanceSource && !cfg.FixResonanceAlpha){
        ++nFreePars;
    }
    if(cfg.UseBackground){
        nFreePars += cfg.FixBackground ? 0 : 1;
    }
    else {
        nFreePars += cfg.FixNorm ? 0 : 1;
    }
    result.Ndf = result.NFitBins - nFreePars;
    result.FitStatus = ierflg;

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
