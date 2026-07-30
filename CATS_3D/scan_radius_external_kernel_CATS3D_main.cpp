#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"
#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include "TFile.h"
#include "TF3.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

double gRejectQminMeV = 100.0;
double gQAxisScaleToMeV = 1000.0;
double gSliceCutMeV = 40.0;

TH3D* CloneTH3D(TFile& f, const TString& name, const TString& cloneName)
{
    TH3D* h = dynamic_cast<TH3D*>(f.Get(name));
    if(!h) throw std::runtime_error(Form("missing TH3D '%s'", name.Data()));
    TH3D* out = dynamic_cast<TH3D*>(h->Clone(cloneName));
    out->SetDirectory(nullptr);
    if(out->GetSumw2N() == 0) out->Sumw2();
    return out;
}

TH3D* RebinClone(TH3D* h, const int rebin, const TString& name)
{
    TH3D* out = dynamic_cast<TH3D*>(h->Clone(name));
    out->SetDirectory(nullptr);
    if(out->GetSumw2N() == 0) out->Sumw2();
    if(rebin > 1) out->Rebin3D(rebin, rebin, rebin);
    return out;
}

double KstarMeV(const TH3D* h, const int ix, const int iy, const int iz,
                const FitCATS3DConfig& cfg)
{
    const double qo = h->GetXaxis()->GetBinCenter(ix)*cfg.qAxisScaleToMeV;
    const double qs = h->GetYaxis()->GetBinCenter(iy)*cfg.qAxisScaleToMeV;
    const double ql = h->GetZaxis()->GetBinCenter(iz)*cfg.qAxisScaleToMeV;
    return cfg.relativeMomentumFactor*std::sqrt(qo*qo + qs*qs + ql*ql);
}

bool InsideFitBox(const TH3D* h, const int ix, const int iy, const int iz,
                  const FitCATS3DConfig& cfg)
{
    const double qo = h->GetXaxis()->GetBinCenter(ix)*cfg.qAxisScaleToMeV;
    const double qs = h->GetYaxis()->GetBinCenter(iy)*cfg.qAxisScaleToMeV;
    const double ql = h->GetZaxis()->GetBinCenter(iz)*cfg.qAxisScaleToMeV;
    return std::fabs(qo) < cfg.fitQmaxMeV &&
           std::fabs(qs) < cfg.fitQmaxMeV &&
           std::fabs(ql) < cfg.fitQmaxMeV;
}

double NormalizeCk(TH3D* h, const FitCATS3DConfig& cfg)
{
    double sumW = 0.0;
    double sumWC = 0.0;
    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                const double k = KstarMeV(h, ix, iy, iz, cfg);
                if(k < cfg.normQminMeV || k > cfg.normQmaxMeV) continue;
                const double c = h->GetBinContent(ix, iy, iz);
                const double e = h->GetBinError(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;
                const double w = (e > 0.0 && std::isfinite(e)) ? 1.0/(e*e) : 1.0;
                sumW += w;
                sumWC += w*c;
            }
        }
    }
    if(sumW <= 0.0) return 1.0;
    const double norm = sumWC/sumW;
    if(norm > 0.0 && std::isfinite(norm)) h->Scale(1.0/norm);
    return norm;
}

double MeanNormWindow(const TH3D& h, const TH3D& weightRef,
                      const FitCATS3DConfig& cfg)
{
    double sumW = 0.0;
    double sumWC = 0.0;
    for(int ix=1; ix<=h.GetNbinsX(); ++ix){
        for(int iy=1; iy<=h.GetNbinsY(); ++iy){
            for(int iz=1; iz<=h.GetNbinsZ(); ++iz){
                const double k = KstarMeV(&h, ix, iy, iz, cfg);
                if(k < cfg.normQminMeV || k > cfg.normQmaxMeV) continue;
                const double c = h.GetBinContent(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;
                const double e = weightRef.GetBinError(ix, iy, iz);
                const double w = (e > 0.0 && std::isfinite(e)) ? 1.0/(e*e) : 1.0;
                sumW += w;
                sumWC += w*c;
            }
        }
    }
    return sumW > 0.0 ? sumWC/sumW : 0.0;
}

double BkgValue(const double qo, const double qs, const double ql,
                const double* p)
{
    return p[0]
         + p[1]*qo + p[2]*qs + p[3]*ql
         + p[4]*qo*qo + p[5]*qs*qs + p[6]*ql*ql
         + p[7]*qo*qs + p[8]*qo*ql + p[9]*qs*ql;
}

double BkgValue(const TH3D& h, const int ix, const int iy, const int iz,
                const double* p)
{
    return BkgValue(h.GetXaxis()->GetBinCenter(ix),
                    h.GetYaxis()->GetBinCenter(iy),
                    h.GetZaxis()->GetBinCenter(iz), p);
}

double BkgRejectTF3(double* x, double* p)
{
    const double qo = std::fabs(x[0])*gQAxisScaleToMeV;
    const double qs = std::fabs(x[1])*gQAxisScaleToMeV;
    const double ql = std::fabs(x[2])*gQAxisScaleToMeV;
    const bool corner = qo >= gRejectQminMeV &&
                        qs >= gRejectQminMeV &&
                        ql >= gRejectQminMeV;
    const bool outSlice = qo >= gRejectQminMeV &&
                          qs <= gSliceCutMeV &&
                          ql <= gSliceCutMeV;
    const bool sideSlice = qs >= gRejectQminMeV &&
                           qo <= gSliceCutMeV &&
                           ql <= gSliceCutMeV;
    const bool longSlice = ql >= gRejectQminMeV &&
                           qo <= gSliceCutMeV &&
                           qs <= gSliceCutMeV;
    if(!corner && !outSlice && !sideSlice && !longSlice) TF1::RejectPoint();
    return BkgValue(x[0], x[1], x[2], p);
}

void PrefitBkg(TH3D* data, const FitCATS3DConfig& cfg, double* bkg)
{
    for(int i=0; i<10; ++i) bkg[i] = 0.0;
    bkg[0] = cfg.NormInit;
    if(!cfg.UseBackground || !cfg.PrefitBackground) return;

    gRejectQminMeV = cfg.BackgroundNonFemtoQminMeV;
    gQAxisScaleToMeV = cfg.qAxisScaleToMeV;
    gSliceCutMeV = cfg.sliceCutAbs*cfg.qAxisScaleToMeV;

    TF3 f("fScanBkg", BkgRejectTF3,
          data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(),
          data->GetYaxis()->GetXmin(), data->GetYaxis()->GetXmax(),
          data->GetZaxis()->GetXmin(), data->GetZaxis()->GetXmax(), 10);
    for(int i=0; i<10; ++i) f.SetParameter(i, bkg[i]);
    f.SetParLimits(0, 0.5, 1.5);
    data->Fit(&f, "RQ0");
    for(int i=0; i<10; ++i) bkg[i] = f.GetParameter(i);
}

double LambdaMain(const TH1* hLambda, const double kstarMeV,
                  const FitCATS3DConfig& cfg)
{
    if(!cfg.UseLambdaFeedDown || !hLambda) return 1.0;
    const double x = kstarMeV/cfg.lambdaAxisScaleToMeV;
    const int bin = std::clamp(hLambda->GetXaxis()->FindBin(x), 1, hLambda->GetNbinsX());
    const double weakScale = (1.0-cfg.weakDecayFraction)*(1.0-cfg.weakDecayFraction);
    double lambda = hLambda->GetBinContent(bin)*weakScale;
    if(!std::isfinite(lambda)) lambda = 1.0;
    return std::clamp(lambda, 0.0, 1.0);
}

void ConfigureKitty(CATS3D& kitty, const TH3D& data, const FitCATS3DConfig& cfg)
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

    CATSparameters src(CATSparameters::tSource3D, 3, true);
    src.SetParameter(0, cfg.RoutInit);
    src.SetParameter(1, cfg.RsideInit);
    src.SetParameter(2, cfg.RlongInit);
    kitty.SetAnaSource(Basics_Source_Gauss3D, src);

    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(cfg.WithCoulomb ? 1 : 0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(cfg.WithCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
}

std::unique_ptr<TH3D> BuildModel(CATS3D& kitty, const TH3D& data,
                                 const TH1* hLambda,
                                 const FitCATS3DConfig& cfg,
                                 const double ro, const double rs, const double rl)
{
    kitty.SetAnaSource(0, ro, true);
    kitty.SetAnaSource(1, rs, true);
    kitty.SetAnaSource(2, rl, true);
    kitty.KillTheCat();

    std::unique_ptr<TH3D> model(dynamic_cast<TH3D*>(data.Clone("hScanModel")));
    model->SetDirectory(nullptr);
    model->Reset("ICES");
    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                const double core = kitty.GetCorrFun(ix-1, iy-1, iz-1);
                const double lambda = LambdaMain(hLambda, KstarMeV(&data, ix, iy, iz, cfg), cfg);
                model->SetBinContent(ix, iy, iz, 1.0 + lambda*(core - 1.0));
            }
        }
    }
    if(cfg.NormalizeModelInNormWindow){
        const double norm = MeanNormWindow(*model, data, cfg);
        if(norm > 0.0 && std::isfinite(norm)) model->Scale(1.0/norm);
    }
    return model;
}

double PmlChi2(const TH3D& modelShape, const TH3D& data,
               const TH3D& same, const TH3D& mixed,
               const double dataNorm, const double* bkg,
               const FitCATS3DConfig& cfg)
{
    double chi2 = 0.0;
    for(int ix=1; ix<=data.GetNbinsX(); ++ix){
        for(int iy=1; iy<=data.GetNbinsY(); ++iy){
            for(int iz=1; iz<=data.GetNbinsZ(); ++iz){
                if(!InsideFitBox(&data, ix, iy, iz, cfg)) continue;
                const double b = cfg.UseBackground ? BkgValue(data, ix, iy, iz, bkg) : bkg[0];
                if(!std::isfinite(b) || b <= 0.0) continue;
                const double model = b*modelShape.GetBinContent(ix, iy, iz);
                if(!std::isfinite(model) || model <= 0.0) continue;

                const double A = same.GetBinContent(ix, iy, iz);
                const double B = mixed.GetBinContent(ix, iy, iz)*dataNorm;
                if(!std::isfinite(A) || !std::isfinite(B) || A < 0.0 || B < 0.0) continue;
                const double total = A+B;
                if(total <= 0.0) continue;
                double logLike = 0.0;
                if(A > 0.0){
                    const double arg = model*total/(A*(model+1.0));
                    if(arg <= 0.0 || !std::isfinite(arg)) continue;
                    logLike += A*std::log(arg);
                }
                if(B > 0.0){
                    const double arg = total/(B*(model+1.0));
                    if(arg <= 0.0 || !std::isfinite(arg)) continue;
                    logLike += B*std::log(arg);
                }
                chi2 += -2.0*logLike;
            }
        }
    }
    return chi2;
}

TGraph* MakeScanGraph(const TString& name,
                      CATS3D& kitty,
                      const TH3D& data,
                      const TH3D& same,
                      const TH3D& mixed,
                      const TH1* hLambda,
                      const double dataNorm,
                      const double* bkg,
                      const FitCATS3DConfig& cfg,
                      const TString& axis,
                      const double ro0,
                      const double rs0,
                      const double rl0,
                      const double rMin,
                      const double rMax,
                      const double rStep)
{
    std::vector<double> xs;
    std::vector<double> ys;
    for(double r=rMin; r<=rMax+1e-9; r+=rStep){
        double ro = ro0;
        double rs = rs0;
        double rl = rl0;
        if(axis=="out") ro = r;
        else if(axis=="side") rs = r;
        else if(axis=="long") rl = r;
        std::unique_ptr<TH3D> model(BuildModel(kitty, data, hLambda, cfg, ro, rs, rl));
        const double chi2 = PmlChi2(*model, data, same, mixed, dataNorm, bkg, cfg);
        xs.push_back(r);
        ys.push_back(chi2);
        std::cout << name << " R=" << r << " chi2=" << chi2 << "\n";
    }
    TGraph* g = new TGraph(int(xs.size()), xs.data(), ys.data());
    g->SetName(name);
    g->SetTitle(name + ";R (fm);#chi^{2}_{PML}");
    return g;
}

TGraph* DeltaGraph(const TGraph* g, const TString& name)
{
    double minY = 1e300;
    for(int i=0; i<g->GetN(); ++i){
        double x, y;
        g->GetPoint(i, x, y);
        if(y < minY) minY = y;
    }
    std::vector<double> xs(g->GetN());
    std::vector<double> ys(g->GetN());
    for(int i=0; i<g->GetN(); ++i){
        double x, y;
        g->GetPoint(i, x, y);
        xs[i] = x;
        ys[i] = y - minY;
    }
    TGraph* out = new TGraph(g->GetN(), xs.data(), ys.data());
    out->SetName(name);
    out->SetTitle(name + ";R (fm);#Delta#chi^{2}_{PML}");
    return out;
}

TH2D* MakeScan2D(const TString& name,
                 CATS3D& kitty,
                 const TH3D& data,
                 const TH3D& same,
                 const TH3D& mixed,
                 const TH1* hLambda,
                 const double dataNorm,
                 const double* bkg,
                 const FitCATS3DConfig& cfg,
                 const TString& axisX,
                 const TString& axisY,
                 const double ro0,
                 const double rs0,
                 const double rl0,
                 const double rMin,
                 const double rMax,
                 const double rStep)
{
    const int nBins = int(std::floor((rMax-rMin)/rStep + 0.5)) + 1;
    const double halfStep = 0.5*rStep;
    TH2D* h = new TH2D(name, name + ";R_{x} (fm);R_{y} (fm);#chi^{2}_{PML}",
                       nBins, rMin-halfStep, rMax+halfStep,
                       nBins, rMin-halfStep, rMax+halfStep);
    h->SetDirectory(nullptr);

    double bestChi2 = 1e300;
    double bestX = 0.0;
    double bestY = 0.0;
    for(int ix=1; ix<=nBins; ++ix){
        const double rx = rMin + (ix-1)*rStep;
        for(int iy=1; iy<=nBins; ++iy){
            const double ry = rMin + (iy-1)*rStep;
            double ro = ro0;
            double rs = rs0;
            double rl = rl0;
            auto setAxis = [&](const TString& axis, const double value) {
                if(axis=="out") ro = value;
                else if(axis=="side") rs = value;
                else if(axis=="long") rl = value;
            };
            setAxis(axisX, rx);
            setAxis(axisY, ry);
            std::unique_ptr<TH3D> model(BuildModel(kitty, data, hLambda, cfg, ro, rs, rl));
            const double chi2 = PmlChi2(*model, data, same, mixed, dataNorm, bkg, cfg);
            h->SetBinContent(ix, iy, chi2);
            if(chi2 < bestChi2){
                bestChi2 = chi2;
                bestX = rx;
                bestY = ry;
            }
        }
        std::cout << name << " row " << ix << "/" << nBins << " done\n";
    }
    std::cout << name << " best: " << axisX << "=" << bestX
              << " " << axisY << "=" << bestY
              << " chi2=" << bestChi2 << "\n";
    return h;
}

TH2D* DeltaHist2D(const TH2D* h, const TString& name)
{
    TH2D* out = dynamic_cast<TH2D*>(h->Clone(name));
    out->SetDirectory(nullptr);
    out->SetTitle(name + ";R_{x} (fm);R_{y} (fm);#Delta#chi^{2}_{PML}");
    double minVal = 1e300;
    for(int ix=1; ix<=out->GetNbinsX(); ++ix){
        for(int iy=1; iy<=out->GetNbinsY(); ++iy){
            const double v = out->GetBinContent(ix, iy);
            if(std::isfinite(v) && v < minVal) minVal = v;
        }
    }
    for(int ix=1; ix<=out->GetNbinsX(); ++ix){
        for(int iy=1; iy<=out->GetNbinsY(); ++iy){
            out->SetBinContent(ix, iy, out->GetBinContent(ix, iy)-minVal);
        }
    }
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.UseModelCache = false;

    double ro0 = 4.0;
    double rs0 = 4.0;
    double rl0 = 4.0;
    double rMin = 2.5;
    double rMax = 6.0;
    double rStep = 0.1;
    if(argc > 1) ro0 = std::atof(argv[1]);
    if(argc > 2) rs0 = std::atof(argv[2]);
    if(argc > 3) rl0 = std::atof(argv[3]);
    if(argc > 4) rMin = std::atof(argv[4]);
    if(argc > 5) rMax = std::atof(argv[5]);
    if(argc > 6) rStep = std::atof(argv[6]);

    TFile fin(cfg.InputFileName, "READ");
    if(fin.IsZombie()) throw std::runtime_error("cannot open input file");
    std::unique_ptr<TH3D> sameRaw(CloneTH3D(fin, cfg.SameHistoName, "hSameRaw"));
    std::unique_ptr<TH3D> mixedRaw(CloneTH3D(fin, cfg.MixedHistoName, "hMixedRaw"));
    std::unique_ptr<TH3D> same(RebinClone(sameRaw.get(), cfg.RebinFactor, "hSameScan"));
    std::unique_ptr<TH3D> mixed(RebinClone(mixedRaw.get(), cfg.RebinFactor, "hMixedScan"));
    std::unique_ptr<TH3D> data(dynamic_cast<TH3D*>(same->Clone("hDataCkScan")));
    data->SetDirectory(nullptr);
    data->Divide(same.get(), mixed.get(), 1.0, 1.0);
    const double dataNorm = NormalizeCk(data.get(), cfg);

    std::unique_ptr<TFile> fLam;
    TH1* hLambda = nullptr;
    if(cfg.UseLambdaFeedDown){
        fLam.reset(TFile::Open(cfg.LambdaFileName, "READ"));
        if(!fLam || fLam->IsZombie()) throw std::runtime_error("cannot open lambda file");
        hLambda = dynamic_cast<TH1*>(fLam->Get(cfg.LambdaHistoName));
        if(!hLambda) throw std::runtime_error("cannot find lambda histogram");
    }

    CATS3D kitty;
    ConfigureKitty(kitty, *data, cfg);
    if(cfg.UseExternalFastKernel){
        const bool loaded = LoadCATS3DFastKernelByQ(
            kitty, cfg.ExternalFastKernelFileName, cfg.ExternalFastKernelTreeName);
        if(!loaded) throw std::runtime_error("failed to load external fast kernel");
    }

    double bkg[10];
    PrefitBkg(data.get(), cfg, bkg);
    std::cout << "Scan fixed point R=(" << ro0 << "," << rs0 << "," << rl0
              << ") range=[" << rMin << "," << rMax << "] step=" << rStep
              << " dataNorm=" << dataNorm
              << " bkgP0=" << bkg[0] << "\n";

    std::unique_ptr<TGraph> gOut(MakeScanGraph("gChi2_Rout", kitty, *data, *same, *mixed,
                                                hLambda, dataNorm, bkg, cfg,
                                                "out", ro0, rs0, rl0, rMin, rMax, rStep));
    std::unique_ptr<TGraph> gSide(MakeScanGraph("gChi2_Rside", kitty, *data, *same, *mixed,
                                                 hLambda, dataNorm, bkg, cfg,
                                                 "side", ro0, rs0, rl0, rMin, rMax, rStep));
    std::unique_ptr<TGraph> gLong(MakeScanGraph("gChi2_Rlong", kitty, *data, *same, *mixed,
                                                 hLambda, dataNorm, bkg, cfg,
                                                 "long", ro0, rs0, rl0, rMin, rMax, rStep));
    std::unique_ptr<TGraph> dOut(DeltaGraph(gOut.get(), "gDeltaChi2_Rout"));
    std::unique_ptr<TGraph> dSide(DeltaGraph(gSide.get(), "gDeltaChi2_Rside"));
    std::unique_ptr<TGraph> dLong(DeltaGraph(gLong.get(), "gDeltaChi2_Rlong"));

    std::unique_ptr<TH2D> hOutSide(MakeScan2D("hChi2_Rout_Rside", kitty, *data, *same, *mixed,
                                               hLambda, dataNorm, bkg, cfg,
                                               "out", "side", ro0, rs0, rl0,
                                               rMin, rMax, rStep));
    std::unique_ptr<TH2D> hOutLong(MakeScan2D("hChi2_Rout_Rlong", kitty, *data, *same, *mixed,
                                               hLambda, dataNorm, bkg, cfg,
                                               "out", "long", ro0, rs0, rl0,
                                               rMin, rMax, rStep));
    std::unique_ptr<TH2D> hSideLong(MakeScan2D("hChi2_Rside_Rlong", kitty, *data, *same, *mixed,
                                                hLambda, dataNorm, bkg, cfg,
                                                "side", "long", ro0, rs0, rl0,
                                                rMin, rMax, rStep));
    std::unique_ptr<TH2D> dOutSide(DeltaHist2D(hOutSide.get(), "hDeltaChi2_Rout_Rside"));
    std::unique_ptr<TH2D> dOutLong(DeltaHist2D(hOutLong.get(), "hDeltaChi2_Rout_Rlong"));
    std::unique_ptr<TH2D> dSideLong(DeltaHist2D(hSideLong.get(), "hDeltaChi2_Rside_Rlong"));

    TFile fout("../OutputFiles/radius_scan_externalKernel.root", "RECREATE");
    gOut->Write();
    gSide->Write();
    gLong->Write();
    dOut->Write();
    dSide->Write();
    dLong->Write();
    hOutSide->Write();
    hOutLong->Write();
    hSideLong->Write();
    dOutSide->Write();
    dOutLong->Write();
    dSideLong->Write();
    data->Write("hDataCkScan");
    fout.Close();

    std::cout << "Wrote ../OutputFiles/radius_scan_externalKernel.root\n";
    return 0;
}
