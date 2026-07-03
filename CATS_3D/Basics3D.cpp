#include "Basics3D.h"
#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATStools.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#ifdef FITCATS3D_ENABLE_RESONANCE
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"
#endif

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void Normalize1D(TH1D* h)
{
    if(!h) return;
    const double integral = h->Integral("width");
    if(integral>0. && std::isfinite(integral)){
        h->Scale(1./integral);
    }
}

double Integral3DWithVolume(const TH3D* h)
{
    if(!h) return 0.;
    double sum = 0.;
    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        const double dx = h->GetXaxis()->GetBinWidth(ix);
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            const double dy = h->GetYaxis()->GetBinWidth(iy);
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                const double dz = h->GetZaxis()->GetBinWidth(iz);
                sum += h->GetBinContent(ix,iy,iz)*dx*dy*dz;
            }
        }
    }
    return sum;
}

TH1D* ProjectSource1D(TH3D* h3, const TString& axis, const TString& name)
{
    if(!h3) return nullptr;
    TH1D* h1 = nullptr;
    if(axis=="out"){
        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, 1, h3->GetNbinsY(), 1, h3->GetNbinsZ(), "e"));
        h1->Scale(h3->GetYaxis()->GetBinWidth(1)*h3->GetZaxis()->GetBinWidth(1));
    }
    else if(axis=="side"){
        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, 1, h3->GetNbinsX(), 1, h3->GetNbinsZ(), "e"));
        h1->Scale(h3->GetXaxis()->GetBinWidth(1)*h3->GetZaxis()->GetBinWidth(1));
    }
    else if(axis=="long"){
        h1 = dynamic_cast<TH1D*>(h3->ProjectionZ(name, 1, h3->GetNbinsX(), 1, h3->GetNbinsY(), "e"));
        h1->Scale(h3->GetXaxis()->GetBinWidth(1)*h3->GetYaxis()->GetBinWidth(1));
    }
    if(h1){
        h1->SetDirectory(nullptr);
        Normalize1D(h1);
    }
    return h1;
}

TH1D* ProjectCkFemtoSlice1D(TH3D* h3,
                            const TString& axis,
                            const TString& name,
                            const double cutAbsMeV)
{
    if(!h3) return nullptr;
    TH1D* h1 = nullptr;
    auto centeredBinRange = [](TAxis* axisIn, const double cutAbs, int& firstBin, int& lastBin) {
        firstBin = 1;
        lastBin = 0;
        for(int iBin=1; iBin<=axisIn->GetNbins(); ++iBin){
            if(std::fabs(axisIn->GetBinCenter(iBin)) <= cutAbs + 1e-12){
                if(lastBin == 0) firstBin = iBin;
                lastBin = iBin;
            }
        }
        if(lastBin == 0){
            const int centerBin = axisIn->FindBin(0.0);
            firstBin = centerBin;
            lastBin = centerBin;
        }
    };
    if(axis=="out"){
        int y1, y2, z1, z2;
        centeredBinRange(h3->GetYaxis(), cutAbsMeV, y1, y2);
        centeredBinRange(h3->GetZaxis(), cutAbsMeV, z1, z2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, y1, y2, z1, z2, "e"));
        h1->Scale(1./double((y2-y1+1)*(z2-z1+1)));
    }
    else if(axis=="side"){
        int x1, x2, z1, z2;
        centeredBinRange(h3->GetXaxis(), cutAbsMeV, x1, x2);
        centeredBinRange(h3->GetZaxis(), cutAbsMeV, z1, z2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, x1, x2, z1, z2, "e"));
        h1->Scale(1./double((x2-x1+1)*(z2-z1+1)));
    }
    else if(axis=="long"){
        int x1, x2, y1, y2;
        centeredBinRange(h3->GetXaxis(), cutAbsMeV, x1, x2);
        centeredBinRange(h3->GetYaxis(), cutAbsMeV, y1, y2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionZ(name, x1, x2, y1, y2, "e"));
        h1->Scale(1./double((x2-x1+1)*(y2-y1+1)));
    }
    if(h1) h1->SetDirectory(nullptr);
    return h1;
}

void DrawOverlay(TH1D* gauss,
                 TH1D* reso,
                 const TString& pngName,
                 const bool logY=false)
{
    if(!gauss || !reso) return;
    TCanvas canvas("c_source_compare", "c_source_compare", 900, 650);
    canvas.SetTicks(1,1);
    if(logY) canvas.SetLogy();

    gauss->SetLineColor(kRed+1);
    gauss->SetLineWidth(2);
    reso->SetLineColor(kBlue+1);
    reso->SetLineWidth(2);

    const double maxVal = std::max(gauss->GetMaximum(), reso->GetMaximum());
    gauss->SetMaximum(logY ? maxVal*8. : maxVal*1.25);
    if(logY) gauss->SetMinimum(std::max(1e-8, 0.1*std::min(gauss->GetMinimum(0), reso->GetMinimum(0))));

    gauss->Draw("hist");
    reso->Draw("hist same");

    TLegend leg(0.62,0.76,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(gauss, "Gaussian R=5", "l");
    leg.AddEntry(reso, "DLM_CleverReso3D R=5", "l");
    leg.Draw();

    canvas.SaveAs(pngName);
}

void DrawSingleCk(TH1D* h,
                  const TString& pngName)
{
    if(!h) return;
    TCanvas canvas("c_ck_check", "c_ck_check", 900, 650);
    canvas.SetTicks(1,1);

    h->SetLineColor(kBlue+1);
    h->SetLineWidth(2);
    h->SetMarkerColor(kBlue+1);
    h->SetMarkerStyle(20);

    const double yMin = std::min(0.85, h->GetMinimum()*0.98);
    const double yMax = std::max(1.15, h->GetMaximum()*1.02);
    h->SetMinimum(yMin);
    h->SetMaximum(yMax);
    h->Draw("hist e");

    TLine line(h->GetXaxis()->GetXmin(), 1.0,
               h->GetXaxis()->GetXmax(), 1.0);
    line.SetLineColor(kBlack);
    line.SetLineStyle(2);
    line.Draw("same");

    canvas.SaveAs(pngName);
}

void ConfigureDiagnosticKitty(CATS3D& kitty,
                              const bool useResonance,
#ifdef FITCATS3D_ENABLE_RESONANCE
                              DLM_CleverReso3D* resonanceSource,
#else
                              void* resonanceSource,
#endif
                              const double Rout,
                              const double Rside,
                              const double Rlong,
                              const bool withCoulomb)
{
    (void)resonanceSource;
    const unsigned nQ = 11;
    const double qMin = -200.;
    const double qMax = 200.;
    kitty.SetMomBinsOut(nQ, qMin, qMax);
    kitty.SetMomBinsSide(nQ, qMin, qMax);
    kitty.SetMomBinsLong(nQ, qMin, qMax);

    const unsigned nR = 15;
    const double rMin = -80.;
    const double rMax = 80.;
    kitty.SetGridX(nR, rMin, rMax);
    kitty.SetGridY(nR, rMin, rMax);
    kitty.SetGridZ(nR, rMin, rMax);
    kitty.SetFastKernelMode(true);

    if(useResonance){
#ifdef FITCATS3D_ENABLE_RESONANCE
        kitty.SetAnaSource(CatsSourceForwarder, resonanceSource, 3);
        kitty.SetAnaSource(0, Rout);
        kitty.SetAnaSource(1, Rside);
        kitty.SetAnaSource(2, Rlong);
#else
        throw std::runtime_error("ConfigureDiagnosticKitty: rebuild with FITCATS3D_ENABLE_RESONANCE for resonance source.");
#endif
    }
    else{
        CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
        sourcePars.SetParameter(0, Rout);
        sourcePars.SetParameter(1, Rside);
        sourcePars.SetParameter(2, Rlong);
        kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);
    }

    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(withCoulomb ? 1 : 0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(withCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
}

TH3D* LoadDataBinningTemplate3D(const FitCATS3DConfig& cfg)
{
    TFile fin(cfg.InputFileName, "READ");
    if(fin.IsZombie()){
        throw std::runtime_error(Form("LoadDataBinningTemplate3D: cannot open '%s'",
                                      cfg.InputFileName.Data()));
    }

    TH3D* hIn = dynamic_cast<TH3D*>(fin.Get(cfg.SameHistoName));
    if(!hIn){
        throw std::runtime_error(Form("LoadDataBinningTemplate3D: missing '%s'",
                                      cfg.SameHistoName.Data()));
    }

    TH3D* h = dynamic_cast<TH3D*>(hIn->Clone("hDataBinningTemplate3D"));
    h->SetDirectory(nullptr);
    if(cfg.RebinFactor > 1){
        h->Rebin3D(cfg.RebinFactor, cfg.RebinFactor, cfg.RebinFactor);
    }
    h->Reset("ICES");
    return h;
}

} // namespace

// --------------------------------------------------
// 3D Gaussian OSL source
// Pars[0] = q_out
// Pars[1] = q_side
// Pars[2] = q_long
// Pars[3] = r_out
// Pars[4] = r_side
// Pars[5] = r_long
// Pars[6] = Rout
// Pars[7] = Rside
// Pars[8] = Rlong
// --------------------------------------------------
double Basics_Source_Gauss3D(double* Pars){
    const double& rOut  = Pars[3];
    const double& rSide = Pars[4];
    const double& rLong = Pars[5];

    const double& Rout  = Pars[6];
    const double& Rside = Pars[7];
    const double& Rlong = Pars[8];

    return 1.0 / (std::pow(4.0*Pi, 1.5) * Rout * Rside * Rlong) *
           std::exp( - rOut*rOut   / (4.0*Rout*Rout)
                    -rSide*rSide   / (4.0*Rside*Rside)
                    -rLong*rLong   / (4.0*Rlong*Rlong) );
}

// --------------------------------------------------
// Setup only
// --------------------------------------------------
void Basics_PiPiCATS3D_Setup(CATS3D& Kitty,
                             const bool& Identical,
                             const bool& WithCoulomb){

    // -----------------------------
    // binning
    // -----------------------------
    const unsigned NumMomBinsOut  = 21;
    const unsigned NumMomBinsSide = 21;
    const unsigned NumMomBinsLong = 21;

    const double qOutMin  = -500.0;
    const double qOutMax  =  500.0;
    const double qSideMin = -500.0;
    const double qSideMax =  500.0;
    const double qLongMin = -500.0;
    const double qLongMax =  500.0;

    // -----------------------------
    // coordinate-space source grid
    // -----------------------------
    const unsigned NumGridBinsX = 21;
    const unsigned NumGridBinsY = 21;
    const unsigned NumGridBinsZ = 21;

    const double xMin = -20.0;
    const double xMax =  20.0;
    const double yMin = -20.0;
    const double yMax =  20.0;
    const double zMin = -20.0;
    const double zMax =  20.0;

    // -----------------------------
    // Gaussian OSL radii
    // -----------------------------
    const double Rout  = 5.0;
    const double Rside = 5.0;
    const double Rlong = 5.0;

    // momentum bins
    Kitty.SetMomBinsOut (NumMomBinsOut , qOutMin , qOutMax );
    Kitty.SetMomBinsSide(NumMomBinsSide, qSideMin, qSideMax);
    Kitty.SetMomBinsLong(NumMomBinsLong, qLongMin, qLongMax);

    // source grid
    Kitty.SetGridX(NumGridBinsX, xMin, xMax);
    Kitty.SetGridY(NumGridBinsY, yMin, yMax);
    Kitty.SetGridZ(NumGridBinsZ, zMin, zMax);

    // analytic source
    CATSparameters SOURCE_PARS(CATSparameters::tSource3D, 3, true);
    SOURCE_PARS.SetParameter(0, Rout);
    SOURCE_PARS.SetParameter(1, Rside);
    SOURCE_PARS.SetParameter(2, Rlong);
    Kitty.SetAnaSource(Basics_Source_Gauss3D, SOURCE_PARS);

    // physics setup
    if(Identical) Kitty.SetPdgId(211, 211);
    else          Kitty.SetPdgId(-211, 211);

    if(Identical) Kitty.SetQuantumStatistics(1);
    else          Kitty.SetQuantumStatistics(0);

    if(WithCoulomb) Kitty.SetQ1Q2(1);
    else            Kitty.SetQ1Q2(0);

    // for now both branches use charged pion mass as the simple test setup
    Kitty.SetRedMass(0.5*Mass_Pic);

    if(WithCoulomb) Kitty.SetWfMode(CATS3D::wfCoulomb);
    else            Kitty.SetWfMode(CATS3D::wfPlaneWave);
}

// --------------------------------------------------
// Run CATS3D and convert result to TH3D
// --------------------------------------------------
TH3D* Basics_CATS3D_GetCkHisto(CATS3D& Kitty,
                               const TString& HistoName){

    Kitty.KillTheCat();

    const unsigned NumMomBinsOut  = Kitty.GetNumMomBinsOut();
    const unsigned NumMomBinsSide = Kitty.GetNumMomBinsSide();
    const unsigned NumMomBinsLong = Kitty.GetNumMomBinsLong();

    if(!NumMomBinsOut || !NumMomBinsSide || !NumMomBinsLong){
        throw std::runtime_error("Basics_CATS3D_GetCkHisto: momentum bins are not configured.");
    }

    double* qOutBins  = new double[NumMomBinsOut+1];
    double* qSideBins = new double[NumMomBinsSide+1];
    double* qLongBins = new double[NumMomBinsLong+1];

    for(unsigned uBin=0; uBin<NumMomBinsOut; uBin++){
        qOutBins[uBin] = Kitty.GetMomBinLowEdgeOut(uBin);
    }
    qOutBins[NumMomBinsOut] = Kitty.GetMomBinUpEdgeOut(NumMomBinsOut-1);

    for(unsigned uBin=0; uBin<NumMomBinsSide; uBin++){
        qSideBins[uBin] = Kitty.GetMomBinLowEdgeSide(uBin);
    }
    qSideBins[NumMomBinsSide] = Kitty.GetMomBinUpEdgeSide(NumMomBinsSide-1);

    for(unsigned uBin=0; uBin<NumMomBinsLong; uBin++){
        qLongBins[uBin] = Kitty.GetMomBinLowEdgeLong(uBin);
    }
    qLongBins[NumMomBinsLong] = Kitty.GetMomBinUpEdgeLong(NumMomBinsLong-1);

    TH3D* hCk3D = new TH3D(HistoName,
                           "C(q_{out},q_{side},q_{long});q_{out} (MeV/c);q_{side} (MeV/c);q_{long} (MeV/c)",
                           NumMomBinsOut,  qOutBins,
                           NumMomBinsSide, qSideBins,
                           NumMomBinsLong, qLongBins);

    for(unsigned uOut=0; uOut<NumMomBinsOut; uOut++){
        for(unsigned uSide=0; uSide<NumMomBinsSide; uSide++){
            for(unsigned uLong=0; uLong<NumMomBinsLong; uLong++){
                hCk3D->SetBinContent(uOut+1, uSide+1, uLong+1,
                                     Kitty.GetCorrFun(uOut,uSide,uLong));
                hCk3D->SetBinError(uOut+1, uSide+1, uLong+1,
                                   Kitty.GetCorrFunErr(uOut,uSide,uLong));
            }
        }
    }

    delete [] qOutBins;
    delete [] qSideBins;
    delete [] qLongBins;

    return hCk3D;
}

// --------------------------------------------------
// Save TH3D to ROOT
// --------------------------------------------------
void Save_TH3D_to_ROOT(const TString& FileName,
                       TH3D* histo,
                       const TString& HistoName){

    if(!histo){
        throw std::runtime_error("Save_TH3D_to_ROOT: input histogram is null.");
    }

    TFile* fout = TFile::Open(FileName, "UPDATE");
    if(!fout || fout->IsZombie()){
        if(fout){
            fout->Close();
            delete fout;
        }
        fout = TFile::Open(FileName, "RECREATE");
    }

    if(!fout || fout->IsZombie()){
        throw std::runtime_error("Save_TH3D_to_ROOT: cannot open output file.");
    }

    fout->cd();

    TH3D* hToWrite = (TH3D*)histo->Clone();
    if(HistoName.Length()){
        hToWrite->SetName(HistoName);
        hToWrite->SetTitle(HistoName);
    }

    hToWrite->Write("", TObject::kOverwrite);

    fout->Close();
    delete fout;
}

void Basics_RunGaussianR20DataRangeCheck3D(const TString& outFileName,
                                           const double qAbsMaxMeV,
                                           const unsigned nQBins,
                                           const bool runCoulomb,
                                           const bool runNoCoulomb)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    const double Rout = 20.0;
    const double Rside = 20.0;
    const double Rlong = 20.0;

    std::unique_ptr<TH3D> hTemplate(LoadDataBinningTemplate3D(cfg));

    auto configureKitty = [&](CATS3D& kitty, const bool withCoulomb) {
        kitty.SetMomBinsOut(nQBins, -qAbsMaxMeV, qAbsMaxMeV);
        kitty.SetMomBinsSide(nQBins, -qAbsMaxMeV, qAbsMaxMeV);
        kitty.SetMomBinsLong(nQBins, -qAbsMaxMeV, qAbsMaxMeV);

        kitty.SetGridX(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
        kitty.SetGridY(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
        kitty.SetGridZ(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
        kitty.SetFastKernelMode(cfg.UseFastKernel);

        CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
        sourcePars.SetParameter(0, Rout);
        sourcePars.SetParameter(1, Rside);
        sourcePars.SetParameter(2, Rlong);
        kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);

        kitty.SetPdgId(211, 211);
        kitty.SetQuantumStatistics(1);
        kitty.SetQ1Q2(withCoulomb ? 1 : 0);
        kitty.SetRedMass(0.5*Mass_Pic);
        kitty.SetWfMode(withCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
    };

    std::unique_ptr<TH3D> hCoulomb;
    if(runCoulomb){
        CATS3D kittyCoulomb;
        configureKitty(kittyCoulomb, true);
        hCoulomb.reset(Basics_CATS3D_GetCkHisto(
            kittyCoulomb, "hCkGauss_R20_dataRange_coulomb"));
        hCoulomb->SetDirectory(nullptr);
    }

    std::unique_ptr<TH3D> hNoCoulomb;
    if(runNoCoulomb){
        CATS3D kittyNoCoulomb;
        configureKitty(kittyNoCoulomb, false);
        hNoCoulomb.reset(Basics_CATS3D_GetCkHisto(
            kittyNoCoulomb, "hCkGauss_R20_dataRange_noCoulomb"));
        hNoCoulomb->SetDirectory(nullptr);
    }

    const double qCutMeV = cfg.sliceCutAbs * cfg.qAxisScaleToMeV;
    std::unique_ptr<TH1D> hCoulombOut;
    std::unique_ptr<TH1D> hCoulombSide;
    std::unique_ptr<TH1D> hCoulombLong;
    if(hCoulomb){
        hCoulombOut.reset(ProjectCkFemtoSlice1D(
            hCoulomb.get(), "out", "hCkGauss_R20_coulomb_proj_out", qCutMeV));
        hCoulombSide.reset(ProjectCkFemtoSlice1D(
            hCoulomb.get(), "side", "hCkGauss_R20_coulomb_proj_side", qCutMeV));
        hCoulombLong.reset(ProjectCkFemtoSlice1D(
            hCoulomb.get(), "long", "hCkGauss_R20_coulomb_proj_long", qCutMeV));
    }

    std::unique_ptr<TH1D> hNoCoulombOut;
    std::unique_ptr<TH1D> hNoCoulombSide;
    std::unique_ptr<TH1D> hNoCoulombLong;
    if(hNoCoulomb){
        hNoCoulombOut.reset(ProjectCkFemtoSlice1D(
            hNoCoulomb.get(), "out", "hCkGauss_R20_noCoulomb_proj_out", qCutMeV));
        hNoCoulombSide.reset(ProjectCkFemtoSlice1D(
            hNoCoulomb.get(), "side", "hCkGauss_R20_noCoulomb_proj_side", qCutMeV));
        hNoCoulombLong.reset(ProjectCkFemtoSlice1D(
            hNoCoulomb.get(), "long", "hCkGauss_R20_noCoulomb_proj_long", qCutMeV));
    }

    TFile fout(outFileName, "RECREATE");
    hTemplate->Write("hDataBinningTemplate3D");
    if(hCoulomb) hCoulomb->Write();
    if(hNoCoulomb) hNoCoulomb->Write();
    if(hCoulombOut) hCoulombOut->Write();
    if(hCoulombSide) hCoulombSide->Write();
    if(hCoulombLong) hCoulombLong->Write();
    if(hNoCoulombOut) hNoCoulombOut->Write();
    if(hNoCoulombSide) hNoCoulombSide->Write();
    if(hNoCoulombLong) hNoCoulombLong->Write();
    fout.Close();

    const TString outDir = outFileName.Contains("/") ?
        outFileName(0, outFileName.Last('/')) : TString(".");
    DrawSingleCk(hCoulombOut.get(), outDir + "/gauss_R20_dataRange_coulomb_proj_out.png");
    DrawSingleCk(hCoulombSide.get(), outDir + "/gauss_R20_dataRange_coulomb_proj_side.png");
    DrawSingleCk(hCoulombLong.get(), outDir + "/gauss_R20_dataRange_coulomb_proj_long.png");
    DrawSingleCk(hNoCoulombOut.get(), outDir + "/gauss_R20_dataRange_noCoulomb_proj_out.png");
    DrawSingleCk(hNoCoulombSide.get(), outDir + "/gauss_R20_dataRange_noCoulomb_proj_side.png");
    DrawSingleCk(hNoCoulombLong.get(), outDir + "/gauss_R20_dataRange_noCoulomb_proj_long.png");

    std::cout << "R20 data-range check wrote " << outFileName << "\n"
              << "q range = [-" << qAbsMaxMeV << ", " << qAbsMaxMeV
              << "] MeV/c with " << nQBins << " bins per axis\n"
              << "runCoulomb=" << runCoulomb
              << " runNoCoulomb=" << runNoCoulomb << "\n";
    if(hCoulombOut){
        std::cout << "Coulomb out min/max = " << hCoulombOut->GetMinimum()
                  << "/" << hCoulombOut->GetMaximum()
                  << ", edge bins = " << hCoulombOut->GetBinContent(1)
                  << ", " << hCoulombOut->GetBinContent(hCoulombOut->GetNbinsX())
                  << "\n";
    }
    if(hNoCoulombOut){
        std::cout << "NoCoulomb out min/max = " << hNoCoulombOut->GetMinimum()
                  << "/" << hNoCoulombOut->GetMaximum()
                  << ", edge bins = " << hNoCoulombOut->GetBinContent(1)
                  << ", " << hNoCoulombOut->GetBinContent(hNoCoulombOut->GetNbinsX())
                  << "\n";
    }
}

void Basics_RunSourceComparison3D(const TString& outFileName)
{
    const double Rout = 5.;
    const double Rside = 5.;
    const double Rlong = 5.;

#ifndef FITCATS3D_ENABLE_RESONANCE
    throw std::runtime_error("Basics_RunSourceComparison3D: use the resonance build target so DLM_CleverReso3D and EPOS input are available.");
#else
    FitCATS3DConfig cfg = MakeFitCATS3DResonanceConfig();
    cfg.UseLambdaFeedDown = false;
    cfg.WithCoulomb = false;
    cfg.RoutInit = Rout;
    cfg.RsideInit = Rside;
    cfg.RlongInit = Rlong;
    cfg.ResonanceRadBins = 101;
    cfg.ResonanceRadMinFm = -80.;
    cfg.ResonanceRadMaxFm = 80.;
    cfg.ResonanceMcIter = 400000;

    std::unique_ptr<DLM_CleverReso3D> resoSource(BuildCATS3DResonanceSource(cfg));
    if(!resoSource){
        throw std::runtime_error("Basics_RunSourceComparison3D: failed to build DLM_CleverReso3D.");
    }

    const unsigned nR = 81;
    const double rMin = -80.;
    const double rMax = 80.;

    std::unique_ptr<TH3D> hSourceGauss(new TH3D(
        "hSourceGauss3D_R5",
        "Gaussian source R=5;r_{out} (fm);r_{side} (fm);r_{long} (fm)",
        nR, rMin, rMax, nR, rMin, rMax, nR, rMin, rMax));
    std::unique_ptr<TH3D> hSourceReso(new TH3D(
        "hSourceReso3D_R5",
        "DLM_CleverReso3D source R=5;r_{out} (fm);r_{side} (fm);r_{long} (fm)",
        nR, rMin, rMax, nR, rMin, rMax, nR, rMin, rMax));
    hSourceGauss->SetDirectory(nullptr);
    hSourceReso->SetDirectory(nullptr);

    double pars[9] = {};
    pars[6] = Rout;
    pars[7] = Rside;
    pars[8] = Rlong;

    for(unsigned ix=1; ix<=nR; ++ix){
        pars[3] = hSourceGauss->GetXaxis()->GetBinCenter(ix);
        for(unsigned iy=1; iy<=nR; ++iy){
            pars[4] = hSourceGauss->GetYaxis()->GetBinCenter(iy);
            for(unsigned iz=1; iz<=nR; ++iz){
                pars[5] = hSourceGauss->GetZaxis()->GetBinCenter(iz);
                hSourceGauss->SetBinContent(ix, iy, iz,
                                            Basics_Source_Gauss3D(pars));
                hSourceReso->SetBinContent(ix, iy, iz,
                                           resoSource->Eval(pars));
            }
        }
    }

    std::cout << "Source integral check, Gaussian="
              << Integral3DWithVolume(hSourceGauss.get())
              << " resonance="
              << Integral3DWithVolume(hSourceReso.get())
              << "\n";

    CATS3D kittyGauss;
    ConfigureDiagnosticKitty(kittyGauss, false, nullptr, Rout, Rside, Rlong, false);
    std::unique_ptr<TH3D> hCkGauss(Basics_CATS3D_GetCkHisto(kittyGauss, "hCkGauss3D_R5"));
    hCkGauss->SetDirectory(nullptr);

    CATS3D kittyReso;
    ConfigureDiagnosticKitty(kittyReso, true, resoSource.get(), Rout, Rside, Rlong, false);
    std::unique_ptr<TH3D> hCkReso(Basics_CATS3D_GetCkHisto(kittyReso, "hCkReso3D_R5"));
    hCkReso->SetDirectory(nullptr);

    CATS3D kittyGaussCoulomb;
    ConfigureDiagnosticKitty(kittyGaussCoulomb, false, nullptr, Rout, Rside, Rlong, true);
    std::unique_ptr<TH3D> hCkGaussCoulomb(
        Basics_CATS3D_GetCkHisto(kittyGaussCoulomb, "hCkGauss3D_R5_coulomb"));
    hCkGaussCoulomb->SetDirectory(nullptr);

    CATS3D kittyResoCoulomb;
    ConfigureDiagnosticKitty(kittyResoCoulomb, true, resoSource.get(), Rout, Rside, Rlong, true);
    std::unique_ptr<TH3D> hCkResoCoulomb(
        Basics_CATS3D_GetCkHisto(kittyResoCoulomb, "hCkReso3D_R5_coulomb"));
    hCkResoCoulomb->SetDirectory(nullptr);

    std::unique_ptr<TH1D> hSrcGaussOut(ProjectSource1D(hSourceGauss.get(), "out", "hSourceGauss_proj_out"));
    std::unique_ptr<TH1D> hSrcGaussSide(ProjectSource1D(hSourceGauss.get(), "side", "hSourceGauss_proj_side"));
    std::unique_ptr<TH1D> hSrcGaussLong(ProjectSource1D(hSourceGauss.get(), "long", "hSourceGauss_proj_long"));
    std::unique_ptr<TH1D> hSrcResoOut(ProjectSource1D(hSourceReso.get(), "out", "hSourceReso_proj_out"));
    std::unique_ptr<TH1D> hSrcResoSide(ProjectSource1D(hSourceReso.get(), "side", "hSourceReso_proj_side"));
    std::unique_ptr<TH1D> hSrcResoLong(ProjectSource1D(hSourceReso.get(), "long", "hSourceReso_proj_long"));

    const double qCut = 60.;
    std::unique_ptr<TH1D> hCkGaussOut(ProjectCkFemtoSlice1D(hCkGauss.get(), "out", "hCkGauss_proj_out", qCut));
    std::unique_ptr<TH1D> hCkGaussSide(ProjectCkFemtoSlice1D(hCkGauss.get(), "side", "hCkGauss_proj_side", qCut));
    std::unique_ptr<TH1D> hCkGaussLong(ProjectCkFemtoSlice1D(hCkGauss.get(), "long", "hCkGauss_proj_long", qCut));
    std::unique_ptr<TH1D> hCkResoOut(ProjectCkFemtoSlice1D(hCkReso.get(), "out", "hCkReso_proj_out", qCut));
    std::unique_ptr<TH1D> hCkResoSide(ProjectCkFemtoSlice1D(hCkReso.get(), "side", "hCkReso_proj_side", qCut));
    std::unique_ptr<TH1D> hCkResoLong(ProjectCkFemtoSlice1D(hCkReso.get(), "long", "hCkReso_proj_long", qCut));
    std::unique_ptr<TH1D> hCkGaussCoulombOut(ProjectCkFemtoSlice1D(hCkGaussCoulomb.get(), "out", "hCkGauss_coulomb_proj_out", qCut));
    std::unique_ptr<TH1D> hCkGaussCoulombSide(ProjectCkFemtoSlice1D(hCkGaussCoulomb.get(), "side", "hCkGauss_coulomb_proj_side", qCut));
    std::unique_ptr<TH1D> hCkGaussCoulombLong(ProjectCkFemtoSlice1D(hCkGaussCoulomb.get(), "long", "hCkGauss_coulomb_proj_long", qCut));
    std::unique_ptr<TH1D> hCkResoCoulombOut(ProjectCkFemtoSlice1D(hCkResoCoulomb.get(), "out", "hCkReso_coulomb_proj_out", qCut));
    std::unique_ptr<TH1D> hCkResoCoulombSide(ProjectCkFemtoSlice1D(hCkResoCoulomb.get(), "side", "hCkReso_coulomb_proj_side", qCut));
    std::unique_ptr<TH1D> hCkResoCoulombLong(ProjectCkFemtoSlice1D(hCkResoCoulomb.get(), "long", "hCkReso_coulomb_proj_long", qCut));

    TFile fout(outFileName, "RECREATE");
    hSourceGauss->Write();
    hSourceReso->Write();
    hCkGauss->Write();
    hCkReso->Write();
    hCkGaussCoulomb->Write();
    hCkResoCoulomb->Write();

    hSrcGaussOut->Write();
    hSrcGaussSide->Write();
    hSrcGaussLong->Write();
    hSrcResoOut->Write();
    hSrcResoSide->Write();
    hSrcResoLong->Write();

    hCkGaussOut->Write();
    hCkGaussSide->Write();
    hCkGaussLong->Write();
    hCkResoOut->Write();
    hCkResoSide->Write();
    hCkResoLong->Write();
    hCkGaussCoulombOut->Write();
    hCkGaussCoulombSide->Write();
    hCkGaussCoulombLong->Write();
    hCkResoCoulombOut->Write();
    hCkResoCoulombSide->Write();
    hCkResoCoulombLong->Write();
    fout.Close();

    const TString outDir = outFileName(0, outFileName.Last('/'));
    DrawOverlay(hSrcGaussOut.get(), hSrcResoOut.get(), outDir + "/source_proj_out_R5.png", true);
    DrawOverlay(hSrcGaussSide.get(), hSrcResoSide.get(), outDir + "/source_proj_side_R5.png", true);
    DrawOverlay(hSrcGaussLong.get(), hSrcResoLong.get(), outDir + "/source_proj_long_R5.png", true);
    DrawOverlay(hCkGaussOut.get(), hCkResoOut.get(), outDir + "/ck_proj_out_R5.png", false);
    DrawOverlay(hCkGaussSide.get(), hCkResoSide.get(), outDir + "/ck_proj_side_R5.png", false);
    DrawOverlay(hCkGaussLong.get(), hCkResoLong.get(), outDir + "/ck_proj_long_R5.png", false);
    DrawOverlay(hCkGaussCoulombOut.get(), hCkResoCoulombOut.get(), outDir + "/ck_coulomb_proj_out_R5.png", false);
    DrawOverlay(hCkGaussCoulombSide.get(), hCkResoCoulombSide.get(), outDir + "/ck_coulomb_proj_side_R5.png", false);
    DrawOverlay(hCkGaussCoulombLong.get(), hCkResoCoulombLong.get(), outDir + "/ck_coulomb_proj_long_R5.png", false);

    std::cout << "Wrote source/CF comparison to " << outFileName << "\n";
#endif
}

void Basics_RunFineGaussianNoCoulomb3D(const TString& outFileName,
                                       const double radiusFm,
                                       const double qAbsMaxMeV,
                                       const unsigned nQBins,
                                       const double sliceCutMeV,
                                       const bool useFastKernel,
                                       const unsigned fastSubdivisions)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();

    CATS3D kitty;
    kitty.SetMomBinsOut(nQBins, -qAbsMaxMeV, qAbsMaxMeV);
    kitty.SetMomBinsSide(nQBins, -qAbsMaxMeV, qAbsMaxMeV);
    kitty.SetMomBinsLong(nQBins, -qAbsMaxMeV, qAbsMaxMeV);

    kitty.SetGridX(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridY(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridZ(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetFastKernelMode(useFastKernel);
    kitty.SetFastKernelSubdivisions(fastSubdivisions);

    CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
    sourcePars.SetParameter(0, radiusFm);
    sourcePars.SetParameter(1, radiusFm);
    sourcePars.SetParameter(2, radiusFm);
    kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);

    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(CATS3D::wfPlaneWave);

    std::unique_ptr<TH3D> hCk(Basics_CATS3D_GetCkHisto(
        kitty, Form("hCkFineGauss_noCoulomb_R%.2f", radiusFm)));
    hCk->SetDirectory(nullptr);

    std::unique_ptr<TH1D> hOut(ProjectCkFemtoSlice1D(
        hCk.get(), "out", "hCkFineGauss_noCoulomb_proj_out", sliceCutMeV));
    std::unique_ptr<TH1D> hSide(ProjectCkFemtoSlice1D(
        hCk.get(), "side", "hCkFineGauss_noCoulomb_proj_side", sliceCutMeV));
    std::unique_ptr<TH1D> hLong(ProjectCkFemtoSlice1D(
        hCk.get(), "long", "hCkFineGauss_noCoulomb_proj_long", sliceCutMeV));

    TFile fout(outFileName, "RECREATE");
    hCk->Write();
    hOut->Write();
    hSide->Write();
    hLong->Write();
    fout.Close();

    const TString outDir = outFileName.Contains("/") ?
        outFileName(0, outFileName.Last('/')) : TString(".");
    DrawSingleCk(hOut.get(), outDir + "/fine_gaussian_noCoulomb_proj_out.png");
    DrawSingleCk(hSide.get(), outDir + "/fine_gaussian_noCoulomb_proj_side.png");
    DrawSingleCk(hLong.get(), outDir + "/fine_gaussian_noCoulomb_proj_long.png");

    std::cout << "Fine pure-Gaussian no-Coulomb CATS3D test wrote "
              << outFileName << "\n"
              << "R = " << radiusFm << " fm, q range = [-"
              << qAbsMaxMeV << ", " << qAbsMaxMeV << "] MeV/c, nQBins = "
              << nQBins << ", bin width = "
              << (2.0*qAbsMaxMeV/double(nQBins)) << " MeV/c, slice |q| <= "
              << sliceCutMeV << " MeV/c, fastKernel = "
              << (useFastKernel ? "true" : "false")
              << ", fastSubdivisions = " << fastSubdivisions << "\n"
              << "out projection min/max = " << hOut->GetMinimum()
              << "/" << hOut->GetMaximum() << "\n";
}
