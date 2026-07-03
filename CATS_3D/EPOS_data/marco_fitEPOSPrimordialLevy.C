// Summary of number of parameters
// normPureGaussian:          R                                 1 param
// normDoubleGaussian:        R_narrow, R_broad, w_broad        3 params
// normSymmetricLevy:         R, alpha                          2 params
// normGaussExpTail:          R_G, L_exp, w_exp                 3 params
// normLogNormalRadiusMixture medianR, sigma_lnR                2 params
// normGammaMixtureQGaussian  R, nu                             2 params
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

void marco_fitEPOSPrimordialLevy(
    TString inputFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/epos_primordial_pion_pair_source_cent_mt_mg.root",
    TString histName = "hRInv_PRF_cent20to30_mT0p40to0p50",
    double fitMin = 0.0,
    double fitMax = 200.0,
    TString outputPrefix = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/epos_primordial_pion_pair_Rside_fit_cent20to30_mT0p40to0p50",
    TString tableFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/source_shape_fit_parameters.tsv")
{
  TFile* inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open " << inputFileName << std::endl;
    return;
  }

  TH1D* hRaw = dynamic_cast<TH1D*>(inputFile->Get(histName));
  if (!hRaw) {
    std::cerr << "ERROR: cannot find histogram " << histName << std::endl;
    inputFile->Close();
    return;
  }

  TH1D* hist = static_cast<TH1D*>(hRaw->Clone("hRInv_PRF_levy_fit_input"));
  hist->SetDirectory(nullptr);
  inputFile->Close();

  const double integral = hist->Integral("width");
  if (integral > 0.) hist->Scale(1. / integral);

  TF1* levy = new TF1("levy_radial_fit",
                      "[0]*x*x*TMath::Exp(-TMath::Power(x/[1],[2]))",
                      fitMin, fitMax);
  levy->SetParNames("norm", "R", "alpha");
  levy->SetParameters(hist->GetMaximum() / 100., 12.0, 1.3);
  levy->SetParLimits(0, 0., 1e9);
  levy->SetParLimits(1, 0.1, 200.);
  levy->SetParLimits(2, 0.1, 2.0);
  levy->SetNpx(1000);

  const bool isRadialDistribution =
      fitMin >= 0.0 ||
      histName.Contains("RInv") ||
      histName.Contains("RAbs") ||
      histName.Contains("FourPiR2S");

  TF1* pureGaussian = new TF1("normPureGaussian",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double radius = p[0];
        if (radius <= 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          double value = TMath::Exp(-0.5 * absX * absX / (radius * radius));
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 600;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(TMath::Abs(x[0])) / norm;
      },
      fitMin, fitMax, 1);
  pureGaussian->SetParNames("R");
  pureGaussian->SetParameter(0, 10.);
  pureGaussian->SetParLimits(0, 0.5, 300.);
  pureGaussian->SetNpx(1000);

  TF1* douGaussian = new TF1("normDoubleGaussian",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double rNarrow = p[0];
        const double rBroad = p[1];
        const double broadWeight = p[2];
        if (rNarrow <= 0. || rBroad <= 0. || broadWeight < 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          double value = TMath::Exp(-0.5 * absX * absX / (rNarrow * rNarrow)) +
                         broadWeight * TMath::Exp(-0.5 * absX * absX / (rBroad * rBroad));
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 600;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(TMath::Abs(x[0])) / norm;
      },
      fitMin, fitMax, 3);
  douGaussian->SetParNames("R_{narrow}", "R_{broad}", "w_{broad}");
  douGaussian->SetParameter(0, 10.);
  douGaussian->SetParameter(1, 45.);
  douGaussian->SetParameter(2, 0.0);
  douGaussian->SetParLimits(0, 0.5, 80.);
  douGaussian->SetParLimits(1, 1.0, 300.);
  douGaussian->SetParLimits(2, 0., 0.);
  douGaussian->SetNpx(1000);

  TF1* symmetricLevy = new TF1("normSymmetricLevy",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double radius = p[0];
        const double alpha = p[1];
        if (radius <= 0. || alpha <= 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          double value = TMath::Exp(-TMath::Power(absX / radius, alpha));
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 600;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(TMath::Abs(x[0])) / norm;
      },
      fitMin, fitMax, 2);
  symmetricLevy->SetParNames("R", "#alpha");
  symmetricLevy->SetParameter(0, 7.);
  symmetricLevy->SetParameter(1, 1.2);
  symmetricLevy->SetParLimits(0, 0.5, 300.);
  symmetricLevy->SetParLimits(1, 1.0, 2.5);
  symmetricLevy->SetNpx(1000);

  TF1* gaussExpTail = new TF1("normGaussExpTail",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double rGauss = p[0];
        const double expLength = p[1];
        const double expWeight = p[2];
        if (rGauss <= 0. || expLength <= 0. || expWeight < 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          double value = TMath::Exp(-0.5 * absX * absX / (rGauss * rGauss)) +
                         expWeight * TMath::Exp(-absX / expLength);
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 600;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(TMath::Abs(x[0])) / norm;
      },
      fitMin, fitMax, 3);
  gaussExpTail->SetParNames("R_{G}", "L_{exp}", "w_{exp}");
  gaussExpTail->SetParameter(0, 7.);
  gaussExpTail->SetParameter(1, 35.);
  gaussExpTail->SetParameter(2, 0.03);
  gaussExpTail->SetParLimits(0, 0.5, 100.);
  gaussExpTail->SetParLimits(1, 1.0, 300.);
  gaussExpTail->SetParLimits(2, 0., 100.);
  gaussExpTail->SetNpx(1000);

  TF1* logNormalRadiusMixture = new TF1("normLogNormalRadiusMixture",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double xx = TMath::Abs(x[0]);
        const double medianR = p[0];
        const double sigmaLogR = p[1];
        if (medianR <= 0. || sigmaLogR <= 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          constexpr int nRadiusSteps = 64;
          const double y0 = TMath::Log(medianR) - 6.0 * sigmaLogR;
          const double y1 = TMath::Log(medianR) + 6.0 * sigmaLogR;
          const double dy = (y1 - y0) / double(nRadiusSteps);
          double weightedGaussianSum = 0.;
          double weightSum = 0.;

          for (int i = 0; i <= nRadiusSteps; ++i) {
            const double y = y0 + dy * double(i);
            const double radius = TMath::Exp(y);
            const double z = (y - TMath::Log(medianR)) / sigmaLogR;
            const double weight = TMath::Exp(-0.5 * z * z) *
                                  ((i == 0 || i == nRadiusSteps) ? 0.5 : 1.0);
            weightedGaussianSum += weight * TMath::Exp(-0.5 * absX * absX / (radius * radius));
            weightSum += weight;
          }

          double value = weightedGaussianSum / weightSum;
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 240;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(xx) / norm;
      },
      fitMin, fitMax, 2);
  logNormalRadiusMixture->SetParNames("median R", "#sigma_{ln R}");
  logNormalRadiusMixture->SetParameter(0, 7.);
  logNormalRadiusMixture->SetParameter(1, 0.6);
  logNormalRadiusMixture->SetParLimits(0, 0.5, 200.);
  logNormalRadiusMixture->SetParLimits(1, 0.02, 10.0);
  logNormalRadiusMixture->SetNpx(1000);

  TF1* gammaMixtureQGaussian = nullptr;
  gammaMixtureQGaussian = new TF1("normGammaMixtureQGaussian",
      [isRadialDistribution, fitMin, fitMax](double* x, double* p) {
        const double radius = p[0];
        const double nu = p[1];
        if (radius <= 0. || nu <= 0.) return 0.0;

        auto rawShape = [&](const double absX) {
          double value = TMath::Power(1.0 + absX * absX / (nu * radius * radius),
                                      -0.5 * (nu + 1.0));
          if (isRadialDistribution) value *= absX * absX;
          return value;
        };

        constexpr int nNormSteps = 600;
        const double dx = (fitMax - fitMin) / double(nNormSteps);
        double norm = 0.0;
        for (int i = 0; i <= nNormSteps; ++i) {
          const double xNorm = fitMin + dx * double(i);
          const double weight = (i == 0 || i == nNormSteps) ? 0.5 : 1.0;
          norm += weight * rawShape(TMath::Abs(xNorm));
        }
        norm *= dx;
        if (norm <= 0.0) return 0.0;
        return rawShape(TMath::Abs(x[0])) / norm;
      },
      fitMin, fitMax, 2);
  gammaMixtureQGaussian->SetParNames("R", "#nu");
  gammaMixtureQGaussian->SetParameter(0, 7.);
  gammaMixtureQGaussian->SetParameter(1, 2.);
  gammaMixtureQGaussian->SetParLimits(0, 0.5, 200.);
  gammaMixtureQGaussian->SetParLimits(1, 0.1, 200.);
  gammaMixtureQGaussian->SetNpx(1000);

  TFitResultPtr fitPureGaussian = hist->Fit(pureGaussian, "RS0");
  TFitResultPtr fitDoubleGaussian = hist->Fit(douGaussian, "RS0+");
  TFitResultPtr fitSymmetricLevy = hist->Fit(symmetricLevy, "RS0+");
  TFitResultPtr fitGaussExpTail = hist->Fit(gaussExpTail, "RS0+");
  TFitResultPtr fitLogNormalMixture = hist->Fit(logNormalRadiusMixture, "RS0+");
  TFitResultPtr fitGammaMixture = hist->Fit(gammaMixtureQGaussian, "RS0+");

  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("cLevyFit", "cLevyFit", 900, 700);
  canvas->SetLogy(true);
  hist->SetTitle("EPOS primordial pion pair source; r^{*}_{inv} (fm); normalized pairs");
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);
  hist->Draw("E");
  pureGaussian->SetLineColor(kCyan + 2);
  pureGaussian->SetLineWidth(2);
  pureGaussian->Draw("same");
  douGaussian->SetLineColor(kRed + 1);
  douGaussian->SetLineWidth(2);
  douGaussian->Draw("same");
  symmetricLevy->SetLineColor(kBlue + 1);
  symmetricLevy->SetLineWidth(2);
  symmetricLevy->Draw("same");
  gaussExpTail->SetLineColor(kGreen + 2);
  gaussExpTail->SetLineWidth(2);
  gaussExpTail->Draw("same");
  logNormalRadiusMixture->SetLineColor(kMagenta + 1);
  logNormalRadiusMixture->SetLineWidth(2);
  logNormalRadiusMixture->Draw("same");
  gammaMixtureQGaussian->SetLineColor(kOrange + 7);
  gammaMixtureQGaussian->SetLineWidth(2);
  gammaMixtureQGaussian->Draw("same");

  TLegend* legend = new TLegend(0.48, 0.68, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist, "EPOS pair source", "lep");
  legend->AddEntry(pureGaussian,
                   isRadialDistribution ?
                   "normalized r^{2}G(R)" :
                   "normalized G(R)",
                   "l");
  legend->AddEntry((TObject*)nullptr, Form("R = %.3f #pm %.3f fm",
                                           pureGaussian->GetParameter(0),
                                           pureGaussian->GetParError(0)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           pureGaussian->GetChisquare(), pureGaussian->GetNDF(),
                                           pureGaussian->GetNDF() > 0 ? pureGaussian->GetChisquare() / pureGaussian->GetNDF() : 0.), "");
  legend->AddEntry(douGaussian,
                   isRadialDistribution ?
                   "normalized r^{2}[G(R_{narrow})+w_{broad}G(R_{broad})]" :
                   "normalized [G(R_{narrow})+w_{broad}G(R_{broad})]",
                   "l");
  legend->AddEntry((TObject*)nullptr, Form("R_{narrow} = %.3f #pm %.3f fm",
                                           douGaussian->GetParameter(0),
                                           douGaussian->GetParError(0)), "");
  legend->AddEntry((TObject*)nullptr, Form("R_{broad} = %.3f #pm %.3f fm",
                                           douGaussian->GetParameter(1),
                                           douGaussian->GetParError(1)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           douGaussian->GetChisquare(), douGaussian->GetNDF(),
                                           douGaussian->GetNDF() > 0 ? douGaussian->GetChisquare() / douGaussian->GetNDF() : 0.), "");
  legend->AddEntry(symmetricLevy, "normalized exp[-(|x|/R)^{#alpha}]", "l");
  legend->AddEntry((TObject*)nullptr, Form("R = %.3f fm, #alpha = %.3f",
                                           symmetricLevy->GetParameter(0),
                                           symmetricLevy->GetParameter(1)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           symmetricLevy->GetChisquare(), symmetricLevy->GetNDF(),
                                           symmetricLevy->GetNDF() > 0 ? symmetricLevy->GetChisquare() / symmetricLevy->GetNDF() : 0.), "");
  legend->AddEntry(gaussExpTail, "normalized [G(R_{G}) + w_{exp}e^{-|x|/L}]", "l");
  legend->AddEntry((TObject*)nullptr, Form("R_{G} = %.3f fm, L = %.3f fm",
                                           gaussExpTail->GetParameter(0),
                                           gaussExpTail->GetParameter(1)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           gaussExpTail->GetChisquare(), gaussExpTail->GetNDF(),
                                           gaussExpTail->GetNDF() > 0 ? gaussExpTail->GetChisquare() / gaussExpTail->GetNDF() : 0.), "");
  legend->AddEntry(logNormalRadiusMixture, "normalized #int dR LogNormal(R) G(x;R)", "l");
  legend->AddEntry((TObject*)nullptr, Form("median R = %.3f fm, #sigma_{lnR} = %.3f",
                                           logNormalRadiusMixture->GetParameter(0),
                                           logNormalRadiusMixture->GetParameter(1)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           logNormalRadiusMixture->GetChisquare(), logNormalRadiusMixture->GetNDF(),
                                           logNormalRadiusMixture->GetNDF() > 0 ? logNormalRadiusMixture->GetChisquare() / logNormalRadiusMixture->GetNDF() : 0.), "");
  legend->AddEntry(gammaMixtureQGaussian, "normalized Gamma-mixture / q-Gaussian", "l");
  legend->AddEntry((TObject*)nullptr, Form("R = %.3f fm, #nu = %.3f",
                                           gammaMixtureQGaussian->GetParameter(0),
                                           gammaMixtureQGaussian->GetParameter(1)), "");
  legend->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                                           gammaMixtureQGaussian->GetChisquare(), gammaMixtureQGaussian->GetNDF(),
                                           gammaMixtureQGaussian->GetNDF() > 0 ? gammaMixtureQGaussian->GetChisquare() / gammaMixtureQGaussian->GetNDF() : 0.), "");
  legend->Draw();

  // canvas->SaveAs(outputPrefix + ".pdf");
  // canvas->SaveAs(outputPrefix + ".png");

  // TFile* outputFile = TFile::Open(outputPrefix + ".root", "RECREATE");
  // hist->Write();
  // levy->Write();
  // pureGaussian->Write();
  // douGaussian->Write();
  // symmetricLevy->Write();
  // gaussExpTail->Write();
  // logNormalRadiusMixture->Write();
  // gammaMixtureQGaussian->Write();

  // outputFile->Close();

  const std::string tablePath = tableFileName.Data();
  const std::string header =
      "centBin\tmTBin\tmodel\t"
      "R_narrow\tR_narrow_err\tR_broad\tR_broad_err\tw_broad\tw_broad_err\t"
      "R\tR_err\talpha\talpha_err\t"
      "R_G\tR_G_err\tL_exp\tL_exp_err\tw_exp\tw_exp_err\t"
      "medianR\tmedianR_err\tsigmaLnR\tsigmaLnR_err\t"
      "nu\tnu_err\t"
      "R_pure\tR_pure_err";

  auto replaceAll = [](std::string value, const std::string& from, const std::string& to) {
    std::size_t pos = 0;
    while ((pos = value.find(from, pos)) != std::string::npos) {
      value.replace(pos, from.length(), to);
      pos += to.length();
    }
    return value;
  };

  auto extractBetween = [](const std::string& text,
                           const std::string& begin,
                           const std::string& end) {
    const std::size_t startPos = text.find(begin);
    if (startPos == std::string::npos) return std::string("unknown");
    const std::size_t valueStart = startPos + begin.size();
    const std::size_t endPos = text.find(end, valueStart);
    if (endPos == std::string::npos) return text.substr(valueStart);
    return text.substr(valueStart, endPos - valueStart);
  };

  auto extractAfter = [](const std::string& text, const std::string& begin) {
    const std::size_t startPos = text.find(begin);
    if (startPos == std::string::npos) return std::string("unknown");
    const std::size_t valueStart = startPos + begin.size();
    const std::size_t endPos = text.find("_", valueStart);
    if (endPos == std::string::npos) return text.substr(valueStart);
    return text.substr(valueStart, endPos - valueStart);
  };

  const std::string histNameStd = histName.Data();
  const std::string centBin = replaceAll(extractBetween(histNameStd, "_cent", "_mT"), "p", ".");
  const std::string mTBin = replaceAll(extractAfter(histNameStd, "_mT"), "p", ".");

  std::vector<std::string> oldLines;
  {
    std::ifstream inputTable(tablePath);
    std::string line;
    while (std::getline(inputTable, line)) {
      if (line.empty()) continue;
      if (line == header) continue;
      if (line.rfind("histName\tfitMin\tfitMax", 0) == 0) continue;
      if (line.rfind("centBin\tmTBin\tmodel", 0) == 0) continue;

      std::stringstream ss(line);
      std::string colCentBin;
      std::string colMtBin;
      std::string colModel;
      std::getline(ss, colCentBin, '\t');
      std::getline(ss, colMtBin, '\t');
      std::getline(ss, colModel, '\t');

      const bool sameCentMtBin =
          colCentBin == centBin &&
          colMtBin == mTBin;
      if (!sameCentMtBin) oldLines.push_back(line);
    }
  }

  std::vector<std::string> newLines;
  const std::string NA = "NA";
  auto fmt = [](const double value) {
    return std::string(Form("%.12g", value));
  };

  auto addWideFitRow = [&](const char* model,
                           const std::vector<std::pair<double,double>>& values) {
    std::vector<std::string> cols(24, NA);
    if (std::string(model) == "normPureGaussian") {
      cols[22] = fmt(values[0].first);
      cols[23] = fmt(values[0].second);
    } else if (std::string(model) == "normDoubleGaussian") {
      cols[0] = fmt(values[0].first);
      cols[1] = fmt(values[0].second);
      cols[2] = fmt(values[1].first);
      cols[3] = fmt(values[1].second);
      cols[4] = fmt(values[2].first);
      cols[5] = fmt(values[2].second);
    } else if (std::string(model) == "normSymmetricLevy") {
      cols[6] = fmt(values[0].first);
      cols[7] = fmt(values[0].second);
      cols[8] = fmt(values[1].first);
      cols[9] = fmt(values[1].second);
    } else if (std::string(model) == "normGaussExpTail") {
      cols[10] = fmt(values[0].first);
      cols[11] = fmt(values[0].second);
      cols[12] = fmt(values[1].first);
      cols[13] = fmt(values[1].second);
      cols[14] = fmt(values[2].first);
      cols[15] = fmt(values[2].second);
    } else if (std::string(model) == "normLogNormalRadiusMixture") {
      cols[16] = fmt(values[0].first);
      cols[17] = fmt(values[0].second);
      cols[18] = fmt(values[1].first);
      cols[19] = fmt(values[1].second);
    } else if (std::string(model) == "normGammaMixtureQGaussian") {
      cols[6] = fmt(values[0].first);
      cols[7] = fmt(values[0].second);
      cols[20] = fmt(values[1].first);
      cols[21] = fmt(values[1].second);
    }

    std::ostringstream row;
    row << centBin << '\t' << mTBin << '\t' << model;
    for (const std::string& col : cols) row << '\t' << col;
    newLines.push_back(row.str());
  };

  // addWideFitRow("normPureGaussian", {
  //     {pureGaussian->GetParameter(0), pureGaussian->GetParError(0)}
  // });
  // addWideFitRow("normDoubleGaussian", {
  //     {douGaussian->GetParameter(0), douGaussian->GetParError(0)},
  //     {douGaussian->GetParameter(1), douGaussian->GetParError(1)},
  //     {douGaussian->GetParameter(2), douGaussian->GetParError(2)}
  // });
  // addWideFitRow("normSymmetricLevy", {
  //     {symmetricLevy->GetParameter(0), symmetricLevy->GetParError(0)},
  //     {symmetricLevy->GetParameter(1), symmetricLevy->GetParError(1)}
  // });
  // addWideFitRow("normGaussExpTail", {
  //     {gaussExpTail->GetParameter(0), gaussExpTail->GetParError(0)},
  //     {gaussExpTail->GetParameter(1), gaussExpTail->GetParError(1)},
  //     {gaussExpTail->GetParameter(2), gaussExpTail->GetParError(2)}
  // });
  // addWideFitRow("normLogNormalRadiusMixture", {
  //     {logNormalRadiusMixture->GetParameter(0), logNormalRadiusMixture->GetParError(0)},
  //     {logNormalRadiusMixture->GetParameter(1), logNormalRadiusMixture->GetParError(1)}
  // });
  // addWideFitRow("normGammaMixtureQGaussian", {
  //     {gammaMixtureQGaussian->GetParameter(0), gammaMixtureQGaussian->GetParError(0)},
  //     {gammaMixtureQGaussian->GetParameter(1), gammaMixtureQGaussian->GetParError(1)}
  // });

  std::ofstream outputTable(tablePath, std::ios::out | std::ios::trunc);
  outputTable << header << '\n';
  for (const std::string& line : oldLines) outputTable << line << '\n';
  for (const std::string& line : newLines) outputTable << line << '\n';
  outputTable.close();
  std::cout << "Updated fit-parameter table: " << tableFileName << std::endl;

  std::cout << "Pure Gaussian fit status = " << int(fitPureGaussian) << std::endl;
  std::cout << "Double Gaussian fit status = " << int(fitDoubleGaussian) << std::endl;
  std::cout << "Symmetric Levy fit status = " << int(fitSymmetricLevy) << std::endl;
  std::cout << "Gaussian + exponential tail fit status = " << int(fitGaussExpTail) << std::endl;
  std::cout << "Log-normal radius mixture fit status = " << int(fitLogNormalMixture) << std::endl;
  std::cout << "Gamma-mixture q-Gaussian fit status = " << int(fitGammaMixture) << std::endl;
  std::cout << "Fit range = [" << fitMin << ", " << fitMax << "] fm" << std::endl;
  std::cout << "Double-Gaussian mode = "
            << (isRadialDistribution ? "radial shell 4pi*r^2*S(r)" : "signed 1D projection")
            << std::endl;
  std::cout << "Normalized pure Gaussian: R = " << pureGaussian->GetParameter(0)
            << " +/- " << pureGaussian->GetParError(0)
            << " fm, chi2/NDF = " << pureGaussian->GetChisquare() << " / " << pureGaussian->GetNDF()
            << " = " << (pureGaussian->GetNDF() > 0 ? pureGaussian->GetChisquare() / pureGaussian->GetNDF() : 0.)
            << std::endl;
  std::cout << "Normalized double Gaussian: R_narrow = " << douGaussian->GetParameter(0)
            << " +/- " << douGaussian->GetParError(0)
            << " fm, R_broad = " << douGaussian->GetParameter(1)
            << " +/- " << douGaussian->GetParError(1)
            << " fm, w_broad = " << douGaussian->GetParameter(2)
            << " +/- " << douGaussian->GetParError(2) << std::endl;
  std::cout << "chi2/NDF = " << douGaussian->GetChisquare() << " / " << douGaussian->GetNDF()
            << " = " << (douGaussian->GetNDF() > 0 ? douGaussian->GetChisquare() / douGaussian->GetNDF() : 0.)
            << std::endl;
  std::cout << "Normalized symmetric Levy: R = " << symmetricLevy->GetParameter(0)
            << " +/- " << symmetricLevy->GetParError(0)
            << " fm, alpha = " << symmetricLevy->GetParameter(1)
            << " +/- " << symmetricLevy->GetParError(1)
            << ", chi2/NDF = " << symmetricLevy->GetChisquare() << " / " << symmetricLevy->GetNDF()
            << " = " << (symmetricLevy->GetNDF() > 0 ? symmetricLevy->GetChisquare() / symmetricLevy->GetNDF() : 0.)
            << std::endl;
  std::cout << "Normalized Gaussian + exp tail: R_G = " << gaussExpTail->GetParameter(0)
            << " +/- " << gaussExpTail->GetParError(0)
            << " fm, L_exp = " << gaussExpTail->GetParameter(1)
            << " +/- " << gaussExpTail->GetParError(1)
            << " fm, w_exp = " << gaussExpTail->GetParameter(2)
            << " +/- " << gaussExpTail->GetParError(2)
            << ", chi2/NDF = " << gaussExpTail->GetChisquare() << " / " << gaussExpTail->GetNDF()
            << " = " << (gaussExpTail->GetNDF() > 0 ? gaussExpTail->GetChisquare() / gaussExpTail->GetNDF() : 0.)
            << std::endl;
  std::cout << "Normalized log-normal radius mixture: median R = " << logNormalRadiusMixture->GetParameter(0)
            << " +/- " << logNormalRadiusMixture->GetParError(0)
            << " fm, sigma_lnR = " << logNormalRadiusMixture->GetParameter(1)
            << " +/- " << logNormalRadiusMixture->GetParError(1)
            << ", chi2/NDF = " << logNormalRadiusMixture->GetChisquare() << " / " << logNormalRadiusMixture->GetNDF()
            << " = " << (logNormalRadiusMixture->GetNDF() > 0 ? logNormalRadiusMixture->GetChisquare() / logNormalRadiusMixture->GetNDF() : 0.)
            << std::endl;
  std::cout << "Normalized gamma-mixture q-Gaussian: R = " << gammaMixtureQGaussian->GetParameter(0)
            << " +/- " << gammaMixtureQGaussian->GetParError(0)
            << " fm, nu = " << gammaMixtureQGaussian->GetParameter(1)
            << " +/- " << gammaMixtureQGaussian->GetParError(1)
            << ", chi2/NDF = " << gammaMixtureQGaussian->GetChisquare() << " / " << gammaMixtureQGaussian->GetNDF()
            << " = " << (gammaMixtureQGaussian->GetNDF() > 0 ? gammaMixtureQGaussian->GetChisquare() / gammaMixtureQGaussian->GetNDF() : 0.)
            << std::endl;
  std::cout << "Saved " << outputPrefix << ".root/.pdf/.png" << std::endl;
}
