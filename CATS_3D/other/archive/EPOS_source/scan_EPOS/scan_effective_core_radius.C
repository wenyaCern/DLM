#include <TFile.h>
#include <TF1.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

TString ExtractBetween(const TString& text, const TString& begin, const TString& end)
{
  const Ssiz_t start = text.Index(begin);
  if (start == kNPOS) return "unknown";
  const Ssiz_t valueStart = start + begin.Length();
  const Ssiz_t stop = text.Index(end, valueStart);
  if (stop == kNPOS) return text(valueStart, text.Length() - valueStart);
  return text(valueStart, stop - valueStart);
}

TString ExtractAfter(const TString& text, const TString& begin)
{
  const Ssiz_t start = text.Index(begin);
  if (start == kNPOS) return "unknown";
  const Ssiz_t valueStart = start + begin.Length();
  const Ssiz_t stop = text.Index(".root", valueStart);
  if (stop == kNPOS) return text(valueStart, text.Length() - valueStart);
  return text(valueStart, stop - valueStart);
}

double MtCenter(const TString& mtBin)
{
  TString clean = mtBin;
  clean.ReplaceAll("p", ".");
  const Ssiz_t pos = clean.Index("to");
  if (pos == kNPOS) return 0.0;
  const double low = TString(clean(0, pos)).Atof();
  const double high = TString(clean(pos + 2, clean.Length() - pos - 2)).Atof();
  return 0.5 * (low + high);
}

struct EffectiveRadiusResult {
  double normFraction = 0.0;
  double meanR2 = 0.0;
  double rRms = 0.0;
  double rCoreEffDiv3 = 0.0;
};

EffectiveRadiusResult ComputeEffectiveRadius(TF1* source,
                                             const double rMin,
                                             const double rMax,
                                             const int nSteps = 6000)
{
  EffectiveRadiusResult result;
  if (!source || rMax <= rMin || nSteps <= 0) return result;

  const double dr = (rMax - rMin) / double(nSteps);
  double norm = 0.0;
  double secondMoment = 0.0;

  for (int i = 0; i <= nSteps; ++i) {
    const double r = rMin + dr * double(i);
    const double weight = (i == 0 || i == nSteps) ? 0.5 : 1.0;
    const double p = std::max(0.0, source->Eval(r));
    norm += weight * p;
    secondMoment += weight * r * r * p;
  }

  norm *= dr;
  secondMoment *= dr;

  result.normFraction = norm;
  if (norm > 0.0) {
    result.meanR2 = secondMoment / norm;
    result.rRms = std::sqrt(result.meanR2);
    result.rCoreEffDiv3 = std::sqrt(result.meanR2 / 3.0);
  }
  return result;
}

void WriteResult(std::ofstream& output,
                 const TString& centBin,
                 const TString& mtBin,
                 const char* model,
                 TF1* source,
                 const double rMin,
                 const double rMax)
{
  if (!source) return;

  const EffectiveRadiusResult result = ComputeEffectiveRadius(source, rMin, rMax);
  output << centBin << '\t'
         << mtBin << '\t'
         << MtCenter(mtBin) << '\t'
         << model << '\t'
         << rMin << '\t'
         << rMax << '\t'
         << Form("%.12g", result.normFraction) << '\t'
         << Form("%.12g", result.meanR2) << '\t'
         << Form("%.12g", result.rRms) << '\t'
         << Form("%.12g", result.rCoreEffDiv3) << '\n';
}

} // namespace

void scan_effective_core_radius(
    TString scanDir = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS",
    TString outputTable = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/effective_core_radius_0to30.tsv",
    double rMin = 0.0,
    double rMax = 30.0)
{
  std::ofstream output(outputTable.Data(), std::ios::out | std::ios::trunc);
  output << "centBin\tmTBin\tmTCenter\tmodel\trMin\trMax\t"
         << "normFractionInRange\tmeanR2InRange\trRMSInRange\tRcoreEffDiv3InRange\n";

  TSystemDirectory dir("scan_EPOS", scanDir);
  TList* files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "ERROR: cannot list directory " << scanDir << std::endl;
    return;
  }

  std::vector<TString> fileNames;
  TIter next(files);
  TSystemFile* file = nullptr;
  while ((file = dynamic_cast<TSystemFile*>(next()))) {
    TString fileName = file->GetName();
    if (file->IsDirectory()) continue;
    if (!fileName.EndsWith(".root")) continue;
    if (!fileName.BeginsWith("epos_primordial_pion_pair_source_fit_cent")) continue;
    fileNames.push_back(fileName);
  }
  std::sort(fileNames.begin(), fileNames.end());

  int nFiles = 0;
  for (const TString& fileName : fileNames) {
    const TString centBin = ExtractBetween(fileName, "_cent", "_mT");
    TString mtBin = ExtractAfter(fileName, "_mT");
    mtBin.ReplaceAll("p", ".");

    TFile inputFile(scanDir + "/" + fileName, "READ");
    if (inputFile.IsZombie()) {
      std::cerr << "WARNING: skip bad ROOT file " << fileName << std::endl;
      continue;
    }

    WriteResult(output,
                centBin,
                mtBin,
                "normLogNormalRadiusMixture",
                dynamic_cast<TF1*>(inputFile.Get("normLogNormalRadiusMixture")),
                rMin,
                rMax);
    WriteResult(output,
                centBin,
                mtBin,
                "normGammaMixtureQGaussian",
                dynamic_cast<TF1*>(inputFile.Get("normGammaMixtureQGaussian")),
                rMin,
                rMax);
    ++nFiles;
  }

  output.close();
  std::cout << "Scanned " << nFiles << " ROOT files.\n";
  std::cout << "Saved " << outputTable << "\n";
  std::cout << "For hRInv radial-shell P(r), RcoreEffDiv3 = sqrt(<r^2>/3).\n";
  std::cout << "rRMS = sqrt(<r^2>) is also stored if you do not want the Gaussian-equivalent /3 convention.\n";
}
