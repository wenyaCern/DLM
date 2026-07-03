#include <TFile.h>
#include <TF1.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>

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

std::string Fmt(const double value)
{
  return std::string(Form("%.12g", value));
}

void WriteRow(std::ofstream& out,
              const TString& centBin,
              const TString& mtBin,
              const char* model,
              TF1* function)
{
  if (!function) return;

  const std::string NA = "NA";
  std::vector<std::string> cols(22, NA);

  const std::string modelName(model);
  if (modelName == "normDoubleGaussian") {
    cols[0] = Fmt(function->GetParameter(0));
    cols[1] = Fmt(function->GetParError(0));
    cols[2] = Fmt(function->GetParameter(1));
    cols[3] = Fmt(function->GetParError(1));
    cols[4] = Fmt(function->GetParameter(2));
    cols[5] = Fmt(function->GetParError(2));
  } else if (modelName == "normSymmetricLevy") {
    cols[6] = Fmt(function->GetParameter(0));
    cols[7] = Fmt(function->GetParError(0));
    cols[8] = Fmt(function->GetParameter(1));
    cols[9] = Fmt(function->GetParError(1));
  } else if (modelName == "normGaussExpTail") {
    cols[10] = Fmt(function->GetParameter(0));
    cols[11] = Fmt(function->GetParError(0));
    cols[12] = Fmt(function->GetParameter(1));
    cols[13] = Fmt(function->GetParError(1));
    cols[14] = Fmt(function->GetParameter(2));
    cols[15] = Fmt(function->GetParError(2));
  } else if (modelName == "normLogNormalRadiusMixture") {
    cols[16] = Fmt(function->GetParameter(0));
    cols[17] = Fmt(function->GetParError(0));
    cols[18] = Fmt(function->GetParameter(1));
    cols[19] = Fmt(function->GetParError(1));
  } else if (modelName == "normGammaMixtureQGaussian") {
    cols[6] = Fmt(function->GetParameter(0));
    cols[7] = Fmt(function->GetParError(0));
    cols[20] = Fmt(function->GetParameter(1));
    cols[21] = Fmt(function->GetParError(1));
  }

  out << centBin << '\t' << mtBin << '\t' << model;
  for (const std::string& col : cols) out << '\t' << col;
  out << '\n';
}

} // namespace

void rebuild_source_shape_fit_table_from_roots(
    TString scanDir = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS",
    TString tableFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/source_shape_fit_parameters.tsv")
{
  const TString header =
      "centBin\tmTBin\tmodel\t"
      "R_narrow\tR_narrow_err\tR_broad\tR_broad_err\tw_broad\tw_broad_err\t"
      "R\tR_err\talpha\talpha_err\t"
      "R_G\tR_G_err\tL_exp\tL_exp_err\tw_exp\tw_exp_err\t"
      "medianR\tmedianR_err\tsigmaLnR\tsigmaLnR_err\t"
      "nu\tnu_err";

  std::ofstream out(tableFileName.Data(), std::ios::out | std::ios::trunc);
  out << header << '\n';

  TSystemDirectory dir("scan_EPOS", scanDir);
  TList* files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "ERROR: cannot list directory " << scanDir << std::endl;
    return;
  }

  TIter next(files);
  TSystemFile* file = nullptr;
  int nFiles = 0;
  while ((file = dynamic_cast<TSystemFile*>(next()))) {
    TString fileName = file->GetName();
    if (file->IsDirectory()) continue;
    if (!fileName.EndsWith(".root")) continue;
    if (!fileName.BeginsWith("epos_primordial_pion_pair_source_fit_cent")) continue;

    const TString centBin = ExtractBetween(fileName, "_cent", "_mT");
    TString mtBin = ExtractAfter(fileName, "_mT");
    mtBin.ReplaceAll("p", ".");

    TFile inputFile(scanDir + "/" + fileName, "READ");
    if (inputFile.IsZombie()) continue;

    WriteRow(out, centBin, mtBin, "normDoubleGaussian",
             dynamic_cast<TF1*>(inputFile.Get("normDoubleGaussian")));
    WriteRow(out, centBin, mtBin, "normSymmetricLevy",
             dynamic_cast<TF1*>(inputFile.Get("normSymmetricLevy")));
    WriteRow(out, centBin, mtBin, "normGaussExpTail",
             dynamic_cast<TF1*>(inputFile.Get("normGaussExpTail")));
    WriteRow(out, centBin, mtBin, "normLogNormalRadiusMixture",
             dynamic_cast<TF1*>(inputFile.Get("normLogNormalRadiusMixture")));
    WriteRow(out, centBin, mtBin, "normGammaMixtureQGaussian",
             dynamic_cast<TF1*>(inputFile.Get("normGammaMixtureQGaussian")));
    ++nFiles;
  }

  out.close();
  std::cout << "Rebuilt " << tableFileName << " from " << nFiles << " ROOT files." << std::endl;
}
