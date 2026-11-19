#include "marco_fitEPOSPrimordialLevy.C"

#include <TFile.h>
#include <TKey.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <vector>

void run_rlong_source_shape_fits()
{
  const TString inputFileName =
      "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/epos_primordial_pion_pair_source_cent_mt.root";
  const TString outputDir =
      "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS";
  const TString tableFileName = outputDir + "/source_shape_fit_parameters3D_Rlong.tsv";

  TFile* inputFile = TFile::Open(inputFileName, "READ");
  if(!inputFile || inputFile->IsZombie()){
    std::cerr << "ERROR: cannot open " << inputFileName << std::endl;
    if(inputFile) inputFile->Close();
    return;
  }

  std::vector<TString> histNames;
  TIter nextKey(inputFile->GetListOfKeys());
  TKey* key = nullptr;
  while((key = dynamic_cast<TKey*>(nextKey()))){
    const TString keyName = key->GetName();
    if(keyName.BeginsWith("hRLong_LCMS_cent") && keyName.Contains("_mT")){
      histNames.push_back(keyName);
    }
  }
  inputFile->Close();

  std::sort(histNames.begin(), histNames.end());

  for(const TString& histName : histNames){
    TString suffix = histName;
    suffix.ReplaceAll("hRLong_LCMS_", "");
    const TString outputPrefix =
        outputDir + "/epos_primordial_pion_pair_3Dsource_Rlong_fit_" + suffix;

    std::cout << "\n=== Running " << histName << " ===\n";
    marco_fitEPOSPrimordialLevy(inputFileName,
                                histName,
                                -200.0,
                                200.0,
                                outputPrefix,
                                tableFileName);
  }

  std::cout << "\nFinished Rlong source-shape scan. Table: "
            << tableFileName << std::endl;
}
