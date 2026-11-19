#include "marco_fitEPOSPrimordialLevy.C"

void run_remaining_source_shape_fits()
{
  const TString inputFileName =
      "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/epos_primordial_pion_pair_source_cent_mt.root";
  const TString outputDir =
      "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS";
  const TString tableFileName = outputDir + "/source_shape_fit_parameters.tsv";

  const TString centralityBins[] = {
      "20to30",
      "30to40",
      "40to50"
  };

  const TString mTBins[] = {
      "0p30to0p40",
      "0p40to0p50",
      "0p50to0p60",
      "0p60to0p70",
      "0p70to0p80"
  };

  for (const TString& centBin : centralityBins) {
    for (const TString& mTBin : mTBins) {
      const TString histName = "hRInv_PRF_cent" + centBin + "_mT" + mTBin;
      const TString outputPrefix =
          outputDir + "/epos_primordial_pion_pair_source_fit_cent" + centBin + "_mT" + mTBin;

      std::cout << "\n=== Running " << histName << " ===\n";
      marco_fitEPOSPrimordialLevy(inputFileName,
                                  histName,
                                  0.0,
                                  200.0,
                                  outputPrefix,
                                  tableFileName);
    }
  }
}
