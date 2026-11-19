#include "FitCATS3D.h"

#include "TROOT.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
const char* CentTag(int cent) {
  switch (cent) {
    case 0: return "0to10";
    case 1: return "10to30";
    case 2: return "30to50";
    // Temporary Pb-Pb proxy: no 50--80 lambda file in the selected folder.
    case 3: return "30to50";
    default: throw std::runtime_error("OO centrality index must be 0..3");
  }
}

FitCATS3DConfig MakeOOConfig(int cent, int qn, int phi) {
  if (cent != 1) throw std::runtime_error("this OO FT0A rough config is prepared for 10--30% (cent=1)");
  if (qn < 0 || qn > 2) throw std::runtime_error("OO qn index must be 0..2");
  if (phi < 0 || phi > 5) throw std::runtime_error("OO phi index must be 0..5");
  const char* centTag = CentTag(cent);

  FitCATS3DConfig cfg;
  cfg.InputFileName =
      "../../inputs/nuclearStructureDatafit/out_3DCF_OO_FT0CqVec_FT0APlane_cent10to30_mt0p244to0p714.root";
  cfg.SameHistoName = Form("H3_same_cent%d_mtMerged_qn%d_phi%d", cent, qn, phi);
  cfg.MixedHistoName = Form("H3_mixed_cent%d_mtMerged_phi%d", cent, phi);
  cfg.OutputFileName = Form(
      "../../fit/nuclearStructureDatafit/results/oo_ft0a_merged_mT/cent%d/qn%d/fit_phi%d.root", cent, qn, phi);

  cfg.UseSameMixed = true;
  cfg.Identical = true;
  cfg.WithCoulomb = true;
  cfg.UseLambdaFeedDown = true;
  cfg.LambdaFileName =
      "../../inputs/nuclearStructureDatafit/lambda_OOproxy_PbPb_cent10to30_mt0p2to0p7.root";
  cfg.LambdaHistoName =
      "hPurityKstar_OOproxy_mt0.244to0.714_cent10to30";
  cfg.weakDecayFraction = 0.12;

  cfg.RebinFactor = 1;
  cfg.fitQmaxMeV = 100.;
  cfg.normQminMeV = 150.;
  cfg.normQmaxMeV = 200.;
  cfg.UseSphericalNormalizationRegion = false;
  cfg.NormalizeDataInNormWindow = true;
  cfg.NormalizeModelInNormWindow = true;

  cfg.UseBackground = true;
  cfg.PrefitBackground = true;
  cfg.UsePmlBackgroundPrefit = true;
  cfg.UseDiagonalQuadraticBackground = true;
  cfg.UseHybridBackgroundRegion = true;
  cfg.BackgroundNonFemtoQminMeV = 100.;
  cfg.BackgroundNonFemtoQmaxMeV = 200.;
  cfg.sliceCutAbs = 0.04;
  cfg.FixBackground = true;
  cfg.FixNorm = true;
  cfg.UsePmlChi2 = true;

  cfg.MinuitMaxCalls = 20000;
  cfg.RunMinuitSimplexFirst = true;
  cfg.RunMinuitHesse = true;
  cfg.MinuitSimplexMaxCalls = 3000;
  cfg.MinuitTolerance = 0.1;
  cfg.MinuitRadiusStep = 0.3;
  cfg.RoutInit = 3.0;
  cfg.RsideInit = 3.0;
  cfg.RlongInit = 3.0;
  cfg.RadiusMinFm = 0.5;
  cfg.RadiusMaxFm = 12.0;
  cfg.ResonanceAlphaInit = 2.0;
  cfg.FixResonanceAlpha = true;

  cfg.UseFastKernel = true;
  cfg.FastKernelSubdivisions = 1;
  // O-O and Pb-Pb use the same 20-bin q grid over [-0.2, 0.2] GeV/c, so the
  // existing local 31^3 Coulomb wave-function kernel is directly compatible.
  cfg.UseExternalFastKernel = true;
  cfg.ExternalFastKernelByQ = true;
  cfg.ExternalFastKernelFileName =
      "../../inputs/OutputFiles/dataKernel_qData_r31_coulomb_byQ.root";
  cfg.ExternalFastKernelTreeName = "CATS3DFastKernel6DByQ";
  cfg.UseModelCache = false;
  cfg.NumGridBins = 31;
  cfg.gridMinFm = -50.0;
  cfg.gridMaxFm = 50.0;
  cfg.relativeMomentumFactor = 0.5;
  return cfg;
}
}  // namespace

int main(int argc, char** argv) {
  gROOT->SetBatch(true);
  if (argc != 4) {
    std::cerr << "Usage: ./fit_CATS3D_oo CENT_INDEX QN_INDEX PHI_INDEX\n";
    return 2;
  }
  try {
    const auto cfg = MakeOOConfig(std::stoi(argv[1]), std::stoi(argv[2]),
                                  std::stoi(argv[3]));
    std::filesystem::create_directories(
        std::filesystem::path(cfg.OutputFileName.Data()).parent_path());
    FitDataCATS3D(cfg);
  } catch (const std::exception& e) {
    std::cerr << "fit_CATS3D_oo: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
