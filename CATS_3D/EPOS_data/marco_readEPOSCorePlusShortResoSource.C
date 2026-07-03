#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

struct SourceCandidate {
  int charge = 0;
  int particleIndex = -1;
  int sourceIndex = -1;
  int pionId = 0;
  int sourceId = 0;
  int originType = 0;
  bool fromShortResonance = false;
  std::string sourceName = "direct";
  TLorentzVector pionP4;
  TLorentzVector pionX4;
  TLorentzVector sourceX4;
};

struct ShortResonanceInfo {
  bool accepted = false;
  double cTau = 0.;
  const char* name = "";
};

double EtaFromMomentum(const double px, const double py, const double pz)
{
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= std::abs(pz)) return (pz >= 0.) ? 1e9 : -1e9;
  return 0.5 * std::log((p + pz) / (p - pz));
}

double P4Energy(const double px, const double py, const double pz, const double mass)
{
  return std::sqrt(px * px + py * py + pz * pz + mass * mass);
}

double PairKStarFromPions(const TLorentzVector& p1, const TLorentzVector& p2)
{
  TLorentzVector p1Prf = p1;
  TLorentzVector p2Prf = p2;
  const TLorentzVector pair = p1 + p2;
  p1Prf.Boost(-pair.BoostVector());
  p2Prf.Boost(-pair.BoostVector());
  return 0.5 * (p1Prf.Vect() - p2Prf.Vect()).Mag();
}

double PairMtFromPions(const TLorentzVector& p1, const TLorentzVector& p2)
{
  constexpr double pionMass = 0.13957;
  const TVector3 kT = 0.5 * (p1.Vect() + p2.Vect());
  return std::sqrt(kT.Perp2() + pionMass * pionMass);
}

ShortResonanceInfo GetShortResonanceInfo(const int eposId, const double cTauMax)
{
  const int absId = std::abs(eposId);

  auto make = [&](const char* name, const double cTau) {
    ShortResonanceInfo info;
    info.accepted = (cTau > 0. && cTau < cTauMax);
    info.cTau = cTau;
    info.name = name;
    return info;
  };

  // Keep only resonance fathers with a physically short decay length.
  // c*tau values are hbar*c / width, rounded to the precision needed here.
  // Long-lived feed-down sources such as eta, omega, phi and weak-decay
  // hyperons are intentionally rejected by cTauMax=5 fm.
  if (absId == 111) return make("rho0", 1.3);
  if (absId == 121) return make("rho+", 1.3);
  if (absId == 131) return make("K*+", 4.2);
  if (absId == 231) return make("K*0", 4.2);
  if (absId == 112) return make("a0(980)", 2.5);
  if (absId == 122) return make("a0(980)+", 2.5);
  if (absId == 332) return make("f0(980)", 4.0);

  // Explicit long-lived mesons that may appear as pion fathers in EPOS.
  if (absId == 220) return make("eta", 1.5e5);
  if (absId == 221) return make("omega782", 23.4);
  if (absId == 331) return make("phi1020", 46.3);
  if (absId == 330) return make("etaPrime", 1.0e3);

  // Open charm mesons have weak/intermediate lifetimes, enormously longer
  // than femtoscopic short-lived resonance scales.
  if (absId == 140 || absId == 240 || absId == 141 || absId == 241) {
    return make("openCharm", 1.0e11);
  }

  // Delta and N* baryonic resonances can produce final pions through strong
  // decays. They are not meson feed-down, but they satisfy c*tau < 5 fm.
  if (absId == 1111 || absId == 1121 || absId == 1221 || absId == 2221) {
    return make("Delta(1232)", 1.7);
  }
  if (absId == 1112 || absId == 1113 || absId == 1114 ||
      absId == 2222 || absId == 2223 || absId == 2224 ||
      absId == 1122 || absId == 1123 || absId == 1124 ||
      absId == 1125 || absId == 1126 || absId == 1127 ||
      absId == 1128 || absId == 1222 || absId == 1223 ||
      absId == 1224 || absId == 1225 || absId == 1226 ||
      absId == 1227 || absId == 1228) {
    return make("NDeltaStar", 1.0);
  }

  return {};
}

bool IsSuppressedNoDecayId(const int id)
{
  static const std::set<int> suppressed = {
      110, 20, 2130, -2130, 2230, -2230, 1130, -1130,
      1330, -1330, 2330, -2330, 3331, -3331};
  return suppressed.count(id) > 0;
}

void FillSourcePair(const SourceCandidate& a,
                    const SourceCandidate& b,
                    TH1D* hRInv,
                    TH1D* hROut,
                    TH1D* hRSide,
                    TH1D* hRLong,
                    TH1D* hRAbsLCMS,
                    TH2D* hROutRSide,
                    TH2D* hRLongROut)
{
  const TLorentzVector pairP = a.pionP4 + b.pionP4;
  const TLorentzVector relX = a.sourceX4 - b.sourceX4;

  TLorentzVector relPrf = relX;
  relPrf.Boost(-pairP.BoostVector());
  hRInv->Fill(relPrf.Vect().Mag());

  TLorentzVector relLcms = relX;
  const double betaLong = pairP.Pz() / pairP.E();
  relLcms.Boost(0., 0., -betaLong);

  TVector3 outHat(pairP.Px(), pairP.Py(), 0.);
  if (outHat.Mag() < 1e-12) outHat.SetXYZ(1., 0., 0.);
  else outHat = outHat.Unit();
  const TVector3 sideHat(-outHat.Y(), outHat.X(), 0.);

  const TVector3 rLcms = relLcms.Vect();
  const double rOut = rLcms.Dot(outHat);
  const double rSide = rLcms.Dot(sideHat);
  const double rLong = rLcms.Z();

  hROut->Fill(rOut);
  hRSide->Fill(rSide);
  hRLong->Fill(rLong);
  hRAbsLCMS->Fill(std::sqrt(rOut * rOut + rSide * rSide + rLong * rLong));
  hROutRSide->Fill(rOut, rSide);
  hRLongROut->Fill(rLong, rOut);
}

void FillPairWithCoordinates(const SourceCandidate& a,
                             const SourceCandidate& b,
                             const TLorentzVector& xA,
                             const TLorentzVector& xB,
                             TH1D* hRInv,
                             TH1D* hROut = nullptr,
                             TH1D* hRSide = nullptr,
                             TH1D* hRLong = nullptr,
                             TH1D* hRAbsLCMS = nullptr,
                             TH2D* hROutRSide = nullptr,
                             TH2D* hRLongROut = nullptr)
{
  const TLorentzVector pairP = a.pionP4 + b.pionP4;
  const TLorentzVector relX = xA - xB;

  TLorentzVector relPrf = relX;
  relPrf.Boost(-pairP.BoostVector());
  if (hRInv) hRInv->Fill(relPrf.Vect().Mag());

  if (!hROut && !hRSide && !hRLong && !hRAbsLCMS &&
      !hROutRSide && !hRLongROut) {
    return;
  }

  TLorentzVector relLcms = relX;
  const double betaLong = pairP.Pz() / pairP.E();
  relLcms.Boost(0., 0., -betaLong);

  TVector3 outHat(pairP.Px(), pairP.Py(), 0.);
  if (outHat.Mag() < 1e-12) outHat.SetXYZ(1., 0., 0.);
  else outHat = outHat.Unit();
  const TVector3 sideHat(-outHat.Y(), outHat.X(), 0.);

  const TVector3 rLcms = relLcms.Vect();
  const double rOut = rLcms.Dot(outHat);
  const double rSide = rLcms.Dot(sideHat);
  const double rLong = rLcms.Z();

  if (hROut) hROut->Fill(rOut);
  if (hRSide) hRSide->Fill(rSide);
  if (hRLong) hRLong->Fill(rLong);
  if (hRAbsLCMS) hRAbsLCMS->Fill(std::sqrt(rOut * rOut + rSide * rSide + rLong * rLong));
  if (hROutRSide) hROutRSide->Fill(rOut, rSide);
  if (hRLongROut) hRLongROut->Fill(rLong, rOut);
}

TH1D* MakeNormalizedShape(const TH1D* sourceHist, const char* name)
{
  TH1D* normalized = dynamic_cast<TH1D*>(sourceHist->Clone(name));
  normalized->SetDirectory(nullptr);
  normalized->SetTitle(Form("%s normalized shape;%s;probability density",
                            sourceHist->GetTitle(),
                            sourceHist->GetXaxis()->GetTitle()));
  const double integral = normalized->Integral("width");
  if (integral > 0.) normalized->Scale(1. / integral);
  return normalized;
}

void marco_readEPOSCorePlusShortResoSource(
    TString inputFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/merged_z-EPOS_20250416.root",
    TString outputFileName = "epos_fluid_core_pion_short_reso_pair_source.root",
    TString treeName = "teposevent",
    Long64_t maxEvents = -1,
    int maxPairsPerChargePerEvent = 200000,
    int finalPionStatus = 0,
    double cTauMax = 5.0,
    double bMin = 0.0,
    double bMax = 3.5,
    double etaMax = 0.8,
    double pTMin = 0.2,
    double pTMax = 2.0,
    double kStarMin = 0.0,
    double kStarMax = 0.2,
    double mTMin = 0.0,
    double mTMax = 10.0,
    double maxAbsCoordinate = 1.0e6,
    bool requireFluidPion = true)
{
  TFile* inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open input file " << inputFileName << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(inputFile->Get(treeName));
  if (!tree) {
    std::cerr << "ERROR: cannot find tree " << treeName << " in " << inputFileName << std::endl;
    inputFile->Close();
    return;
  }

  constexpr int kMaxParticles = 300000;
  Int_t np = 0;
  Float_t bim = 0.;
  static Float_t px[kMaxParticles], py[kMaxParticles], pz[kMaxParticles], mass[kMaxParticles];
  static Float_t x[kMaxParticles], y[kMaxParticles], z[kMaxParticles], t[kMaxParticles];
  static Int_t id[kMaxParticles], ist[kMaxParticles], ity[kMaxParticles], ior[kMaxParticles], jor[kMaxParticles];

  tree->SetBranchAddress("np", &np);
  tree->SetBranchAddress("bim", &bim);
  tree->SetBranchAddress("px", px);
  tree->SetBranchAddress("py", py);
  tree->SetBranchAddress("pz", pz);
  tree->SetBranchAddress("e", mass);
  tree->SetBranchAddress("x", x);
  tree->SetBranchAddress("y", y);
  tree->SetBranchAddress("z", z);
  tree->SetBranchAddress("t", t);
  tree->SetBranchAddress("id", id);
  tree->SetBranchAddress("ist", ist);
  tree->SetBranchAddress("ity", ity);
  tree->SetBranchAddress("ior", ior);
  tree->SetBranchAddress("jor", jor);

  TH1D* hImpactParameter = new TH1D("hImpactParameter", "impact parameter;b (fm);events", 120, 0, 24);
  TH1D* hCandidateType = new TH1D("hCandidateType", "source candidate type;0 direct/string, 1 accepted resonance daughter;counts", 2, -0.5, 1.5);
  TH1D* hMotherIdAll = new TH1D("hMotherIdAll", "single-father pion mother id;mother EPOS id;counts", 20000, -10000, 10000);
  TH1D* hMotherIdAccepted = new TH1D("hMotherIdAccepted", "accepted rho/K* mother id;mother EPOS id;counts", 20000, -10000, 10000);
  TH1D* hSuppressedMotherId = new TH1D("hSuppressedMotherId", "suppressed nodecay id appearing as pion father;mother EPOS id;counts", 20000, -10000, 10000);
  TH1D* hPairKStar = new TH1D("hPairKStar", "selected measured pion pairs;k* (GeV/c);pairs", 200, 0, 1.0);
  TH1D* hPairMt = new TH1D("hPairMt", "selected measured pion pairs;m_{T} (GeV/c^{2});pairs", 200, 0, 2.0);

  auto makeH1 = [](const char* name, const char* title, const char* xTitle, double min, double max) {
    TH1D* h = new TH1D(name, Form("%s;%s;pairs", title, xTitle), 800, min, max);
    h->Sumw2();
    return h;
  };
  auto makeH2 = [](const char* name, const char* title) {
    TH2D* h = new TH2D(name, title, 300, -150, 150, 300, -150, 150);
    h->Sumw2();
    return h;
  };

  TH1D* hRInvAll = makeH1("hRInv_PRF_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in PRF", "r*_{inv} (fm)", 0, 400);
  TH1D* hROutAll = makeH1("hROut_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS", "r_{out} (fm)", -200, 200);
  TH1D* hRSideAll = makeH1("hRSide_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS", "r_{side} (fm)", -200, 200);
  TH1D* hRLongAll = makeH1("hRLong_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS", "r_{long} (fm)", -200, 200);
  TH1D* hRAbsAll = makeH1("hRAbs_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS", "|r|_{LCMS} (fm)", 0, 400);
  TH2D* hROutRSideAll = makeH2("hROutRSide_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS;r_{out} (fm);r_{side} (fm)");
  TH2D* hRLongROutAll = makeH2("hRLongROut_LCMS_primordial_plus_shortResoFather_all", "primordial plus short-resonance father source in LCMS;r_{long} (fm);r_{out} (fm)");

  TH1D* hRInvPP = makeH1("hRInv_PRF_primordial_primordial", "primordial pion-primordial pion source in PRF", "r*_{inv} (fm)", 0, 400);
  TH1D* hRInvPF = makeH1("hRInv_PRF_primordial_shortResoFather", "primordial pion-short-resonance father source in PRF", "r*_{inv} (fm)", 0, 400);
  TH1D* hRInvFF = makeH1("hRInv_PRF_shortResoFather_shortResoFather", "short-resonance father-short-resonance father source in PRF", "r*_{inv} (fm)", 0, 400);

  TH1D* hRInvPPPlusPDecayDaughter = makeH1("hRInv_PRF_primordial_primordial_plus_primordial_decayPion",
                                           "primordial-primordial plus primordial-decay pion source in PRF",
                                           "r*_{inv} (fm)", 0, 400);
  TH1D* hRInvFinalPionAll = makeH1("hRInv_PRF_finalPion_primordial_plus_shortResoDaughter_all",
                                   "final pion source: primordial plus short-resonance daughters in PRF",
                                   "r*_{inv} (fm)", 0, 400);
  TH1D* hROutFinalPionAll = makeH1("hROut_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                   "final pion source: primordial plus short-resonance daughters in LCMS",
                                   "r_{out} (fm)", -200, 200);
  TH1D* hRSideFinalPionAll = makeH1("hRSide_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                    "final pion source: primordial plus short-resonance daughters in LCMS",
                                    "r_{side} (fm)", -200, 200);
  TH1D* hRLongFinalPionAll = makeH1("hRLong_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                    "final pion source: primordial plus short-resonance daughters in LCMS",
                                    "r_{long} (fm)", -200, 200);
  TH1D* hRAbsFinalPionAll = makeH1("hRAbs_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                   "final pion source: primordial plus short-resonance daughters in LCMS",
                                   "|r|_{LCMS} (fm)", 0, 400);
  TH2D* hROutRSideFinalPionAll = makeH2("hROutRSide_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                        "final pion source: primordial plus short-resonance daughters in LCMS;r_{out} (fm);r_{side} (fm)");
  TH2D* hRLongROutFinalPionAll = makeH2("hRLongROut_LCMS_finalPion_primordial_plus_shortResoDaughter_all",
                                        "final pion source: primordial plus short-resonance daughters in LCMS;r_{long} (fm);r_{out} (fm)");

  TRandom3 rng(0);
  const Long64_t nEntries = tree->GetEntries();
  const Long64_t nEventsToRead = (maxEvents < 0) ? nEntries : std::min(maxEvents, nEntries);

  Long64_t acceptedEvents = 0;
  Long64_t rejectedLongOrUnknownFather = 0;
  Long64_t rejectedSuppressedFather = 0;
  Long64_t rejectedBadCoordinates = 0;
  Long64_t testedPairs = 0;
  Long64_t acceptedPairs = 0;

  auto goodCoordinates = [&](const int i) {
    return std::isfinite(x[i]) && std::isfinite(y[i]) &&
           std::isfinite(z[i]) && std::isfinite(t[i]) &&
           std::abs(x[i]) < maxAbsCoordinate &&
           std::abs(y[i]) < maxAbsCoordinate &&
           std::abs(z[i]) < maxAbsCoordinate &&
           std::abs(t[i]) < maxAbsCoordinate;
  };

  for (Long64_t iEvent = 0; iEvent < nEventsToRead; ++iEvent) {
    tree->GetEntry(iEvent);
    if (bim < bMin || bim >= bMax) continue;
    if (np > kMaxParticles) {
      std::cerr << "WARNING: event " << iEvent << " has np = " << np
                << ", larger than buffer. Skipping." << std::endl;
      continue;
    }

    std::vector<SourceCandidate> plus;
    std::vector<SourceCandidate> minus;

    for (int i = 0; i < np; ++i) {
      if (std::abs(id[i]) != 120) continue;
      if (finalPionStatus >= 0 && ist[i] != finalPionStatus) continue;
      if (requireFluidPion && ity[i] != 60) continue;

      const double pT = std::sqrt(px[i] * px[i] + py[i] * py[i]);
      if (pT < pTMin || pT >= pTMax) continue;

      const double eta = EtaFromMomentum(px[i], py[i], pz[i]);
      if (std::abs(eta) > etaMax) continue;

      const double pionEnergy = P4Energy(px[i], py[i], pz[i], mass[i]);

      SourceCandidate candidate;
      candidate.charge = (id[i] > 0) ? +1 : -1;
      candidate.particleIndex = i;
      candidate.pionId = id[i];
      candidate.originType = ity[i];
      candidate.pionP4.SetPxPyPzE(px[i], py[i], pz[i], pionEnergy);
      candidate.pionX4.SetXYZT(x[i], y[i], z[i], t[i]);

      const bool isSingleFatherDecay = (ior[i] >= 0 && jor[i] < 0);
      if (!isSingleFatherDecay) {
        if (!goodCoordinates(i)) {
          ++rejectedBadCoordinates;
          continue;
        }
        candidate.sourceIndex = i;
        candidate.sourceId = id[i];
        candidate.fromShortResonance = false;
        candidate.sourceName = "direct";
        candidate.sourceX4.SetXYZT(x[i], y[i], z[i], t[i]);
        hCandidateType->Fill(0);
      } else {
        const int father = ior[i];
        if (father < 0 || father >= np) continue;
        hMotherIdAll->Fill(id[father]);

        if (IsSuppressedNoDecayId(id[father])) {
          hSuppressedMotherId->Fill(id[father]);
          ++rejectedSuppressedFather;
          continue;
        }

        const ShortResonanceInfo reso = GetShortResonanceInfo(id[father], cTauMax);
        if (!reso.accepted) {
          ++rejectedLongOrUnknownFather;
          continue;
        }
        if (!goodCoordinates(i)) {
          ++rejectedBadCoordinates;
          continue;
        }
        if (!goodCoordinates(father)) {
          ++rejectedBadCoordinates;
          continue;
        }

        candidate.sourceIndex = father;
        candidate.sourceId = id[father];
        candidate.fromShortResonance = true;
        candidate.sourceName = reso.name;
        candidate.sourceX4.SetXYZT(x[father], y[father], z[father], t[father]);
        hMotherIdAccepted->Fill(id[father]);
        hCandidateType->Fill(1);
      }

      if (candidate.charge > 0) plus.push_back(candidate);
      else minus.push_back(candidate);
    }

    hImpactParameter->Fill(bim);
    ++acceptedEvents;

    auto processSameCharge = [&](const std::vector<SourceCandidate>& list) {
      const Long64_t n = static_cast<Long64_t>(list.size());
      if (n < 2) return;
      const Long64_t possiblePairs = n * (n - 1) / 2;
      const Long64_t trials = maxPairsPerChargePerEvent > 0
                                  ? std::min<Long64_t>(maxPairsPerChargePerEvent, possiblePairs)
                                  : possiblePairs;

      for (Long64_t trial = 0; trial < trials; ++trial) {
        Long64_t i = 0;
        Long64_t j = 0;
        if (maxPairsPerChargePerEvent > 0 && possiblePairs > maxPairsPerChargePerEvent) {
          i = rng.Integer(n);
          j = rng.Integer(n - 1);
          if (j >= i) ++j;
        } else {
          Long64_t count = 0;
          bool found = false;
          for (Long64_t a = 0; a < n && !found; ++a) {
            for (Long64_t b = a + 1; b < n; ++b) {
              if (count == trial) {
                i = a;
                j = b;
                found = true;
                break;
              }
              ++count;
            }
          }
        }

        ++testedPairs;
        const double kStar = PairKStarFromPions(list[i].pionP4, list[j].pionP4);
        if (kStar < kStarMin || kStar >= kStarMax) continue;
        const double mT = PairMtFromPions(list[i].pionP4, list[j].pionP4);
        if (mT < mTMin || mT >= mTMax) continue;

        hPairKStar->Fill(kStar);
        hPairMt->Fill(mT);
        FillSourcePair(list[i], list[j], hRInvAll, hROutAll, hRSideAll, hRLongAll,
                       hRAbsAll, hROutRSideAll, hRLongROutAll);

        const bool iReso = list[i].fromShortResonance;
        const bool jReso = list[j].fromShortResonance;
        if (!iReso && !jReso) {
          FillPairWithCoordinates(list[i], list[j], list[i].sourceX4, list[j].sourceX4, hRInvPP);
          FillPairWithCoordinates(list[i], list[j], list[i].pionX4, list[j].pionX4, hRInvPPPlusPDecayDaughter);
        }
        if (!iReso && jReso) {
          FillPairWithCoordinates(list[i], list[j], list[i].sourceX4, list[j].sourceX4, hRInvPF);
          FillPairWithCoordinates(list[i], list[j], list[i].pionX4, list[j].pionX4, hRInvPPPlusPDecayDaughter);
        }
        if (iReso && !jReso) {
          FillPairWithCoordinates(list[i], list[j], list[i].sourceX4, list[j].sourceX4, hRInvPF);
          FillPairWithCoordinates(list[i], list[j], list[i].pionX4, list[j].pionX4, hRInvPPPlusPDecayDaughter);
        }
        if (iReso && jReso) {
          FillPairWithCoordinates(list[i], list[j], list[i].sourceX4, list[j].sourceX4, hRInvFF);
        }
        FillPairWithCoordinates(list[i], list[j], list[i].pionX4, list[j].pionX4,
                                hRInvFinalPionAll, hROutFinalPionAll,
                                hRSideFinalPionAll, hRLongFinalPionAll,
                                hRAbsFinalPionAll, hROutRSideFinalPionAll,
                                hRLongROutFinalPionAll);

        ++acceptedPairs;
      }
    };

    processSameCharge(plus);
    processSameCharge(minus);

    if ((iEvent + 1) % 50 == 0) {
      std::cout << "Processed " << iEvent + 1 << " / " << nEventsToRead
                << " events, accepted pairs = " << acceptedPairs << std::endl;
    }
  }

  TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
  hImpactParameter->Write();
  hCandidateType->Write();
  hMotherIdAll->Write();
  hMotherIdAccepted->Write();
  hPairKStar->Write();
  hPairMt->Write();
  hRInvPP->Write();
  hRInvPF->Write();
  hRInvPPPlusPDecayDaughter->Write();
  TH1D* hRInvPPNorm = MakeNormalizedShape(hRInvPP, "hRInv_PRF_primordial_primordial_normShape");
  TH1D* hRInvPFNorm = MakeNormalizedShape(hRInvPF, "hRInv_PRF_primordial_shortResoFather_normShape");
  TH1D* hRInvPPPlusPDecayDaughterNorm =
      MakeNormalizedShape(hRInvPPPlusPDecayDaughter,
                          "hRInv_PRF_primordial_primordial_plus_primordial_decayPion_normShape");
  hRInvPPNorm->Write();
  hRInvPFNorm->Write();
  hRInvPPPlusPDecayDaughterNorm->Write();
  outputFile->Close();

  std::cout << "Saved " << outputFileName << std::endl;
  std::cout << "Require fluid/core pion ity==60: " << (requireFluidPion ? "yes" : "no") << std::endl;
  std::cout << "Accepted events: " << acceptedEvents << std::endl;
  std::cout << "Rejected long/unknown single-father pions: " << rejectedLongOrUnknownFather << std::endl;
  std::cout << "Rejected suppressed nodecay fathers: " << rejectedSuppressedFather << std::endl;
  std::cout << "Rejected bad-coordinate candidates: " << rejectedBadCoordinates << std::endl;
  std::cout << "Tested same-charge pairs: " << testedPairs << std::endl;
  std::cout << "Accepted pairs after k*/mT cuts: " << acceptedPairs << std::endl;

  inputFile->Close();
}
