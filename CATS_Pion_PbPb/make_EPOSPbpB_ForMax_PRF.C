#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace EPOSPbPbForMaxPRF {

constexpr int kMaxParticles = 250000;
constexpr double kPionMassMax = 0.25;

struct ParticleIdInfo {
  int pdgId = 0;
};

struct Particle {
  int entry = -1;
  int rawId = 0;
  int pdg = 0;
  TLorentzVector p4;
  TLorentzVector x4;
};

struct DecayedPion {
  int pionPDG = 0;
  int motherEntry = -1;
  int motherPDG = 0;
  TLorentzVector pionP4;
  TLorentzVector motherP4;
  TLorentzVector motherX4;
};

double GetPDGMass(const int pdg)
{
  const int apdg = std::abs(pdg);
  if (apdg == 211)  return 0.13957039;
  if (apdg == 111)  return 0.1349768;
  if (apdg == 2212) return 0.938272088;
  if (apdg == 2112) return 0.939565420;
  if (apdg == 321)  return 0.493677;
  if (apdg == 311)  return 0.497611;
  if (apdg == 113 || apdg == 213) return 0.77526;
  if (apdg == 223)  return 0.78265;
  if (apdg == 3222) return 1.18937;
  if (apdg == 3212) return 1.192642;
  if (apdg == 3112) return 1.197449;
  return -1.0;
}

ParticleIdInfo IdentifyParticleCode(const int rawId, const double mass)
{
  ParticleIdInfo out;
  const int sign = (rawId >= 0) ? 1 : -1;
  const int apdg = std::abs(rawId);
  auto set = [&](const int pdg) { out.pdgId = sign * pdg; };

  if (apdg == 120) set(211);
  else if (apdg == 121) set(213);
  else if (apdg == 131) set(323);
  else if (apdg == 231) set(313);
  else if (apdg == 221) set(223);
  else if (apdg == 111 && mass > 0.5) set(113);
  else if (apdg == 1111 || apdg == 1112 || apdg == 1113) set(2224);
  else if (apdg == 1121 || apdg == 1124 || apdg == 1126 || apdg == 1128 ||
           apdg == 12214) set(2214);
  else if (apdg == 1221 || apdg == 1224 || apdg == 1226 || apdg == 1228) set(2114);
  else if (apdg == 2221 || apdg == 2222 || apdg == 2223) set(1114);
  else if (apdg == 1114 && mass > 1.5) set(2224);
  else if (apdg == 2224 && mass > 1.5) set(1114);
  else if (apdg == 1122 || apdg == 1125 || apdg == 1127 ||
           apdg == 12212 || apdg == 22212 || apdg == 32212) set(12212);
  else if (apdg == 12112) set(12112);
  else if (apdg == 12224) set(12224);
  else set(apdg);
  return out;
}

bool GetDecayChannelWithChargedPion(const int motherPDG,
                                    const double motherMass,
                                    const bool rejectOmega,
                                    const bool includeStrangeBaryonSurrogates,
                                    int& pionPDG,
                                    int& otherPDG,
                                    TRandom3& rng)
{
  pionPDG = 0;
  otherPDG = 0;
  const int sign = (motherPDG >= 0) ? +1 : -1;
  const int apdg = std::abs(motherPDG);
  int pion = 0;
  int other = 0;

  if (apdg == 211 && motherMass < kPionMassMax) return false;
  if (rejectOmega && apdg == 223) return false;

  auto setDeltaPP = [&]() { pion = 211; other = 2212; };
  auto setDeltaP = [&]() {
    if (rng.Uniform() < 2.0/3.0) { pion = 111; other = 2212; }
    else { pion = 211; other = 2112; }
  };
  auto setDelta0 = [&]() {
    if (rng.Uniform() < 2.0/3.0) { pion = 111; other = 2112; }
    else { pion = -211; other = 2212; }
  };
  auto setDeltaM = [&]() { pion = -211; other = 2112; };
  auto setNstarP = [&]() {
    if (rng.Uniform() < 2.0/3.0) { pion = 211; other = 2112; }
    else { pion = 111; other = 2212; }
  };
  auto setNstar0 = [&]() {
    if (rng.Uniform() < 2.0/3.0) { pion = -211; other = 2212; }
    else { pion = 111; other = 2112; }
  };

  if (apdg == 2224) setDeltaPP();
  else if (apdg == 2214) setDeltaP();
  else if (apdg == 2114) setDelta0();
  else if (apdg == 1114) setDeltaM();
  else if (apdg == 12212 || apdg == 22212 || apdg == 32212) setNstarP();
  else if (apdg == 12112) setNstar0();
  else if (apdg == 12224) setDeltaPP();
  else if (apdg == 113) {
    if (rng.Uniform() < 0.5) { pion = 211; other = -211; }
    else { pion = -211; other = 211; }
  }
  else if (apdg == 213) { pion = 211; other = 111; }
  else if (apdg == 223) {
    if (rng.Uniform() < 0.5) { pion = 211; other = -211; }
    else { pion = -211; other = 211; }
  }
  else if (apdg == 313) {
    if (rng.Uniform() < 2.0/3.0) { pion = -211; other = 321; }
    else { pion = 111; other = 311; }
  }
  else if (apdg == 323) {
    if (rng.Uniform() < 2.0/3.0) { pion = 211; other = 311; }
    else { pion = 111; other = 321; }
  }
  else if (apdg == 20113) {
    if (rng.Uniform() < 0.5) { pion = -211; other = 213; }
    else { pion = 211; other = -213; }
  }
  else if (includeStrangeBaryonSurrogates && apdg == 2122) {
    if (rng.Uniform() < 0.5) { pion = -211; other = 3222; }
    else { pion = 211; other = 3112; }
  }
  else return false;

  if (sign < 0) {
    if (pion != 111 && pion != 113) pion = -pion;
    if (other != 111 && other != 113) other = -other;
  }
  if (std::abs(pion) != 211) return false;
  pionPDG = pion;
  otherPDG = other;
  return true;
}

bool DecayTwoBody(const Particle& parent,
                  const bool rejectOmega,
                  const bool includeStrangeBaryonSurrogates,
                  DecayedPion& outPion,
                  TRandom3& rng)
{
  const int motherPDG = IdentifyParticleCode(parent.rawId, parent.p4.M()).pdgId;
  if (motherPDG == 0) return false;

  int pionPDG = 0;
  int otherPDG = 0;
  if (!GetDecayChannelWithChargedPion(motherPDG, parent.p4.M(), rejectOmega,
                                      includeStrangeBaryonSurrogates,
                                      pionPDG, otherPDG, rng)) return false;

  int dau1PDG = pionPDG;
  int dau2PDG = otherPDG;
  if (rng.Uniform() < 0.5) std::swap(dau1PDG, dau2PDG);
  const double m1 = GetPDGMass(dau1PDG);
  const double m2 = GetPDGMass(dau2PDG);
  if (m1 <= 0. || m2 <= 0.) return false;

  const double mParent = parent.p4.M();
  if (mParent < m1 + m2) return false;
  const double term1 = mParent*mParent - (m1 + m2)*(m1 + m2);
  const double term2 = mParent*mParent - (m1 - m2)*(m1 - m2);
  const double pstar2 = term1 * term2 / (4.0 * mParent * mParent);
  if (pstar2 < 0.) return false;

  const double pstar = std::sqrt(pstar2);
  const double cosTheta = rng.Uniform(-1.0, 1.0);
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));
  const double phi = rng.Uniform(0.0, 2.0*TMath::Pi());
  TLorentzVector d1(pstar*sinTheta*std::cos(phi),
                    pstar*sinTheta*std::sin(phi),
                    pstar*cosTheta,
                    std::sqrt(m1*m1 + pstar*pstar));
  TLorentzVector d2(-d1.Px(), -d1.Py(), -d1.Pz(),
                    std::sqrt(m2*m2 + pstar*pstar));
  d1.Boost(parent.p4.BoostVector());
  d2.Boost(parent.p4.BoostVector());

  const bool firstIsPion = std::abs(dau1PDG) == 211;
  outPion.pionPDG = firstIsPion ? dau1PDG : dau2PDG;
  outPion.pionP4 = firstIsPion ? d1 : d2;
  outPion.motherEntry = parent.entry;
  outPion.motherPDG = motherPDG;
  outPion.motherP4 = parent.p4;
  outPion.motherX4 = parent.x4;
  return true;
}

double SafeAngle(const TVector3& a, const TVector3& b)
{
  const double ma = a.Mag();
  const double mb = b.Mag();
  if (ma <= 0. || mb <= 0.) return std::numeric_limits<double>::quiet_NaN();
  double c = a.Dot(b) / (ma * mb);
  c = std::max(-1.0, std::min(1.0, c));
  return std::acos(c);
}

bool FillTupleForPair(TNtuple* tuple,
                      TH1F* hMom1,
                      TH1F* hMom2,
                      TH2F* hMom12,
                      const TLorentzVector& source1P4,
                      const TLorentzVector& source1X4,
                      const TLorentzVector& daughter1P4,
                      const TLorentzVector& source2P4,
                      const TLorentzVector& source2X4,
                      const TLorentzVector& daughter2P4,
                      const double source1MassMeV,
                      const double source2MassMeV,
                      const double tau1,
                      const double tau2,
                      const double kStarMaxMeV)
{
  const TLorentzVector pair = daughter1P4 + daughter2P4;
  if (pair.E() <= 0. || pair.M2() <= 0.) return false;

  TLorentzVector d1PRF = daughter1P4;
  TLorentzVector source1PRF = source1P4;
  TLorentzVector source2PRF = source2P4;
  TLorentzVector rPRF = source2X4 - source1X4;
  const TVector3 boost = -pair.BoostVector();
  d1PRF.Boost(boost);
  source1PRF.Boost(boost);
  source2PRF.Boost(boost);
  rPRF.Boost(boost);

  const double kStarMeV = 1000.0 * d1PRF.P();
  if (kStarMeV > kStarMaxMeV) return false;

  const double p1 = 1000.0 * source1PRF.P();
  const double p2 = 1000.0 * source2PRF.P();
  const double m1 = source1MassMeV;
  const double m2 = source2MassMeV;
  if (p1 <= 0. || p2 <= 0. || m1 <= 0. || m2 <= 0.) return false;

  const TVector3 rVec = rPRF.Vect();
  const double angleRcP1 = SafeAngle(rVec, source1PRF.Vect());
  const double angleRcP2 = SafeAngle(rVec, source2PRF.Vect());
  const double angleP1P2 = SafeAngle(source1PRF.Vect(), source2PRF.Vect());
  if (!std::isfinite(angleRcP1) || !std::isfinite(angleRcP2) ||
      !std::isfinite(angleP1P2)) return false;

  tuple->Fill(kStarMeV, p1, p2, m1, m2, tau1, tau2,
              angleRcP1, angleRcP2, angleP1P2);
  hMom1->Fill(p1);
  hMom2->Fill(p2);
  hMom12->Fill(p1, p2);
  return true;
}

void AddCompatibilityObjects(TFile& file, TNtuple* tuple,
                             TH1F* hMom1, TH1F* hMom2, TH2F* hMom12)
{
  file.cd();
  auto* h1_AP_D = new TH1F("h1_AP_D", "h1_AP_D", 200, 0, 200);
  auto* hResoMass = new TH1F("hResoMass", "hResoMass", 240, 0, 2400);
  auto* hTau1 = new TH1F("h_TauDist_P1", "h_TauDist_P1", 100, 0, 10);
  auto* hTau2 = new TH1F("h_TauDist_P2", "h_TauDist_P2", 100, 0, 10);
  auto* hTau12 = new TH2F("h_TauDist_P1P2", "h_TauDist_P1P2", 100, 0, 10, 100, 0, 10);
  auto* hAngle1 = new TH1F("h_AngleRCP1_D", "h_AngleRCP1_D", 128, 0, TMath::Pi());
  auto* hAngle2 = new TH1F("h_AngleRCP2_D", "h_AngleRCP2_D", 128, 0, TMath::Pi());
  auto* hAngle12 = new TH1F("h_AngleP1P2_D", "h_AngleP1P2_D", 128, 0, TMath::Pi());
  auto* hAPLab = new TH1F("h_AP_LAB", "h_AP_LAB", 240, 0, 6);

  const std::vector<const char*> empty2DNames = {
      "h_AD_D", "h_AD_D1", "h_AD_D2", "h_AD_P1", "h_AD_P2",
      "h_AngleRCP1_P1P2_D", "h_AngleRCP1_RCP2_D",
      "h_AngleRCP2_P1P2_D", "h_AP_D", "h_AP_D1", "h_AP_D2",
      "h_AP_D_Large_AP_P1", "h_AP_D_Large_AP_P2",
      "h_AP_D_LargeAngles", "h_AP_LAB_AR_LAB",
      "h_AP_LAB_AR_LAB_Large_AP_P1", "h_AP_LAB_dR_LAB",
      "h_AP_P1", "h_AP_P1_AR_P1P2", "h_AP_P1_dR_P1P2",
      "h_AP_P1_LargeAngles", "h_AP_P2"
  };

  Float_t k_D = 0, P1 = 0, P2 = 0, M1 = 0, M2 = 0;
  Float_t Tau1 = 0, Tau2 = 0, AngleRcP1 = 0, AngleRcP2 = 0, AngleP1P2 = 0;
  tuple->SetBranchAddress("k_D", &k_D);
  tuple->SetBranchAddress("P1", &P1);
  tuple->SetBranchAddress("P2", &P2);
  tuple->SetBranchAddress("M1", &M1);
  tuple->SetBranchAddress("M2", &M2);
  tuple->SetBranchAddress("Tau1", &Tau1);
  tuple->SetBranchAddress("Tau2", &Tau2);
  tuple->SetBranchAddress("AngleRcP1", &AngleRcP1);
  tuple->SetBranchAddress("AngleRcP2", &AngleRcP2);
  tuple->SetBranchAddress("AngleP1P2", &AngleP1P2);
  for (Long64_t i = 0; i < tuple->GetEntries(); ++i) {
    tuple->GetEntry(i);
    h1_AP_D->Fill(k_D);
    hResoMass->Fill(M1);
    hResoMass->Fill(M2);
    hTau1->Fill(Tau1);
    hTau2->Fill(Tau2);
    hTau12->Fill(Tau1, Tau2);
    hAngle1->Fill(AngleRcP1);
    hAngle2->Fill(AngleRcP2);
    hAngle12->Fill(AngleP1P2);
  }
  hMom1->Write("h_MomDist_P1");
  hMom2->Write("h_MomDist_P2");
  hMom12->Write("h_MomDist_P1P2");
  h1_AP_D->Write();
  hResoMass->Write();
  hTau1->Write();
  hTau2->Write();
  hTau12->Write();
  hAngle1->Write();
  hAngle2->Write();
  hAngle12->Write();
  hAPLab->Write();
  for (const char* name : empty2DNames) {
    TH2F empty(name, name, 1, 0, 1, 1, 0, 1);
    empty.Write();
  }
  tuple->Write();
}

} // namespace EPOSPbPbForMaxPRF

void make_EPOSPbpB_ForMax_PRF(
    const char* inputName = "../CATS_3D/EPOS_data/merged_z-EPOS_20250416.root",
    const char* outPRName = "Files/EPOSPbPb_pi_piReso_noOmega_PRF.root",
    const char* outRRName = "Files/EPOSPbPb_piReso_piReso_noOmega_PRF.root",
    const double kStarMaxMeV = 200.0,
    const Long64_t maxEvents = -1,
    const bool rejectOmega = true,
    const bool includeStrangeBaryonSurrogates = false,
    const bool symmetrizeRR = false)
{
  using namespace EPOSPbPbForMaxPRF;

  TFile input(inputName, "READ");
  if (input.IsZombie()) {
    std::cerr << "ERROR: cannot open " << inputName << std::endl;
    return;
  }
  TTree* tree = dynamic_cast<TTree*>(input.Get("teposevent"));
  if (!tree) tree = dynamic_cast<TTree*>(input.Get("teposevent0"));
  if (!tree) {
    std::cerr << "ERROR: cannot find teposevent or teposevent0 in " << inputName << std::endl;
    return;
  }

  int np = 0;
  int nev = 0;
  static float px[kMaxParticles], py[kMaxParticles], pz[kMaxParticles], e[kMaxParticles];
  static float x[kMaxParticles], y[kMaxParticles], z[kMaxParticles], t[kMaxParticles];
  static int id[kMaxParticles], ist[kMaxParticles];

  tree->SetBranchAddress("np", &np);
  tree->SetBranchAddress("nev", &nev);
  tree->SetBranchAddress("px", px);
  tree->SetBranchAddress("py", py);
  tree->SetBranchAddress("pz", pz);
  tree->SetBranchAddress("e", e);
  tree->SetBranchAddress("x", x);
  tree->SetBranchAddress("y", y);
  tree->SetBranchAddress("z", z);
  tree->SetBranchAddress("t", t);
  tree->SetBranchAddress("id", id);
  tree->SetBranchAddress("ist", ist);

  TFile outPR(outPRName, "RECREATE");
  TFile outRR(outRRName, "RECREATE");
  outPR.cd();
  TNtuple* tuplePR = new TNtuple("InfoTuple_ClosePairs", "InfoTuple_ClosePairs",
                                 "k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
  outRR.cd();
  TNtuple* tupleRR = new TNtuple("InfoTuple_ClosePairs", "InfoTuple_ClosePairs",
                                 "k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
  auto* hPRP1 = new TH1F("hPRP1", "hPRP1", 250, 0, 10000);
  auto* hPRP2 = new TH1F("hPRP2", "hPRP2", 250, 0, 10000);
  auto* hPRP12 = new TH2F("hPRP12", "hPRP12", 125, 0, 10000, 125, 0, 10000);
  auto* hRRP1 = new TH1F("hRRP1", "hRRP1", 250, 0, 10000);
  auto* hRRP2 = new TH1F("hRRP2", "hRRP2", 250, 0, 10000);
  auto* hRRP12 = new TH2F("hRRP12", "hRRP12", 125, 0, 10000, 125, 0, 10000);

  TRandom3 rng(0);
  const Long64_t nEntries = (maxEvents > 0) ? std::min(maxEvents, tree->GetEntries())
                                            : tree->GetEntries();
  Long64_t nPR = 0;
  Long64_t nRR = 0;

  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry);
    if (np > kMaxParticles) {
      std::cerr << "WARNING: event " << iEntry << " np=" << np
                << " exceeds kMaxParticles=" << kMaxParticles << ", skipped\n";
      continue;
    }

    std::vector<Particle> primaryPlus;
    std::vector<Particle> primaryMinus;
    std::vector<DecayedPion> resoPlus;
    std::vector<DecayedPion> resoMinus;

    for (int i = 0; i < np; ++i) {
      Particle p;
      p.entry = i;
      p.rawId = id[i];
      p.p4.SetPxPyPzE(px[i], py[i], pz[i], e[i]);
      p.x4.SetXYZT(x[i], y[i], z[i], t[i]);
      p.pdg = IdentifyParticleCode(p.rawId, p.p4.M()).pdgId;

      if (std::abs(p.pdg) == 211 && p.p4.M() < kPionMassMax) {
        if (p.pdg > 0) primaryPlus.push_back(p);
        else primaryMinus.push_back(p);
      }

      DecayedPion decayed;
      if (DecayTwoBody(p, rejectOmega, includeStrangeBaryonSurrogates, decayed, rng)) {
        if (decayed.pionPDG > 0) resoPlus.push_back(decayed);
        else resoMinus.push_back(decayed);
      }
    }

    auto processPR = [&](const std::vector<Particle>& primaries,
                         const std::vector<DecayedPion>& resonances) {
      for (const auto& primary : primaries) {
        for (const auto& reso : resonances) {
          if (primary.entry == reso.motherEntry) continue;
          if (FillTupleForPair(tuplePR, hPRP1, hPRP2, hPRP12,
                               primary.p4, primary.x4, primary.p4,
                               reso.motherP4, reso.motherX4, reso.pionP4,
                               1000.0 * GetPDGMass(211),
                               1000.0 * reso.motherP4.M(),
                               0.0, 1.0, kStarMaxMeV)) ++nPR;
        }
      }
    };

    auto processRR = [&](const std::vector<DecayedPion>& resonances) {
      for (size_t i = 0; i < resonances.size(); ++i) {
        for (size_t j = i + 1; j < resonances.size(); ++j) {
          const auto& r1 = resonances[i];
          const auto& r2 = resonances[j];
          if (r1.motherEntry == r2.motherEntry) continue;
          if (FillTupleForPair(tupleRR, hRRP1, hRRP2, hRRP12,
                               r1.motherP4, r1.motherX4, r1.pionP4,
                               r2.motherP4, r2.motherX4, r2.pionP4,
                               1000.0 * r1.motherP4.M(),
                               1000.0 * r2.motherP4.M(),
                               1.0, 1.0, kStarMaxMeV)) ++nRR;
          if (symmetrizeRR &&
              FillTupleForPair(tupleRR, hRRP1, hRRP2, hRRP12,
                               r2.motherP4, r2.motherX4, r2.pionP4,
                               r1.motherP4, r1.motherX4, r1.pionP4,
                               1000.0 * r2.motherP4.M(),
                               1000.0 * r1.motherP4.M(),
                               1.0, 1.0, kStarMaxMeV)) ++nRR;
        }
      }
    };

    processPR(primaryPlus, resoPlus);
    processPR(primaryMinus, resoMinus);
    processRR(resoPlus);
    processRR(resoMinus);

    if ((iEntry + 1) % 50 == 0) {
      std::cout << "[INFO] processed " << (iEntry + 1) << "/" << nEntries
                << " events, PR=" << nPR << ", RR=" << nRR << std::endl;
    }
  }

  outPR.cd();
  AddCompatibilityObjects(outPR, tuplePR, hPRP1, hPRP2, hPRP12);
  outPR.Close();

  outRR.cd();
  AddCompatibilityObjects(outRR, tupleRR, hRRP1, hRRP2, hRRP12);
  outRR.Close();

  std::cout << "[INFO] Done. Wrote " << outPRName << " with " << nPR
            << " PR pairs and " << outRRName << " with " << nRR
            << " RR pairs." << std::endl;
}
