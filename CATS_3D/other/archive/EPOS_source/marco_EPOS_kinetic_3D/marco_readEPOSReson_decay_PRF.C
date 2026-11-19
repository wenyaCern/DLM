#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct ParticleIdInfo;

double GetPDGMass(const int pdg);
ParticleIdInfo IdentifyParticleCode(const int rawId, const double mass);
bool GetDecayChannelWithChargedPion(const int motherPDG,
                                    const double motherMass,
                                    int& pionPDG,
                                    int& otherPDG,
                                    TRandom3& rng);

struct MomentumOSLPRF {
  double pOut = 0.;
  double pSide = 0.;
  double pLong = 0.;
  double ePRF = 0.;
  bool valid = false;
};

struct PRFKernel {
  double oS_dot_oD = 0.;
  double sS_dot_oD = 0.;
  double oS_dot_sD = 0.;
  double sS_dot_sD = 0.;
  double longScale = 0.;
  bool valid = false;
};

struct ParticleIdInfo {
  int rawId = 0;
  int pdgId = 0;
  bool isEposId = false;
  bool isPdgLike = false;
  const char* name = "";
};

struct DecayedPionWithMother {
  int pionPDG = 0;
  TLorentzVector pionP4;

  int motherPDG = 0;
  int motherEntry = -1;   // index of mother in an event，to avoid same mother pairing
  TLorentzVector motherP4;

  double motherX = 0.;
  double motherY = 0.;
  double motherZ = 0.;
  double motherT = 0.;
};

struct EposEventHeader {
  int eventIndex = 0;
  int nParticles = 0;
  double weight = 0.;
  double impactPar = 0.;
};

struct EposParticle {
  int idx = 0;
  int pdg = 0;

  double px = 0.;
  double py = 0.;
  double pz = 0.;
  double e  = 0.;
  double m  = 0.;

  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;

  int tag = 0;
};

struct DecayParticle {
  int pdg = 0;
  TLorentzVector p4;   // four-momentum in lab frame
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;
};

void FillResonTree(const EposEventHeader& evh,
                   TTree* treeReson,
                   int& out_eventIndex,
                   int& out_nPairs,
                   std::vector<int>* out_mother_pdg1,
                   std::vector<double>* out_motherInDaughterPRF_pOut1,
                   std::vector<double>* out_motherInDaughterPRF_pSide1,
                   std::vector<double>* out_motherInDaughterPRF_pLong1,
                   std::vector<double>* out_motherInDaughterPRF_e1,
                   std::vector<double>* out_m1,
                   std::vector<double>* out_motherInCorePRF_pOut1,
                   std::vector<double>* out_motherInCorePRF_pSide1,
                   std::vector<double>* out_motherInCorePRF_pLong1,
                   std::vector<double>* out_motherInCorePRF_e1,
                   const std::vector<MomentumOSLPRF>& boostedMothers1,
                   const std::vector<MomentumOSLPRF>& coreMothers1,
                   const std::vector<int>& pdgMothers1,
                   const std::vector<double>& motherMasses1,
                   std::vector<int>* out_mother_pdg2,
                   std::vector<double>* out_motherInDaughterPRF_pOut2,
                   std::vector<double>* out_motherInDaughterPRF_pSide2,
                   std::vector<double>* out_motherInDaughterPRF_pLong2,
                   std::vector<double>* out_motherInDaughterPRF_e2,
                   std::vector<double>* out_m2,
                   std::vector<double>* out_motherInCorePRF_pOut2,
                   std::vector<double>* out_motherInCorePRF_pSide2,
                   std::vector<double>* out_motherInCorePRF_pLong2,
                   std::vector<double>* out_motherInCorePRF_e2,
                   std::vector<double>* out_oS_dot_oD,
                   std::vector<double>* out_sS_dot_oD,
                   std::vector<double>* out_oS_dot_sD,
                   std::vector<double>* out_sS_dot_sD,
                   std::vector<double>* out_longScale,
                   const std::vector<MomentumOSLPRF>& boostedMothers2,
                   const std::vector<MomentumOSLPRF>& coreMothers2,
                   const std::vector<PRFKernel>& kernels,
                   const std::vector<int>& pdgMothers2,
                   const std::vector<double>& motherMasses2)
{
  out_mother_pdg1->clear();
  out_motherInDaughterPRF_pOut1->clear();
  out_motherInDaughterPRF_pSide1->clear();
  out_motherInDaughterPRF_pLong1->clear();
  out_motherInDaughterPRF_e1->clear();
  out_m1->clear();
  out_motherInCorePRF_pOut1->clear();
  out_motherInCorePRF_pSide1->clear();
  out_motherInCorePRF_pLong1->clear();
  out_motherInCorePRF_e1->clear();

  out_mother_pdg2->clear();
  out_motherInDaughterPRF_pOut2->clear();
  out_motherInDaughterPRF_pSide2->clear();
  out_motherInDaughterPRF_pLong2->clear();
  out_motherInDaughterPRF_e2->clear();
  out_m2->clear();
  out_motherInCorePRF_pOut2->clear();
  out_motherInCorePRF_pSide2->clear();
  out_motherInCorePRF_pLong2->clear();
  out_motherInCorePRF_e2->clear();

  out_oS_dot_oD->clear();
  out_sS_dot_oD->clear();
  out_oS_dot_sD->clear();
  out_sS_dot_sD->clear();
  out_longScale->clear();

  if (boostedMothers1.size() != pdgMothers1.size() ||
      boostedMothers1.size() != motherMasses1.size() ||
      boostedMothers1.size() != coreMothers1.size() ||
      boostedMothers1.size() != kernels.size()) {
    std::cerr << "FillResonTree ERROR: size mismatch in mother1 containers" << std::endl;
    return;
  }

  if (boostedMothers2.size() != pdgMothers2.size() ||
      boostedMothers2.size() != motherMasses2.size() ||
      boostedMothers2.size() != coreMothers2.size() ||
      boostedMothers2.size() != kernels.size()) {
    std::cerr << "FillResonTree ERROR: size mismatch in mother2 containers" << std::endl;
    return;
  }

  out_eventIndex = evh.eventIndex;
  out_nPairs = static_cast<int>(boostedMothers1.size());

  for (size_t i = 0; i < boostedMothers1.size(); ++i) {
    out_mother_pdg1->push_back(pdgMothers1[i]);
    out_motherInDaughterPRF_pOut1->push_back(boostedMothers1[i].pOut);
    out_motherInDaughterPRF_pSide1->push_back(boostedMothers1[i].pSide);
    out_motherInDaughterPRF_pLong1->push_back(boostedMothers1[i].pLong);
    out_motherInDaughterPRF_e1->push_back(boostedMothers1[i].ePRF);
    out_m1->push_back(motherMasses1[i]);
    out_motherInCorePRF_pOut1->push_back(coreMothers1[i].pOut);
    out_motherInCorePRF_pSide1->push_back(coreMothers1[i].pSide);
    out_motherInCorePRF_pLong1->push_back(coreMothers1[i].pLong);
    out_motherInCorePRF_e1->push_back(coreMothers1[i].ePRF);
  }

  for (size_t i = 0; i < boostedMothers2.size(); ++i) {
    out_mother_pdg2->push_back(pdgMothers2[i]);
    out_motherInDaughterPRF_pOut2->push_back(boostedMothers2[i].pOut);
    out_motherInDaughterPRF_pSide2->push_back(boostedMothers2[i].pSide);
    out_motherInDaughterPRF_pLong2->push_back(boostedMothers2[i].pLong);
    out_motherInDaughterPRF_e2->push_back(boostedMothers2[i].ePRF);
    out_m2->push_back(motherMasses2[i]);
    out_motherInCorePRF_pOut2->push_back(coreMothers2[i].pOut);
    out_motherInCorePRF_pSide2->push_back(coreMothers2[i].pSide);
    out_motherInCorePRF_pLong2->push_back(coreMothers2[i].pLong);
    out_motherInCorePRF_e2->push_back(coreMothers2[i].ePRF);
  }

  for (size_t i = 0; i < kernels.size(); ++i) {
    out_oS_dot_oD->push_back(kernels[i].oS_dot_oD);
    out_sS_dot_oD->push_back(kernels[i].sS_dot_oD);
    out_oS_dot_sD->push_back(kernels[i].oS_dot_sD);
    out_sS_dot_sD->push_back(kernels[i].sS_dot_sD);
    out_longScale->push_back(kernels[i].longScale);
  }

  treeReson->Fill();
}

bool ParseEventHeader(const std::string& line, EposEventHeader& h) {
  std::istringstream iss(line);
  if (!(iss >> h.eventIndex >> h.nParticles >> h.weight >> h.impactPar)) {
    return false;
  }
  return true;
}

bool ParseParticleLine(const std::string& line, EposParticle& p) {
  std::istringstream iss(line);
  if (!(iss >> p.idx >> p.pdg
            >> p.px >> p.py >> p.pz >> p.e >> p.m
            >> p.x >> p.y >> p.z >> p.t
            >> p.tag)) {
    return false;
  }
  return true;
}

// decay functions
bool DecayTwoBodyGivenDaughters(const TLorentzVector& parentP4,
                                const double parentX,
                                const double parentY,
                                const double parentZ,
                                const double parentT,
                                const int dau1PDG,
                                const int dau2PDG,
                                const double m1,
                                const double m2,
                                DecayParticle& d1,
                                DecayParticle& d2,
                                TRandom3& rng)
{
  const double M2 = parentP4.M2();
  if (M2 <= 0.) {
    std::cerr << "DecayTwoBodyGivenDaughters ERROR: parent M2 <= 0" << std::endl;
    return false;
  }

  const double M = std::sqrt(M2);

  if (M < m1 + m2) {
    std::cerr << "DecayTwoBodyGivenDaughters ERROR: threshold violated. "
              << "M = " << M << ", m1 + m2 = " << (m1 + m2) << std::endl;
    return false;
  }

  // two-body momentum in parent rest frame
  const double term1 = M*M - (m1 + m2)*(m1 + m2);
  const double term2 = M*M - (m1 - m2)*(m1 - m2);
  const double pstar2 = term1 * term2 / (4.0 * M * M);

  if (pstar2 < 0.) {
    std::cerr << "DecayTwoBodyGivenDaughters ERROR: pstar2 < 0" << std::endl;
    return false;
  }

  const double pstar = std::sqrt(pstar2);

  // isotropic direction in parent rest frame
  const double cosTheta = rng.Uniform(-1.0, 1.0);
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));
  const double phi = rng.Uniform(0.0, 2.0*TMath::Pi());

  const double px = pstar * sinTheta * std::cos(phi);
  const double py = pstar * sinTheta * std::sin(phi);
  const double pz = pstar * cosTheta;

  const double e1 = std::sqrt(m1*m1 + pstar*pstar);
  const double e2 = std::sqrt(m2*m2 + pstar*pstar);

  // rest-frame daughters
  TLorentzVector p4d1_rest(px,  py,  pz,  e1);
  TLorentzVector p4d2_rest(-px, -py, -pz, e2);

  // boost to lab
  TVector3 beta = parentP4.BoostVector();
  p4d1_rest.Boost(beta);
  p4d2_rest.Boost(beta);

  d1.pdg = dau1PDG;
  d1.p4 = p4d1_rest;
  d1.x = parentX;
  d1.y = parentY;
  d1.z = parentZ;
  d1.t = parentT;

  d2.pdg = dau2PDG;
  d2.p4 = p4d2_rest;
  d2.x = parentX;
  d2.y = parentY;
  d2.z = parentZ;
  d2.t = parentT;

  return true;
}

MomentumOSLPRF GetMotherMomentumInPairPRF(const TLorentzVector& daughter1,
                                            const TLorentzVector& daughter2,
                                            const TLorentzVector& mother)
{
  MomentumOSLPRF out;

  const TLorentzVector pair = daughter1 + daughter2;

  const double tPx = pair.Px();
  const double tPy = pair.Py();
  const double tPtSq = tPx * tPx + tPy * tPy;
  const double tPt = std::sqrt(tPtSq);

  if (tPt < 1e-12) {
    return out;
  }

  if (pair.E() <= 0. || pair.M2() <= 0.) {
    return out;
  }

  TLorentzVector motherPRF = mother;
  motherPRF.Boost(-pair.BoostVector());

  const double px = motherPRF.Px();
  const double py = motherPRF.Py();
  out.pOut  = (px * tPx + py * tPy) / tPt;
  out.pSide = (-px * tPy + py * tPx) / tPt;
  out.pLong = motherPRF.Pz();
  out.ePRF = motherPRF.E();

  out.valid = true;
  return out;
}

PRFKernel GetSourceToDaughterPRFKernel(const TLorentzVector& source1,
                                         const TLorentzVector& source2,
                                         const TLorentzVector& daughter1,
                                         const TLorentzVector& daughter2)
{
  PRFKernel out;

  const TLorentzVector pairS = source1 + source2;
  const TLorentzVector pairD = daughter1 + daughter2;

  const double ptS = std::sqrt(pairS.Px()*pairS.Px() + pairS.Py()*pairS.Py());
  const double ptD = std::sqrt(pairD.Px()*pairD.Px() + pairD.Py()*pairD.Py());
  if (ptS < 1e-12 || ptD < 1e-12) {
    return out;
  }
  if (pairS.E() <= 0. || pairD.E() <= 0.) {
    return out;
  }

  if (pairS.M2() <= 0. || pairD.M2() <= 0.) {
    return out;
  }

  TVector3 oS(pairS.Px()/ptS, pairS.Py()/ptS, 0.);
  TVector3 sS(-pairS.Py()/ptS, pairS.Px()/ptS, 0.);
  TVector3 oD(pairD.Px()/ptD, pairD.Py()/ptD, 0.);
  TVector3 sD(-pairD.Py()/ptD, pairD.Px()/ptD, 0.);

  out.oS_dot_oD = oS.Dot(oD);
  out.sS_dot_oD = sS.Dot(oD);
  out.oS_dot_sD = oS.Dot(sD);
  out.sS_dot_sD = sS.Dot(sD);

  out.longScale = 1.0;

  out.valid = true;
  return out;
}

bool GetTwoMothersInPairPRF(const TLorentzVector& daughter1,
                             const TLorentzVector& daughter2,
                             const TLorentzVector& mother1,
                             const TLorentzVector& mother2,
                             MomentumOSLPRF& mom1PRF,
                             MomentumOSLPRF& mom2PRF)
{
  mom1PRF = GetMotherMomentumInPairPRF(daughter1, daughter2, mother1);
  mom2PRF = GetMotherMomentumInPairPRF(daughter1, daughter2, mother2);

  if (!mom1PRF.valid || !mom2PRF.valid) {
    return false;
  }
  return true;
}

bool DecayTwoBody(const int motherPDG,
                  const TLorentzVector& parentP4,
                  const double parentX,
                  const double parentY,
                  const double parentZ,
                  const double parentT,
                  DecayParticle& d1,
                  DecayParticle& d2,
                  TRandom3& rng)
{
  const ParticleIdInfo motherId = IdentifyParticleCode(motherPDG, parentP4.M());
  if (motherId.pdgId == 0) {
    return false;
  }

  int pionPDG = 0;
  int otherPDG = 0;

  if (!GetDecayChannelWithChargedPion(motherId.pdgId, parentP4.M(),
                                      pionPDG, otherPDG, rng)) {
    return false;
  }

  int dau1PDG = pionPDG;
  int dau2PDG = otherPDG;
  if (rng.Uniform() < 0.5) {
    dau1PDG = otherPDG;
    dau2PDG = pionPDG;
  }

  const double m1 = GetPDGMass(dau1PDG);
  const double m2 = GetPDGMass(dau2PDG);

  if (m1 <= 0. || m2 <= 0.) {
    std::cerr << "DecayTwoBody ERROR: unknown daughter mass. "
              << "dau1PDG = " << dau1PDG
              << ", dau2PDG = " << dau2PDG << std::endl;
    return false;
  }

  return DecayTwoBodyGivenDaughters(parentP4,
                                    parentX, parentY, parentZ, parentT,
                                    dau1PDG, dau2PDG,
                                    m1, m2,
                                    d1, d2,
                                    rng);
}

void CheckDecayConservation(const TLorentzVector& parent,
                            const DecayParticle& d1,
                            const DecayParticle& d2)
{
  TLorentzVector sum = d1.p4 + d2.p4;

  std::cout << "---- Conservation check ----" << std::endl;
  std::cout << "Parent : "
            << parent.Px() << "  "
            << parent.Py() << "  "
            << parent.Pz() << "  "
            << parent.E()  << std::endl;

  std::cout << "Dausum : "
            << sum.Px() << "  "
            << sum.Py() << "  "
            << sum.Pz() << "  "
            << sum.E()  << std::endl;

  std::cout << "Diff   : "
            << parent.Px() - sum.Px() << "  "
            << parent.Py() - sum.Py() << "  "
            << parent.Pz() - sum.Pz() << "  "
            << parent.E()  - sum.E()  << std::endl;
}

void ProcessEvent(const EposEventHeader& evh,
                  const std::vector<EposParticle>& particles,
                  TTree* treeReson,
                  int& out_eventIndex,
                  int& out_nPairs,
                  std::vector<int>* out_mother_pdg1,
                  std::vector<double>* out_motherInDaughterPRF_pOut1,
                  std::vector<double>* out_motherInDaughterPRF_pSide1,
                  std::vector<double>* out_motherInDaughterPRF_pLong1,
                  std::vector<double>* out_motherInDaughterPRF_e1,
                  std::vector<double>* out_m1,
                  std::vector<double>* out_motherInCorePRF_pOut1,
                  std::vector<double>* out_motherInCorePRF_pSide1,
                  std::vector<double>* out_motherInCorePRF_pLong1,
                  std::vector<double>* out_motherInCorePRF_e1,
                  std::vector<int>* out_mother_pdg2,
                  std::vector<double>* out_motherInDaughterPRF_pOut2,
                  std::vector<double>* out_motherInDaughterPRF_pSide2,
                  std::vector<double>* out_motherInDaughterPRF_pLong2,
                  std::vector<double>* out_motherInDaughterPRF_e2,
                  std::vector<double>* out_m2,
                  std::vector<double>* out_motherInCorePRF_pOut2,
                  std::vector<double>* out_motherInCorePRF_pSide2,
                  std::vector<double>* out_motherInCorePRF_pLong2,
                  std::vector<double>* out_motherInCorePRF_e2,
                  std::vector<double>* out_oS_dot_oD,
                  std::vector<double>* out_sS_dot_oD,
                  std::vector<double>* out_oS_dot_sD,
                  std::vector<double>* out_sS_dot_sD,
                  std::vector<double>* out_longScale)
{
  static TRandom3 rng(0);
  const bool printPairs = false;

  std::vector<DecayedPionWithMother> piPlusList;
  std::vector<DecayedPionWithMother> piMinusList;

  // 1) decay all mothers, keep charged pions
  for (size_t iMom = 0; iMom < particles.size(); ++iMom) {
    const auto& p = particles[iMom];

    TLorentzVector parentP4(p.px, p.py, p.pz, p.e);

    DecayParticle d1;
    DecayParticle d2;

    bool ok = DecayTwoBody(p.pdg, parentP4, p.x, p.y, p.z, p.t, d1, d2, rng);
    if (!ok) {
      continue;
    }

    if (d1.pdg == 211 || d1.pdg == -211) {
      DecayedPionWithMother onePion;
      onePion.pionPDG = d1.pdg;
      onePion.pionP4 = d1.p4;
      onePion.motherPDG = p.pdg;
      onePion.motherEntry = static_cast<int>(iMom);
      onePion.motherP4 = parentP4;
      onePion.motherX = p.x;
      onePion.motherY = p.y;
      onePion.motherZ = p.z;
      onePion.motherT = p.t;

      if (d1.pdg == 211) piPlusList.push_back(onePion);
      else piMinusList.push_back(onePion);
    }

    if (d2.pdg == 211 || d2.pdg == -211) {
      DecayedPionWithMother onePion;
      onePion.pionPDG = d2.pdg;
      onePion.pionP4 = d2.p4;
      onePion.motherPDG = p.pdg;
      onePion.motherEntry = static_cast<int>(iMom);
      onePion.motherP4 = parentP4;
      onePion.motherX = p.x;
      onePion.motherY = p.y;
      onePion.motherZ = p.z;
      onePion.motherT = p.t;

      if (d2.pdg == 211) piPlusList.push_back(onePion);
      else piMinusList.push_back(onePion);
    }
  }

  // event-level containers for resonance tree
  std::vector<MomentumOSLPRF> boostedMothers1;
  std::vector<MomentumOSLPRF> boostedMothers2;
  std::vector<MomentumOSLPRF> coreMothers1;
  std::vector<MomentumOSLPRF> coreMothers2;
  std::vector<PRFKernel> kernels;
  std::vector<int> pdgMothers1;
  std::vector<int> pdgMothers2;
  std::vector<double> motherMasses1;
  std::vector<double> motherMasses2;

  auto ProcessSameSignPairs = [&](const std::vector<DecayedPionWithMother>& pionList,
                                  const char* pairName)
  {
    for (size_t i = 0; i < pionList.size(); ++i) {
      for (size_t j = i + 1; j < pionList.size(); ++j) {
        const auto& pi1 = pionList[i];
        const auto& pi2 = pionList[j];

        if (pi1.motherEntry == pi2.motherEntry) {
          continue;
        }

        const TLorentzVector& pion1 = pi1.pionP4;
        const TLorentzVector& pion2 = pi2.pionP4;
        const TLorentzVector& mother1 = pi1.motherP4;
        const TLorentzVector& mother2 = pi2.motherP4;

        MomentumOSLPRF mom1PRF = GetMotherMomentumInPairPRF(pion1, pion2, mother1);
        MomentumOSLPRF mom2PRF = GetMotherMomentumInPairPRF(pion1, pion2, mother2);
        MomentumOSLPRF mom1CorePRF = GetMotherMomentumInPairPRF(mother1, mother2, mother1);
        MomentumOSLPRF mom2CorePRF = GetMotherMomentumInPairPRF(mother1, mother2, mother2);
        PRFKernel kernel = GetSourceToDaughterPRFKernel(mother1, mother2, pion1, pion2);

        if (!mom1PRF.valid || !mom2PRF.valid ||
            !mom1CorePRF.valid || !mom2CorePRF.valid ||
            !kernel.valid) {
          continue;
        }

        boostedMothers1.push_back(mom1PRF);
        boostedMothers2.push_back(mom2PRF);
        coreMothers1.push_back(mom1CorePRF);
        coreMothers2.push_back(mom2CorePRF);
        kernels.push_back(kernel);

        pdgMothers1.push_back(pi1.motherPDG);
        pdgMothers2.push_back(pi2.motherPDG);

        motherMasses1.push_back(pi1.motherP4.M());
        motherMasses2.push_back(pi2.motherP4.M());

        if (printPairs) {
          std::cout << "[PAIR " << pairName << "] "
                    << "event=" << evh.eventIndex
                    << " motherPDG=(" << pi1.motherPDG << "," << pi2.motherPDG << ") "
                    << " pionPDG=(" << pi1.pionPDG << "," << pi2.pionPDG << ") "
                    << "\n   mother1 OSL-PRF = ("
                    << mom1PRF.pOut << ", "
                    << mom1PRF.pSide << ", "
                    << mom1PRF.pLong << "), E=" << mom1PRF.ePRF
                    << "\n   mother2 OSL-PRF = ("
                    << mom2PRF.pOut << ", "
                    << mom2PRF.pSide << ", "
                    << mom2PRF.pLong << "), E=" << mom2PRF.ePRF
                    << std::endl;
        }
      }
    }
  };

  ProcessSameSignPairs(piPlusList, "pi+pi+");
  ProcessSameSignPairs(piMinusList, "pi-pi-");

  FillResonTree(evh,
                treeReson,
                out_eventIndex,
                out_nPairs,
                out_mother_pdg1,
                out_motherInDaughterPRF_pOut1,
                out_motherInDaughterPRF_pSide1,
                out_motherInDaughterPRF_pLong1,
                out_motherInDaughterPRF_e1,
                out_m1,
                out_motherInCorePRF_pOut1,
                out_motherInCorePRF_pSide1,
                out_motherInCorePRF_pLong1,
                out_motherInCorePRF_e1,
                boostedMothers1,
                coreMothers1,
                pdgMothers1,
                motherMasses1,
                out_mother_pdg2,
                out_motherInDaughterPRF_pOut2,
                out_motherInDaughterPRF_pSide2,
                out_motherInDaughterPRF_pLong2,
                out_motherInDaughterPRF_e2,
                out_m2,
                out_motherInCorePRF_pOut2,
                out_motherInCorePRF_pSide2,
                out_motherInCorePRF_pLong2,
                out_motherInCorePRF_e2,
                out_oS_dot_oD,
                out_sS_dot_oD,
                out_oS_dot_sD,
                out_sS_dot_sD,
                out_longScale,
                boostedMothers2,
                coreMothers2,
                kernels,
                pdgMothers2,
                motherMasses2);
}

void marco_readEPOSReson_decay_PRF(const char* inputFileName = "epos_input.dat",
                                   const char* outputResonanceName = "boosted_resonance_prf_output.root")
{
  ifstream fin(inputFileName);
  if (!fin.is_open()) {
    cerr << "ERROR: cannot open input file: " << inputFileName << endl;
    return;
  }

  // -------------------------
  // Event-level branches
  // -------------------------
  int out_eventIndex = 0;
  int out_nPairs = 0;

  TFile* foutReson = new TFile(outputResonanceName, "RECREATE");
  if (!foutReson || foutReson->IsZombie()) {
    cerr << "ERROR: cannot create output ROOT file: " << outputResonanceName << endl;
    fin.close();
    return;
  }

  TTree* treeReson = new TTree("EPOS_Reson_PRF",
                               "EPOS resonance mothers in daughter-pair-defined OSL-PRF");
  
  treeReson->Branch("eventIndex", &out_eventIndex, "eventIndex/I");
  treeReson->Branch("nPairs", &out_nPairs, "nPairs/I");
  
  // -------------------------
  // mother 1 branches
  // -------------------------
  std::vector<int>* out_mother_pdg1       = new std::vector<int>();
  std::vector<double>* out_motherInDaughterPRF_pOut1 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_pSide1 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_pLong1 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_e1    = new std::vector<double>();
  std::vector<double>* out_m1             = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pOut1     = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pSide1    = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pLong1    = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_e1        = new std::vector<double>();
  
  treeReson->Branch("mother_pdg1",       &out_mother_pdg1);
  treeReson->Branch("motherInDaughterPRF_pOut1",    &out_motherInDaughterPRF_pOut1);
  treeReson->Branch("motherInDaughterPRF_pSide1",   &out_motherInDaughterPRF_pSide1);
  treeReson->Branch("motherInDaughterPRF_pLong1",   &out_motherInDaughterPRF_pLong1);
  treeReson->Branch("motherInDaughterPRF_e1",       &out_motherInDaughterPRF_e1);
  treeReson->Branch("m1",             &out_m1);
  treeReson->Branch("motherInCorePRF_pOut1",     &out_motherInCorePRF_pOut1);
  treeReson->Branch("motherInCorePRF_pSide1",    &out_motherInCorePRF_pSide1);
  treeReson->Branch("motherInCorePRF_pLong1",    &out_motherInCorePRF_pLong1);
  treeReson->Branch("motherInCorePRF_e1",        &out_motherInCorePRF_e1);
  
  // -------------------------
  // mother 2 branches
  // -------------------------
  std::vector<int>* out_mother_pdg2       = new std::vector<int>();
  std::vector<double>* out_motherInDaughterPRF_pOut2 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_pSide2 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_pLong2 = new std::vector<double>();
  std::vector<double>* out_motherInDaughterPRF_e2    = new std::vector<double>();
  std::vector<double>* out_m2             = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pOut2     = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pSide2    = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_pLong2    = new std::vector<double>();
  std::vector<double>* out_motherInCorePRF_e2        = new std::vector<double>();
  
  treeReson->Branch("mother_pdg2",       &out_mother_pdg2);
  treeReson->Branch("motherInDaughterPRF_pOut2",    &out_motherInDaughterPRF_pOut2);
  treeReson->Branch("motherInDaughterPRF_pSide2",   &out_motherInDaughterPRF_pSide2);
  treeReson->Branch("motherInDaughterPRF_pLong2",   &out_motherInDaughterPRF_pLong2);
  treeReson->Branch("motherInDaughterPRF_e2",       &out_motherInDaughterPRF_e2);
  treeReson->Branch("m2",             &out_m2);
  treeReson->Branch("motherInCorePRF_pOut2",     &out_motherInCorePRF_pOut2);
  treeReson->Branch("motherInCorePRF_pSide2",    &out_motherInCorePRF_pSide2);
  treeReson->Branch("motherInCorePRF_pLong2",    &out_motherInCorePRF_pLong2);
  treeReson->Branch("motherInCorePRF_e2",        &out_motherInCorePRF_e2);

  std::vector<double>* out_oS_dot_oD     = new std::vector<double>();
  std::vector<double>* out_sS_dot_oD     = new std::vector<double>();
  std::vector<double>* out_oS_dot_sD     = new std::vector<double>();
  std::vector<double>* out_sS_dot_sD     = new std::vector<double>();
  std::vector<double>* out_longScale     = new std::vector<double>();

  treeReson->Branch("oS_dot_oD",    &out_oS_dot_oD);
  treeReson->Branch("sS_dot_oD",    &out_sS_dot_oD);
  treeReson->Branch("oS_dot_sD",    &out_oS_dot_sD);
  treeReson->Branch("sS_dot_sD",    &out_sS_dot_sD);
  treeReson->Branch("longScale",    &out_longScale);

  string line;

  // -------------------------
  // 读前 3 行全局 header
  // -------------------------
  string globalHeader1, globalHeader2, globalHeader3;

  if (!std::getline(fin, globalHeader1)) {
    cerr << "ERROR: failed to read first global header line." << endl;
    foutReson->Close();
    fin.close();
    return;
  }

  if (!std::getline(fin, globalHeader2)) {
    cerr << "ERROR: failed to read second global header line." << endl;
    foutReson->Close();
    fin.close();
    return;
  }

  if (!std::getline(fin, globalHeader3)) {
    cerr << "ERROR: failed to read third global header line." << endl;
    foutReson->Close();
    fin.close();
    return;
  }

  cout << "[INFO] Global header 1: " << globalHeader1 << endl;
  cout << "[INFO] Global header 2: " << globalHeader2 << endl;
  cout << "[INFO] Global header 3: " << globalHeader3 << endl;

  // -------------------------
  // 逐 event 读取
  // -------------------------
  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    EposEventHeader evh;
    if (!ParseEventHeader(line, evh)) {
      cerr << "WARNING: failed to parse event header line: " << line << endl;
      continue;
    }

    vector<EposParticle> particles;
    particles.reserve(evh.nParticles);

    for (int i = 0; i < evh.nParticles; ++i) {
      if (!std::getline(fin, line)) {
        cerr << "ERROR: unexpected end of file while reading particles for event "
             << evh.eventIndex << endl;
        break;
      }

      if (line.empty()) {
        --i;
        continue;
      }

      EposParticle p;
      if (!ParseParticleLine(line, p)) {
        cerr << "WARNING: failed to parse particle line in event "
             << evh.eventIndex << ": " << line << endl;
        continue;
      }

      particles.push_back(p);
    }

    ProcessEvent(evh,
             particles,
             treeReson,
             out_eventIndex,
             out_nPairs,
             out_mother_pdg1,
             out_motherInDaughterPRF_pOut1,
             out_motherInDaughterPRF_pSide1,
             out_motherInDaughterPRF_pLong1,
             out_motherInDaughterPRF_e1,
             out_m1,
             out_motherInCorePRF_pOut1,
             out_motherInCorePRF_pSide1,
             out_motherInCorePRF_pLong1,
             out_motherInCorePRF_e1,
             out_mother_pdg2,
             out_motherInDaughterPRF_pOut2,
             out_motherInDaughterPRF_pSide2,
             out_motherInDaughterPRF_pLong2,
             out_motherInDaughterPRF_e2,
             out_m2,
             out_motherInCorePRF_pOut2,
             out_motherInCorePRF_pSide2,
             out_motherInCorePRF_pLong2,
             out_motherInCorePRF_e2,
             out_oS_dot_oD,
             out_sS_dot_oD,
             out_oS_dot_sD,
             out_sS_dot_sD,
             out_longScale);
  }

  foutReson->cd();
  treeReson->Write();
  foutReson->Close();

  fin.close();
  
  delete out_mother_pdg1;
  delete out_motherInDaughterPRF_pOut1;
  delete out_motherInDaughterPRF_pSide1;
  delete out_motherInDaughterPRF_pLong1;
  delete out_motherInDaughterPRF_e1;
  delete out_m1;
  delete out_motherInCorePRF_pOut1;
  delete out_motherInCorePRF_pSide1;
  delete out_motherInCorePRF_pLong1;
  delete out_motherInCorePRF_e1;
  
  delete out_mother_pdg2;
  delete out_motherInDaughterPRF_pOut2;
  delete out_motherInDaughterPRF_pSide2;
  delete out_motherInDaughterPRF_pLong2;
  delete out_motherInDaughterPRF_e2;
  delete out_m2;
  delete out_motherInCorePRF_pOut2;
  delete out_motherInCorePRF_pSide2;
  delete out_motherInCorePRF_pLong2;
  delete out_motherInCorePRF_e2;
  delete out_oS_dot_oD;
  delete out_sS_dot_oD;
  delete out_oS_dot_sD;
  delete out_sS_dot_sD;
  delete out_longScale;
  
  std::cout << "[INFO] Done. Output written to " << outputResonanceName << std::endl;
}

// ---------------------------------------------------------------------------
// EPOS internal particle IDs used in the OSC-style input, and their closest
// standard PDG IDs. The EPOS table can be inspected in EPOS4 src/KWt/idt.dt
// or in the EPOS4 particle-id documentation.
//
// Basic mesons:
//   EPOS ID   EPOS name       closest PDG ID      PDG name
//      110    PI0                  111            pi0
//      120    PI+                  211            pi+
//     -120    PI-                 -211            pi-
//      111    RHO0                 113            rho0
//      121    RHO+                 213            rho+
//     -121    RHO-                -213            rho-
//      221    OMEG                 223            omega
//      131    K*+                  323            K*(892)+
//     -131    K*-                 -323            K*(892)-
//      231    K*0                  313            K*(892)0
//     -231    AK*0                -313            anti-K*(892)0
//
// Delta(1232) baryons:
//   EPOS ID   EPOS name       closest PDG ID      PDG name
//     1111    DL++                2224            Delta++
//     1121    DL+                 2214            Delta+
//     1221    DL0                 2114            Delta0
//     2221    DL-                 1114            Delta-
//
// Higher EPOS baryon resonances currently mapped to N pi or Delta pi
// two-body charged-pion surrogates in this macro:
//   EPOS ID   EPOS name       approximate treatment
//     1112    DL++(1620)      Delta++ -> p pi+
//     1113    DL++(1700)      Delta++ -> p pi+
//     1114    DL++(1925)      Delta++ -> p pi+   (if EPOS mass is high)
//     2222    DL-(1620)       Delta-  -> n pi-
//     2223    DL-(1700)       Delta-  -> n pi-
//     2224    DL-(1925)       Delta-  -> n pi-   (if EPOS mass is high)
//     1122    N*+(1440)       N*+ -> n pi+ / p pi0
//     1125    N*+(1665)       N*+ -> n pi+ / p pi0
//     1127    N*+(1710)       N*+ -> n pi+ / p pi0
//     1124    DL+(1620)       Delta+ -> n pi+ / p pi0
//     1126    DL+(1700)       Delta+ -> n pi+ / p pi0
//     1128    DL+(1925)       Delta+ -> n pi+ / p pi0
//     1224    DL0(1620)       Delta0 -> p pi- / n pi0
//     1226    DL0(1700)       Delta0 -> p pi- / n pi0
//     1228    DL0(1925)       Delta0 -> p pi- / n pi0
//
// Standard PDG-like IDs also present in the input and handled here:
//   PDG ID    PDG name        approximate treatment
//      113    rho0            pi+ pi-
//      213    rho+            pi+ pi0
//      223    omega           pi+ pi- two-body surrogate
//      313    K*(892)0        K+ pi- / K0 pi0
//      323    K*(892)+        K0 pi+ / K+ pi0
//     1114    Delta-          n pi-   (if mass is Delta(1232)-like)
//     2114    Delta0          p pi- / n pi0
//     2214    Delta+          n pi+ / p pi0
//     2224    Delta++         p pi+   (if mass is Delta(1232)-like)
//    12112    N(1440)0        p pi- / n pi0
//    12212    N*+ surrogate   n pi+ / p pi0
//    12214    Delta+ surrogate
//    12224    Delta++ surrogate
//    20113    a1(1260)0       rho+ pi- / rho- pi+ two-body surrogate
//    22212    N*+ surrogate   n pi+ / p pi0
//    32212    N*+ surrogate   n pi+ / p pi0
// ---------------------------------------------------------------------------

ParticleIdInfo IdentifyParticleCode(const int rawId, const double mass)
{
  ParticleIdInfo out;
  out.rawId = rawId;

  const int sign = (rawId >= 0) ? 1 : -1;
  const int apdg = std::abs(rawId);

  auto setEpos = [&](const int pdgId, const char* name) {
    out.pdgId = sign * pdgId;
    out.isEposId = true;
    out.isPdgLike = false;
    out.name = name;
  };

  auto setPdg = [&](const int pdgId, const char* name) {
    out.pdgId = sign * pdgId;
    out.isEposId = false;
    out.isPdgLike = true;
    out.name = name;
  };

  // EPOS basic hadrons.
  if (apdg == 120) setEpos(211, "EPOS PI+");
  else if (apdg == 121) setEpos(213, "EPOS RHO+");
  else if (apdg == 131) setEpos(323, "EPOS K*+");
  else if (apdg == 231) setEpos(313, "EPOS K*0");
  else if (apdg == 221) setEpos(223, "EPOS OMEG");
  else if (apdg == 111 && mass > 0.5) setEpos(113, "EPOS RHO0");

  // EPOS Delta(1232) multiplet and higher Delta-like states.
  else if (apdg == 1111 || apdg == 1112 || apdg == 1113) setEpos(2224, "EPOS DL++");
  else if (apdg == 1121 || apdg == 1124 || apdg == 1126 || apdg == 1128 ||
           apdg == 12214) setEpos(2214, "EPOS DL+");
  else if (apdg == 1221 || apdg == 1224 || apdg == 1226 || apdg == 1228) setEpos(2114, "EPOS DL0");
  else if (apdg == 2221 || apdg == 2222 || apdg == 2223) setEpos(1114, "EPOS DL-");

  // Ambiguous with standard PDG Delta codes; use the stored mass.
  else if (apdg == 1114 && mass > 1.5) setEpos(2224, "EPOS DL++ high");
  else if (apdg == 2224 && mass > 1.5) setEpos(1114, "EPOS DL- high");

  // EPOS/PDG-like N* surrogates used as N pi channels below.
  else if (apdg == 1122 || apdg == 1125 || apdg == 1127 ||
           apdg == 12212 || apdg == 22212 || apdg == 32212) setEpos(12212, "N*+ surrogate");
  else if (apdg == 12112) setPdg(12112, "N(1440)0");
  else if (apdg == 12224) setEpos(12224, "Delta++ surrogate");

  // No EPOS remapping needed. Keep PDG-like codes as-is.
  else setPdg(apdg, "PDG-like");

  return out;
}

bool GetDecayChannelWithChargedPion(const int motherPDG,
                                    const double motherMass,
                                    int& pionPDG,
                                    int& otherPDG,
                                    TRandom3& rng)
{
  pionPDG = 0;
  otherPDG = 0;

  // ----------------------------------------
  // Handle anti-particles by charge conjugation
  // We first decide channel for abs(PDG), then flip signs if mother < 0
  // ----------------------------------------
  const int sign = (motherPDG >= 0) ? +1 : -1;
  const int apdg = std::abs(motherPDG);

  int pion = 0;
  int other = 0;

  auto setDeltaPP = [&]() {
    pion = 211;
    other = 2212;
  };

  auto setDeltaP = [&]() {
    if (rng.Uniform() < 2.0/3.0) {
      pion = 111;
      other = 2212;
    } else {
      pion = 211;
      other = 2112;
    }
  };

  auto setDelta0 = [&]() {
    if (rng.Uniform() < 2.0/3.0) {
      pion = 111;
      other = 2112;
    } else {
      pion = -211;
      other = 2212;
    }
  };

  auto setDeltaM = [&]() {
    pion = -211;
    other = 2112;
  };

  auto setNstarP = [&]() {
    if (rng.Uniform() < 2.0/3.0) {
      pion = 211;
      other = 2112;
    } else {
      pion = 111;
      other = 2212;
    }
  };

  auto setNstar0 = [&]() {
    if (rng.Uniform() < 2.0/3.0) {
      pion = -211;
      other = 2212;
    } else {
      pion = 111;
      other = 2112;
    }
  };

  if (apdg == 211 && motherMass < 0.25) {
    return false;
  }

  // =========================================================
  // Delta(1232) baryons
  // PDG:
  //  2224 = Delta++
  //  2214 = Delta+
  //  2114 = Delta0
  //  1114 = Delta-
  //
  // Standard strong decays:
  // Delta++ -> p pi+
  // Delta+  -> p pi0   or n pi+
  // Delta0  -> n pi0   or p pi-
  // Delta-  -> n pi-
  //
  // Relative weights from isospin CG:
  // Delta+  : p pi0 (2/3), n pi+ (1/3)
  // Delta0  : n pi0 (2/3), p pi- (1/3)
  // =========================================================
  if (apdg == 2224) {
    setDeltaPP();
  }
  else if (apdg == 2214) {
    setDeltaP();
  }
  else if (apdg == 2114) {
    setDelta0();
  }
  else if (apdg == 1114) {
    setDeltaM();
  }
  else if (apdg == 12212 || apdg == 22212 || apdg == 32212) {
    setNstarP();
  }
  else if (apdg == 12112) {
    setNstar0();
  }
  else if (apdg == 12224) {
    setDeltaPP();
  }

  // =========================================================
  // rho(770)
  //  113  = rho0    -> pi+ pi-
  //  213  = rho+    -> pi+ pi0
  // -213  = rho-    -> pi- pi0
  // =========================================================
  else if (apdg == 113) {
    if (rng.Uniform() < 0.5) {
      pion = 211;
      other = -211;
    } else {
      pion = -211;
      other = 211;
    }
  }
  else if (apdg == 213) {
    pion = 211;
    other = 111;
  }

  // omega -> pi+ pi- is a small real branch; here it is used as a two-body
  // charged-pion surrogate because 3-body decays are intentionally disabled.
  else if (apdg == 223) {
    if (rng.Uniform() < 0.5) {
      pion = 211;
      other = -211;
    } else {
      pion = -211;
      other = 211;
    }
  }

  // =========================================================
  // K*(892)
  //  313  = K*0     -> K+ pi-  or K0 pi0
  // -313  = anti-K*0-> K- pi+  or anti-K0 pi0
  //  323  = K*+     -> K+ pi0  or K0 pi+
  // -323  = K*-     -> K- pi0  or anti-K0 pi-
  //
  // Here we pick one standard 2-body pion channel with rough isospin weights:
  // neutral/charged mode ratios ~ 2:1 depending on state.
  // =========================================================
  else if (apdg == 313) {
    double r = rng.Uniform();
    if (r < 2.0/3.0) {
      pion = -211;
      other = 321;
    } else {
      pion = 111;
      other = 311;
    }
  }
  else if (apdg == 323) {
    double r = rng.Uniform();
    if (r < 2.0/3.0) {
      pion = 211;
      other = 311;
    } else {
      pion = 111;
      other = 321;
    }
  }

  // a1(1260) -> rho pi surrogate. The rho is not decayed further here.
  else if (apdg == 20113) {
    if (rng.Uniform() < 0.5) {
      pion = -211;
      other = 213;
    } else {
      pion = 211;
      other = -213;
    }
  }

  // Lambda*/Sigma* two-body pion surrogates.
  else if (apdg == 2122) {
    if (rng.Uniform() < 0.5) {
      pion = -211;
      other = 3222;
    } else {
      pion = 211;
      other = 3112;
    }
  }

  else {
    return false;
  }

  // ----------------------------------------
  // Charge conjugation for anti-particles
  // Note: pi0 and rho0 stay unchanged because PDG = 111,113 are self-conjugate.
  // For K0 / anti-K0, PDG sign flips: 311 <-> -311
  // ----------------------------------------
  if (sign < 0) {
    if (pion != 111 && pion != 113) pion = -pion;
    if (other != 111 && other != 113) other = -other;
  }

  if (std::abs(pion) != 211) {
    return false;
  }

  pionPDG = pion;
  otherPDG = other;
  return true;
}

double GetPDGMass(const int pdg)
{
  const int apdg = std::abs(pdg);

  if (apdg == 211)  return 0.13957039;   // pi+/-
  if (apdg == 111)  return 0.1349768;    // pi0
  if (apdg == 2212) return 0.938272088;  // proton
  if (apdg == 2112) return 0.939565420;  // neutron
  if (apdg == 321)  return 0.493677;     // K+/-
  if (apdg == 311)  return 0.497611;     // K0 / anti-K0
  if (apdg == 113 || apdg == 213) return 0.77526; // rho
  if (apdg == 223)  return 0.78265;       // omega
  if (apdg == 3222) return 1.18937;       // Sigma+
  if (apdg == 3212) return 1.192642;      // Sigma0
  if (apdg == 3112) return 1.197449;      // Sigma-

  return -1.0; // unknown
}
