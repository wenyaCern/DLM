#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

struct EposPion {
  int charge = 0;
  TLorentzVector p4;
  TLorentzVector x4;
};

double Pseudorapidity(const double px, const double py, const double pz)
{
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= std::abs(pz)) return (pz >= 0.) ? 1e9 : -1e9;
  return 0.5 * std::log((p + pz) / (p - pz));
}

double PairKStar(const TLorentzVector& p1, const TLorentzVector& p2)
{
  TLorentzVector p1Prf = p1;
  TLorentzVector p2Prf = p2;
  const TLorentzVector pair = p1 + p2;
  p1Prf.Boost(-pair.BoostVector());
  p2Prf.Boost(-pair.BoostVector());
  return 0.5 * (p1Prf.Vect() - p2Prf.Vect()).Mag();
}

double PairMt(const TLorentzVector& p1, const TLorentzVector& p2)
{
  const TVector3 kT = 0.5 * (p1.Vect() + p2.Vect());
  const double pionMass = 0.13957;
  return std::sqrt(kT.Perp2() + pionMass * pionMass);
}

void FillPairSource(const EposPion& pi1,
                    const EposPion& pi2,
                    TH1D* hRInv,
                    TH1D* hRLong,
                    TH1D* hROut,
                    TH1D* hRSide,
                    TH1D* hRAbsLCMS,
                    TH2D* hROutRSide,
                    TH2D* hRLongROut)
{
  const TLorentzVector pairP = pi1.p4 + pi2.p4;
  const TLorentzVector relX = pi1.x4 - pi2.x4;

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

void marco_readEPOSPionSource(
    TString inputFileName = "./merged_z-EPOS_20250416.root",
    TString outputFileName = "epos_pion_source.root",
    TString treeName = "teposevent",
    Long64_t maxEvents = -1,
    int maxPairsPerEvent = 200000,
    double bMin = 0.0,
    double bMax = 3.5,
    double etaMax = 0.8,
    double pTMin = 0.2,
    double pTMax = 2.0,
    double kStarMin = 0.0,
    double kStarMax = 0.2,
    double mTMin = 0.0,
    double mTMax = 10.0,
    bool sameChargeOnly = true)
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
  static Int_t id[kMaxParticles], ist[kMaxParticles];

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

  TH1D* hNumPions = new TH1D("hNumPions", "selected pion multiplicity;N_{#pi};events", 400, 0, 4000);
  TH1D* hB = new TH1D("hImpactParameter", "impact parameter;b (fm);events", 120, 0, 24);
  TH1D* hKStar = new TH1D("hPairKStar", "selected same-event pion pairs;k* (GeV/c);pairs", 200, 0, 1.0);
  TH1D* hMt = new TH1D("hPairMt", "selected same-event pion pairs;m_{T} (GeV/c^{2});pairs", 200, 0, 2.0);

  TH1D* hRInv = new TH1D("hRInv_PRF", "pair source in PRF;r*_{inv} (fm);pairs", 800, 0, 400);
  TH1D* hRLong = new TH1D("hRLong_LCMS", "pair source in LCMS;r_{long} (fm);pairs", 800, -200, 200);
  TH1D* hROut = new TH1D("hROut_LCMS", "pair source in LCMS;r_{out} (fm);pairs", 800, -200, 200);
  TH1D* hRSide = new TH1D("hRSide_LCMS", "pair source in LCMS;r_{side} (fm);pairs", 800, -200, 200);
  TH1D* hRAbsLCMS = new TH1D("hRAbs_LCMS", "pair source in LCMS;|r|_{LCMS} (fm);pairs", 800, 0, 400);
  TH2D* hROutRSide = new TH2D("hROutRSide_LCMS", "pair source in LCMS;r_{out} (fm);r_{side} (fm)", 300, -150, 150, 300, -150, 150);
  TH2D* hRLongROut = new TH2D("hRLongROut_LCMS", "pair source in LCMS;r_{long} (fm);r_{out} (fm)", 300, -150, 150, 300, -150, 150);

  TRandom3 rng(0);
  const Long64_t nEntries = tree->GetEntries();
  const Long64_t nEventsToRead = (maxEvents < 0) ? nEntries : std::min(maxEvents, nEntries);

  Long64_t totalPairsTested = 0;
  Long64_t totalPairsAccepted = 0;
  Long64_t totalEventsAccepted = 0;

  for (Long64_t iEvent = 0; iEvent < nEventsToRead; ++iEvent) {
    tree->GetEntry(iEvent);

    if (bim < bMin || bim >= bMax) continue;
    if (np > kMaxParticles) {
      std::cerr << "WARNING: event " << iEvent << " has np = " << np
                << ", larger than buffer. Skipping." << std::endl;
      continue;
    }

    std::vector<EposPion> piPlus;
    std::vector<EposPion> piMinus;
    piPlus.reserve(2000);
    piMinus.reserve(2000);

    for (int i = 0; i < np; ++i) {
      if (std::abs(id[i]) != 120) continue;
      if (ist[i] != 7) continue;

      const double pT = std::sqrt(px[i] * px[i] + py[i] * py[i]);
      if (pT < pTMin || pT >= pTMax) continue;

      const double eta = Pseudorapidity(px[i], py[i], pz[i]);
      if (std::abs(eta) > etaMax) continue;

      const double energy = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] + mass[i] * mass[i]);

      EposPion pion;
      pion.charge = (id[i] > 0) ? +1 : -1;
      pion.p4.SetPxPyPzE(px[i], py[i], pz[i], energy);
      pion.x4.SetXYZT(x[i], y[i], z[i], t[i]);

      if (pion.charge > 0) piPlus.push_back(pion);
      else piMinus.push_back(pion);
    }

    hB->Fill(bim);
    hNumPions->Fill(piPlus.size() + piMinus.size());
    ++totalEventsAccepted;

    auto processList = [&](const std::vector<EposPion>& pions) {
      const Long64_t nPions = static_cast<Long64_t>(pions.size());
      if (nPions < 2) return;
      const Long64_t nPossiblePairs = nPions * (nPions - 1) / 2;
      const Long64_t nTrials = (maxPairsPerEvent > 0)
                                   ? std::min<Long64_t>(maxPairsPerEvent, nPossiblePairs)
                                   : nPossiblePairs;

      for (Long64_t iTrial = 0; iTrial < nTrials; ++iTrial) {
        Long64_t i = 0;
        Long64_t j = 0;
        if (maxPairsPerEvent > 0 && nPossiblePairs > maxPairsPerEvent) {
          i = rng.Integer(nPions);
          j = rng.Integer(nPions - 1);
          if (j >= i) ++j;
        } else {
          Long64_t count = 0;
          bool found = false;
          for (Long64_t a = 0; a < nPions && !found; ++a) {
            for (Long64_t b = a + 1; b < nPions; ++b) {
              if (count == iTrial) {
                i = a;
                j = b;
                found = true;
                break;
              }
              ++count;
            }
          }
        }

        ++totalPairsTested;
        const double kStar = PairKStar(pions[i].p4, pions[j].p4);
        if (kStar < kStarMin || kStar >= kStarMax) continue;

        const double mT = PairMt(pions[i].p4, pions[j].p4);
        if (mT < mTMin || mT >= mTMax) continue;

        hKStar->Fill(kStar);
        hMt->Fill(mT);
        FillPairSource(pions[i], pions[j], hRInv, hRLong, hROut, hRSide,
                       hRAbsLCMS, hROutRSide, hRLongROut);
        ++totalPairsAccepted;
      }
    };

    processList(piPlus);
    if (sameChargeOnly) processList(piMinus);

    if ((iEvent + 1) % 50 == 0) {
      std::cout << "Processed " << iEvent + 1 << " / " << nEventsToRead
                << " events, accepted pairs = " << totalPairsAccepted << std::endl;
    }
  }

  TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
  hNumPions->Write();
  hB->Write();
  hKStar->Write();
  hMt->Write();
  hRInv->Write();
  hRLong->Write();
  hROut->Write();
  hRSide->Write();
  hRAbsLCMS->Write();
  hROutRSide->Write();
  hRLongROut->Write();
  outputFile->Close();

  std::cout << "Saved " << outputFileName << std::endl;
  std::cout << "Accepted events: " << totalEventsAccepted << std::endl;
  std::cout << "Tested pairs: " << totalPairsTested << std::endl;
  std::cout << "Accepted pairs after k*/mT cuts: " << totalPairsAccepted << std::endl;

  inputFile->Close();
}
