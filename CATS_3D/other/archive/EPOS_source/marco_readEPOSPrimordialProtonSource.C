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
#include <vector>

struct PrimordialProton {
  int charge = 0;
  int originType = 0;
  TLorentzVector p4;
  TLorentzVector x4;
};

double EposProtonEta(const double px, const double py, const double pz)
{
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= std::abs(pz)) return (pz >= 0.) ? 1e9 : -1e9;
  return 0.5 * std::log((p + pz) / (p - pz));
}

double EposProtonPairKStar(const TLorentzVector& p1, const TLorentzVector& p2)
{
  TLorentzVector p1Prf = p1;
  TLorentzVector p2Prf = p2;
  const TLorentzVector pair = p1 + p2;
  p1Prf.Boost(-pair.BoostVector());
  p2Prf.Boost(-pair.BoostVector());
  return 0.5 * (p1Prf.Vect() - p2Prf.Vect()).Mag();
}

double EposProtonPairMt(const TLorentzVector& p1, const TLorentzVector& p2)
{
  const TVector3 kT = 0.5 * (p1.Vect() + p2.Vect());
  constexpr double protonMass = 0.938272088;
  return std::sqrt(kT.Perp2() + protonMass * protonMass);
}

void FillPrimordialProtonPairSource(const PrimordialProton& pr1,
                                    const PrimordialProton& pr2,
                                    TH1D* hRInv,
                                    TH1D* hROut,
                                    TH1D* hRSide,
                                    TH1D* hRLong,
                                    TH1D* hRAbsLCMS,
                                    TH2D* hROutRSide,
                                    TH2D* hRLongROut)
{
  const TLorentzVector pairP = pr1.p4 + pr2.p4;
  const TLorentzVector relX = pr1.x4 - pr2.x4;

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

void BuildProtonRadialSourceOutputs(const TH1D* hRInv,
                                    TH1D*& hFourPiR2S,
                                    TH1D*& hSourceDensity)
{
  hFourPiR2S = dynamic_cast<TH1D*>(hRInv->Clone("hFourPiR2S_PRF"));
  hFourPiR2S->SetTitle("normalized PRF radial source; r*_{inv} (fm);4#pi r^{2} S(r*)");

  hSourceDensity = dynamic_cast<TH1D*>(hRInv->Clone("hSourceDensityS_PRF"));
  hSourceDensity->SetTitle("PRF source density; r*_{inv} (fm);S(r*) (fm^{-3})");

  const double norm = hRInv->Integral();
  for (int iBin = 1; iBin <= hRInv->GetNbinsX(); ++iBin) {
    const double counts = hRInv->GetBinContent(iBin);
    const double countsErr = hRInv->GetBinError(iBin);
    const double r = hRInv->GetBinCenter(iBin);
    const double dr = hRInv->GetBinWidth(iBin);

    double shellDensity = 0.;
    double shellDensityErr = 0.;
    double sourceDensity = 0.;
    double sourceDensityErr = 0.;

    if (norm > 0. && dr > 0.) {
      shellDensity = counts / (norm * dr);
      shellDensityErr = countsErr / (norm * dr);

      constexpr double pi = 3.14159265358979323846;
      const double shellVolumeFactor = 4. * pi * r * r;
      if (shellVolumeFactor > 0.) {
        sourceDensity = shellDensity / shellVolumeFactor;
        sourceDensityErr = shellDensityErr / shellVolumeFactor;
      }
    }

    hFourPiR2S->SetBinContent(iBin, shellDensity);
    hFourPiR2S->SetBinError(iBin, shellDensityErr);
    hSourceDensity->SetBinContent(iBin, sourceDensity);
    hSourceDensity->SetBinError(iBin, sourceDensityErr);
  }
}

void marco_readEPOSPrimordialProtonSource(
    TString inputFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/merged_z-EPOS_20250416.root",
    TString outputFileName = "epos_primordial_proton_pair_source_kStarCutoff20.root",
    TString treeName = "teposevent",
    Long64_t maxEvents = -1,
    int maxPairsPerChargePerEvent = 200000,
    int finalProtonStatus = 0,
    double bMin = 0.0,
    double bMax = 3.5,
    double etaMax = 0.8,
    double pTMin = 0.2,
    double pTMax = 5.0,
    double kStarMin = 0.0,
    double kStarMax = 0.02,
    double mTMin = 0.0,
    double mTMax = 10.0,
    double maxAbsCoordinate = 1.0e6)
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
  TH1D* hSelectedProtons = new TH1D("hSelectedPrimordialProtons", "selected primordial proton multiplicity;N_{p};events", 400, 0, 4000);
  TH1D* hOriginType = new TH1D("hPrimordialOriginType", "primordial proton origin type;ity;counts", 100, 0, 100);
  TH1D* hPairKStar = new TH1D("hPairKStar", "primordial proton pairs;k* (GeV/c);pairs", 200, 0, 1.0);
  TH1D* hPairMt = new TH1D("hPairMt", "primordial proton pairs;m_{T} (GeV/c^{2});pairs", 200, 0, 3.0);

  TH1D* hRInv = new TH1D("hRInv_PRF", "primordial proton pair source in PRF;r*_{inv} (fm);pairs", 800, 0, 400);
  TH1D* hROut = new TH1D("hROut_LCMS", "primordial proton pair source in LCMS;r_{out} (fm);pairs", 800, -200, 200);
  TH1D* hRSide = new TH1D("hRSide_LCMS", "primordial proton pair source in LCMS;r_{side} (fm);pairs", 800, -200, 200);
  TH1D* hRLong = new TH1D("hRLong_LCMS", "primordial proton pair source in LCMS;r_{long} (fm);pairs", 800, -200, 200);
  TH1D* hRAbsLCMS = new TH1D("hRAbs_LCMS", "primordial proton pair source in LCMS;|r|_{LCMS} (fm);pairs", 800, 0, 400);
  TH2D* hROutRSide = new TH2D("hROutRSide_LCMS", "primordial proton pair source in LCMS;r_{out} (fm);r_{side} (fm)", 300, -150, 150, 300, -150, 150);
  TH2D* hRLongROut = new TH2D("hRLongROut_LCMS", "primordial proton pair source in LCMS;r_{long} (fm);r_{out} (fm)", 300, -150, 150, 300, -150, 150);

  TRandom3 rng(0);
  const Long64_t nEntries = tree->GetEntries();
  const Long64_t nEventsToRead = (maxEvents < 0) ? nEntries : std::min(maxEvents, nEntries);

  Long64_t totalEventsAccepted = 0;
  Long64_t totalPairsTested = 0;
  Long64_t totalPairsAccepted = 0;
  Long64_t rejectedDecayProtons = 0;
  Long64_t rejectedBadCoordinates = 0;

  auto hasGoodCoordinates = [&](const int i) {
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

    std::vector<PrimordialProton> protons;
    std::vector<PrimordialProton> antiprotons;

    for (int i = 0; i < np; ++i) {
      if (std::abs(id[i]) != 1120) continue;
      if (finalProtonStatus >= 0 && ist[i] != finalProtonStatus) continue;

      const bool isSingleFatherDecay = (ior[i] >= 0 && jor[i] < 0);
      if (isSingleFatherDecay) {
        ++rejectedDecayProtons;
        continue;
      }

      if (!hasGoodCoordinates(i)) {
        ++rejectedBadCoordinates;
        continue;
      }

      const double pT = std::sqrt(px[i] * px[i] + py[i] * py[i]);
      if (pT < pTMin || pT >= pTMax) continue;

      const double eta = EposProtonEta(px[i], py[i], pz[i]);
      if (std::abs(eta) > etaMax) continue;

      const double energy =
          std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] + mass[i] * mass[i]);

      PrimordialProton proton;
      proton.charge = (id[i] > 0) ? +1 : -1;
      proton.originType = ity[i];
      proton.p4.SetPxPyPzE(px[i], py[i], pz[i], energy);
      proton.x4.SetXYZT(x[i], y[i], z[i], t[i]);

      hOriginType->Fill(ity[i]);
      if (proton.charge > 0) protons.push_back(proton);
      else antiprotons.push_back(proton);
    }

    hImpactParameter->Fill(bim);
    hSelectedProtons->Fill(protons.size() + antiprotons.size());
    ++totalEventsAccepted;

    auto processSameChargePairs = [&](const std::vector<PrimordialProton>& particles) {
      const Long64_t nParticles = static_cast<Long64_t>(particles.size());
      if (nParticles < 2) return;

      const Long64_t nPossiblePairs = nParticles * (nParticles - 1) / 2;
      const Long64_t nTrials = maxPairsPerChargePerEvent > 0
                                   ? std::min<Long64_t>(maxPairsPerChargePerEvent, nPossiblePairs)
                                   : nPossiblePairs;

      for (Long64_t iTrial = 0; iTrial < nTrials; ++iTrial) {
        Long64_t i = 0;
        Long64_t j = 0;

        if (maxPairsPerChargePerEvent > 0 && nPossiblePairs > maxPairsPerChargePerEvent) {
          i = rng.Integer(nParticles);
          j = rng.Integer(nParticles - 1);
          if (j >= i) ++j;
        } else {
          Long64_t count = 0;
          bool found = false;
          for (Long64_t a = 0; a < nParticles && !found; ++a) {
            for (Long64_t b = a + 1; b < nParticles; ++b) {
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

        const double kStar = EposProtonPairKStar(particles[i].p4, particles[j].p4);
        if (kStar < kStarMin || kStar >= kStarMax) continue;

        const double mT = EposProtonPairMt(particles[i].p4, particles[j].p4);
        if (mT < mTMin || mT >= mTMax) continue;

        hPairKStar->Fill(kStar);
        hPairMt->Fill(mT);
        FillPrimordialProtonPairSource(particles[i], particles[j],
                                       hRInv, hROut, hRSide, hRLong,
                                       hRAbsLCMS, hROutRSide, hRLongROut);
        ++totalPairsAccepted;
      }
    };

    processSameChargePairs(protons);
    processSameChargePairs(antiprotons);

    if ((iEvent + 1) % 50 == 0) {
      std::cout << "Processed " << iEvent + 1 << " / " << nEventsToRead
                << " events, accepted pairs = " << totalPairsAccepted << std::endl;
    }
  }

  TH1D* hFourPiR2S = nullptr;
  TH1D* hSourceDensity = nullptr;
  BuildProtonRadialSourceOutputs(hRInv, hFourPiR2S, hSourceDensity);

  TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
  hImpactParameter->Write();
  hSelectedProtons->Write();
  hOriginType->Write();
  hPairKStar->Write();
  hPairMt->Write();
  hRInv->Write();
  hFourPiR2S->Write();
  hSourceDensity->Write();
  hROut->Write();
  hRSide->Write();
  hRLong->Write();
  hRAbsLCMS->Write();
  hROutRSide->Write();
  hRLongROut->Write();
  outputFile->Close();

  std::cout << "Saved " << outputFileName << std::endl;
  std::cout << "Accepted events: " << totalEventsAccepted << std::endl;
  std::cout << "Rejected single-father decay protons: " << rejectedDecayProtons << std::endl;
  std::cout << "Rejected bad-coordinate protons: " << rejectedBadCoordinates << std::endl;
  std::cout << "Tested same-charge pairs: " << totalPairsTested << std::endl;
  std::cout << "Accepted pairs after k*/mT cuts: " << totalPairsAccepted << std::endl;

  inputFile->Close();
}
