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
#include <sstream>
#include <vector>

struct PrimPionCentMt {
  int charge = 0;
  TLorentzVector p4;
  TLorentzVector x4;
};

double EtaCentMt(const double px, const double py, const double pz)
{
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= std::abs(pz)) return (pz >= 0.) ? 1e9 : -1e9;
  return 0.5 * std::log((p + pz) / (p - pz));
}

double PairKStarCentMt(const TLorentzVector& p1, const TLorentzVector& p2)
{
  TLorentzVector p1Prf = p1;
  TLorentzVector p2Prf = p2;
  const TLorentzVector pair = p1 + p2;
  p1Prf.Boost(-pair.BoostVector());
  p2Prf.Boost(-pair.BoostVector());
  return 0.5 * (p1Prf.Vect() - p2Prf.Vect()).Mag();
}

double PairMtCentMt(const TLorentzVector& p1, const TLorentzVector& p2)
{
  constexpr double pionMass = 0.13957;
  const TVector3 kT = 0.5 * (p1.Vect() + p2.Vect());
  return std::sqrt(kT.Perp2() + pionMass * pionMass);
}

void FillPairSourceCentMt(const PrimPionCentMt& pi1,
                          const PrimPionCentMt& pi2,
                          TH1D* hRInv,
                          TH1D* hROut,
                          TH1D* hRSide,
                          TH1D* hRLong,
                          TH1D* hRAbs)
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
  hRAbs->Fill(std::sqrt(rOut * rOut + rSide * rSide + rLong * rLong));
}

std::vector<double> GetImpactParameterQuantileEdges(TTree* tree,
                                                    const std::vector<double>& centEdges,
                                                    const Long64_t maxEvents)
{
  Float_t bim = 0.;
  tree->SetBranchAddress("bim", &bim);

  const Long64_t nEntries = tree->GetEntries();
  const Long64_t nToRead = (maxEvents < 0) ? nEntries : std::min(maxEvents, nEntries);

  std::vector<double> bValues;
  bValues.reserve(nToRead);
  for (Long64_t iEvent = 0; iEvent < nToRead; ++iEvent) {
    tree->GetEntry(iEvent);
    bValues.push_back(bim);
  }
  std::sort(bValues.begin(), bValues.end());

  std::vector<double> bEdges;
  bEdges.reserve(centEdges.size());
  for (double cent : centEdges) {
    const double q = cent / 100.0;
    size_t idx = static_cast<size_t>(std::floor(q * (bValues.size() - 1)));
    if (idx >= bValues.size()) idx = bValues.size() - 1;
    bEdges.push_back(bValues[idx]);
  }

  if (!bEdges.empty()) {
    bEdges.front() = std::max(0.0, bEdges.front() - 1e-6);
    bEdges.back() = bEdges.back() + 1e-6;
  }

  tree->ResetBranchAddresses();
  return bEdges;
}

int FindBinIndex(const double value, const std::vector<double>& edges)
{
  if (edges.size() < 2) return -1;
  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    if (value >= edges[i] && value < edges[i + 1]) return static_cast<int>(i);
  }
  return -1;
}

TString BinLabel(const char* prefix, const double low, const double high)
{
  TString label;
  label.Form("%s%.0fto%.0f", prefix, low * 100.0, high * 100.0);
  label.ReplaceAll(".", "p");
  return label;
}

void marco_readEPOSPrimordialPionSource_CentMt(
    TString inputFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/merged_z-EPOS_20250416.root",
    TString outputFileName = "epos_primordial_pion_pair_source_cent_mt_1.root",
    TString treeName = "teposevent",
    Long64_t maxEvents = -1,
    int maxPairsPerChargePerEvent = 100000,
    int finalPionStatus = 0,
    double etaMax = 0.8,
    double pTMin = 0.2,
    double pTMax = 2.0,
    double kStarMin = 0.0,
    double kStarMax = 0.2,
    double maxAbsCoordinate = 1.0e6)
{
  std::vector<double> centEdges = {0, 10, 20, 30, 40, 50};
  std::vector<double> mtEdges = {0.30, 0.40, 0.50, 0.60, 0.70, 0.80};

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

  const std::vector<double> bEdges = GetImpactParameterQuantileEdges(tree, centEdges, maxEvents);
  std::cout << "Impact-parameter quantile edges:" << std::endl;
  for (size_t i = 0; i < centEdges.size(); ++i) {
    std::cout << "  " << centEdges[i] << "% -> b = " << bEdges[i] << " fm" << std::endl;
  }

  constexpr int kMaxParticles = 300000;
  Int_t np = 0;
  Float_t bim = 0.;
  static Float_t px[kMaxParticles], py[kMaxParticles], pz[kMaxParticles], mass[kMaxParticles];
  static Float_t x[kMaxParticles], y[kMaxParticles], z[kMaxParticles], t[kMaxParticles];
  static Int_t id[kMaxParticles], ist[kMaxParticles], ior[kMaxParticles], jor[kMaxParticles];

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
  tree->SetBranchAddress("ior", ior);
  tree->SetBranchAddress("jor", jor);

  const int nCentBins = static_cast<int>(centEdges.size()) - 1;
  const int nMtBins = static_cast<int>(mtEdges.size()) - 1;

  std::vector<std::vector<TH1D*>> hRInv(nCentBins), hROut(nCentBins), hRSide(nCentBins), hRLong(nCentBins), hRAbs(nCentBins);
  std::vector<std::vector<TH1D*>> hKStar(nCentBins), hMt(nCentBins);

  for (int iCent = 0; iCent < nCentBins; ++iCent) {
    hRInv[iCent].resize(nMtBins);
    hROut[iCent].resize(nMtBins);
    hRSide[iCent].resize(nMtBins);
    hRLong[iCent].resize(nMtBins);
    hRAbs[iCent].resize(nMtBins);
    hKStar[iCent].resize(nMtBins);
    hMt[iCent].resize(nMtBins);

    for (int iMt = 0; iMt < nMtBins; ++iMt) {
      TString centLabel;
      centLabel.Form("cent%.0fto%.0f", centEdges[iCent], centEdges[iCent + 1]);
      TString mtLabel;
      mtLabel.Form("mT%.2fto%.2f", mtEdges[iMt], mtEdges[iMt + 1]);
      mtLabel.ReplaceAll(".", "p");
      const TString suffix = centLabel + "_" + mtLabel;

      hRInv[iCent][iMt] = new TH1D("hRInv_PRF_" + suffix, "primordial source in PRF;r*_{inv} (fm);pairs", 800, 0, 400);
      hROut[iCent][iMt] = new TH1D("hROut_LCMS_" + suffix, "primordial source in LCMS;r_{out} (fm);pairs", 800, -200, 200);
      hRSide[iCent][iMt] = new TH1D("hRSide_LCMS_" + suffix, "primordial source in LCMS;r_{side} (fm);pairs", 800, -200, 200);
      hRLong[iCent][iMt] = new TH1D("hRLong_LCMS_" + suffix, "primordial source in LCMS;r_{long} (fm);pairs", 800, -200, 200);
      hRAbs[iCent][iMt] = new TH1D("hRAbs_LCMS_" + suffix, "primordial source in LCMS;|r|_{LCMS} (fm);pairs", 800, 0, 400);
      hKStar[iCent][iMt] = new TH1D("hPairKStar_" + suffix, "selected pairs;k* (GeV/c);pairs", 200, 0, 1.0);
      hMt[iCent][iMt] = new TH1D("hPairMt_" + suffix, "selected pairs;m_{T} (GeV/c^{2});pairs", 200, 0, 2.0);
    }
  }

  TH1D* hImpactParameter = new TH1D("hImpactParameter", "impact parameter;b (fm);events", 120, 0, 24);
  TH2D* hCentMtPairs = new TH2D("hCentMtPairs", "accepted pairs;centrality bin;m_{T} bin", nCentBins, 0, nCentBins, nMtBins, 0, nMtBins);

  TRandom3 rng(0);
  const Long64_t nEntries = tree->GetEntries();
  const Long64_t nEventsToRead = (maxEvents < 0) ? nEntries : std::min(maxEvents, nEntries);

  Long64_t acceptedPairs = 0;
  Long64_t testedPairs = 0;
  Long64_t rejectedDecayPions = 0;
  Long64_t rejectedBadCoordinates = 0;

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
    const int iCent = FindBinIndex(bim, bEdges);
    if (iCent < 0) continue;
    hImpactParameter->Fill(bim);

    std::vector<PrimPionCentMt> plus;
    std::vector<PrimPionCentMt> minus;

    for (int i = 0; i < np; ++i) {
      if (std::abs(id[i]) != 120) continue;
      if (finalPionStatus >= 0 && ist[i] != finalPionStatus) continue;

      if (ior[i] >= 0 && jor[i] < 0) {
        ++rejectedDecayPions;
        continue;
      }
      if (!(ior[i] < 0 && jor[i] < 0)) continue;
      if (!goodCoordinates(i)) {
        ++rejectedBadCoordinates;
        continue;
      }

      const double pT = std::sqrt(px[i] * px[i] + py[i] * py[i]);
      if (pT < pTMin || pT >= pTMax) continue;
      const double eta = EtaCentMt(px[i], py[i], pz[i]);
      if (std::abs(eta) > etaMax) continue;

      const double energy =
          std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] + mass[i] * mass[i]);

      PrimPionCentMt pion;
      pion.charge = (id[i] > 0) ? +1 : -1;
      pion.p4.SetPxPyPzE(px[i], py[i], pz[i], energy);
      pion.x4.SetXYZT(x[i], y[i], z[i], t[i]);
      if (pion.charge > 0) plus.push_back(pion);
      else minus.push_back(pion);
    }

    auto processSameCharge = [&](const std::vector<PrimPionCentMt>& pions) {
      const Long64_t n = static_cast<Long64_t>(pions.size());
      if (n < 2) return;
      const Long64_t nPossiblePairs = n * (n - 1) / 2;
      const Long64_t nTrials = maxPairsPerChargePerEvent > 0
                                   ? std::min<Long64_t>(maxPairsPerChargePerEvent, nPossiblePairs)
                                   : nPossiblePairs;

      for (Long64_t trial = 0; trial < nTrials; ++trial) {
        Long64_t i = 0;
        Long64_t j = 0;
        if (maxPairsPerChargePerEvent > 0 && nPossiblePairs > maxPairsPerChargePerEvent) {
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
        const double kStar = PairKStarCentMt(pions[i].p4, pions[j].p4);
        if (kStar < kStarMin || kStar >= kStarMax) continue;

        const double mT = PairMtCentMt(pions[i].p4, pions[j].p4);
        const int iMt = FindBinIndex(mT, mtEdges);
        if (iMt < 0) continue;

        hKStar[iCent][iMt]->Fill(kStar);
        hMt[iCent][iMt]->Fill(mT);
        FillPairSourceCentMt(pions[i], pions[j],
                             hRInv[iCent][iMt],
                             hROut[iCent][iMt],
                             hRSide[iCent][iMt],
                             hRLong[iCent][iMt],
                             hRAbs[iCent][iMt]);
        hCentMtPairs->Fill(iCent, iMt);
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
  hCentMtPairs->Write();
  for (int iCent = 0; iCent < nCentBins; ++iCent) {
    for (int iMt = 0; iMt < nMtBins; ++iMt) {
      hRInv[iCent][iMt]->Write();
      hROut[iCent][iMt]->Write();
      hRSide[iCent][iMt]->Write();
      hRLong[iCent][iMt]->Write();
      hRAbs[iCent][iMt]->Write();
      hKStar[iCent][iMt]->Write();
      hMt[iCent][iMt]->Write();
    }
  }
  outputFile->Close();

  std::cout << "Saved " << outputFileName << std::endl;
  std::cout << "Tested same-charge pairs: " << testedPairs << std::endl;
  std::cout << "Accepted pairs: " << acceptedPairs << std::endl;
  std::cout << "Rejected single-father decay pions: " << rejectedDecayPions << std::endl;
  std::cout << "Rejected bad-coordinate pions: " << rejectedBadCoordinates << std::endl;

  inputFile->Close();
}
