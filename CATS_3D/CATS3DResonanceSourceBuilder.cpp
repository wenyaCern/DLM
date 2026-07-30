#include "FitCATS3D.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

bool HasBranch(TTree* tree, const char* name)
{
    return tree && tree->GetBranch(name);
}

void RequireBranch(TTree* tree, const char* name)
{
    if(!HasBranch(tree, name)){
        throw std::runtime_error(Form("LoadCATS3DResonanceSource: missing branch '%s'", name));
    }
}

unsigned LoadResonancePR(DLM_CleverReso3D& source,
                         const FitCATS3DConfig& cfg,
                         const TString& fileName,
                         const TString& treeName)
{
    TFile fin(fileName, "READ");
    if(fin.IsZombie()){
        throw std::runtime_error(Form("LoadResonancePR: cannot open '%s'",
                                      fileName.Data()));
    }

    TTree* tree = dynamic_cast<TTree*>(fin.Get(treeName));
    if(!tree){
        throw std::runtime_error(Form("LoadResonancePR: missing TTree '%s'",
                                      treeName.Data()));
    }

    const char* needed[] = {
        "mother_pdg2", "motherInDaughterLCMS_pOut2",
        "motherInDaughterLCMS_pSide2", "motherInDaughterLCMS_pLong2",
        "m2", "oS_dot_oD", "sS_dot_oD", "oS_dot_sD",
        "sS_dot_sD", "longScale"
    };
    for(const char* branch : needed) RequireBranch(tree, branch);

    std::vector<int>* motherPdg2 = nullptr;
    std::vector<double>* pOut2 = nullptr;
    std::vector<double>* pSide2 = nullptr;
    std::vector<double>* pLong2 = nullptr;
    std::vector<double>* m2 = nullptr;
    std::vector<double>* oS_dot_oD = nullptr;
    std::vector<double>* sS_dot_oD = nullptr;
    std::vector<double>* oS_dot_sD = nullptr;
    std::vector<double>* sS_dot_sD = nullptr;
    std::vector<double>* longScale = nullptr;

    tree->SetBranchAddress("mother_pdg2", &motherPdg2);
    tree->SetBranchAddress("motherInDaughterLCMS_pOut2", &pOut2);
    tree->SetBranchAddress("motherInDaughterLCMS_pSide2", &pSide2);
    tree->SetBranchAddress("motherInDaughterLCMS_pLong2", &pLong2);
    tree->SetBranchAddress("m2", &m2);
    tree->SetBranchAddress("oS_dot_oD", &oS_dot_oD);
    tree->SetBranchAddress("sS_dot_oD", &sS_dot_oD);
    tree->SetBranchAddress("oS_dot_sD", &oS_dot_sD);
    tree->SetBranchAddress("sS_dot_sD", &sS_dot_sD);
    tree->SetBranchAddress("longScale", &longScale);

    unsigned loaded = 0;
    const Long64_t nEntries = tree->GetEntries();

    for(Long64_t iEntry=0; iEntry<nEntries; ++iEntry){
        tree->GetEntry(iEntry);
        if(!motherPdg2 || !pOut2 || !pSide2 || !pLong2 || !m2 ||
           !oS_dot_oD || !sS_dot_oD || !oS_dot_sD || !sS_dot_sD || !longScale) continue;

        const size_t nPairs = std::min({
            motherPdg2->size(), pOut2->size(), pSide2->size(), pLong2->size(),
            m2->size(), oS_dot_oD->size(), sS_dot_oD->size(),
            oS_dot_sD->size(), sS_dot_sD->size(), longScale->size()
        });

        for(size_t iPair=0; iPair<nPairs; ++iPair){
            source.AddBGT_PR(TVector3(pOut2->at(iPair),
                                      pSide2->at(iPair),
                                      pLong2->at(iPair)),
                             m2->at(iPair),
                             cfg.ResonanceDefaultTauFm,
                             oS_dot_oD->at(iPair),
                             sS_dot_oD->at(iPair),
                             oS_dot_sD->at(iPair),
                             sS_dot_sD->at(iPair),
                             longScale->at(iPair));
            ++loaded;
            if(cfg.ResonanceMaxPairs && loaded>=cfg.ResonanceMaxPairs) return loaded;
        }
    }

    return loaded;
}

unsigned LoadResonanceRR(DLM_CleverReso3D& source,
                         const FitCATS3DConfig& cfg,
                         const TString& fileName,
                         const TString& treeName)
{
    TFile fin(fileName, "READ");
    if(fin.IsZombie()){
        throw std::runtime_error(Form("LoadResonanceRR: cannot open '%s'",
                                      fileName.Data()));
    }

    TTree* tree = dynamic_cast<TTree*>(fin.Get(treeName));
    if(!tree){
        throw std::runtime_error(Form("LoadResonanceRR: missing TTree '%s'",
                                      treeName.Data()));
    }

    const char* needed[] = {
        "mother_pdg1", "motherInDaughterLCMS_pOut1",
        "motherInDaughterLCMS_pSide1", "motherInDaughterLCMS_pLong1",
        "m1", "mother_pdg2", "motherInDaughterLCMS_pOut2",
        "motherInDaughterLCMS_pSide2", "motherInDaughterLCMS_pLong2",
        "m2", "oS_dot_oD", "sS_dot_oD", "oS_dot_sD",
        "sS_dot_sD", "longScale"
    };
    for(const char* branch : needed) RequireBranch(tree, branch);

    std::vector<int>* motherPdg1 = nullptr;
    std::vector<double>* pOut1 = nullptr;
    std::vector<double>* pSide1 = nullptr;
    std::vector<double>* pLong1 = nullptr;
    std::vector<double>* m1 = nullptr;
    std::vector<int>* motherPdg2 = nullptr;
    std::vector<double>* pOut2 = nullptr;
    std::vector<double>* pSide2 = nullptr;
    std::vector<double>* pLong2 = nullptr;
    std::vector<double>* m2 = nullptr;
    std::vector<double>* oS_dot_oD = nullptr;
    std::vector<double>* sS_dot_oD = nullptr;
    std::vector<double>* oS_dot_sD = nullptr;
    std::vector<double>* sS_dot_sD = nullptr;
    std::vector<double>* longScale = nullptr;

    tree->SetBranchAddress("mother_pdg1", &motherPdg1);
    tree->SetBranchAddress("motherInDaughterLCMS_pOut1", &pOut1);
    tree->SetBranchAddress("motherInDaughterLCMS_pSide1", &pSide1);
    tree->SetBranchAddress("motherInDaughterLCMS_pLong1", &pLong1);
    tree->SetBranchAddress("m1", &m1);
    tree->SetBranchAddress("mother_pdg2", &motherPdg2);
    tree->SetBranchAddress("motherInDaughterLCMS_pOut2", &pOut2);
    tree->SetBranchAddress("motherInDaughterLCMS_pSide2", &pSide2);
    tree->SetBranchAddress("motherInDaughterLCMS_pLong2", &pLong2);
    tree->SetBranchAddress("m2", &m2);
    tree->SetBranchAddress("oS_dot_oD", &oS_dot_oD);
    tree->SetBranchAddress("sS_dot_oD", &sS_dot_oD);
    tree->SetBranchAddress("oS_dot_sD", &oS_dot_sD);
    tree->SetBranchAddress("sS_dot_sD", &sS_dot_sD);
    tree->SetBranchAddress("longScale", &longScale);

    unsigned loaded = 0;
    const Long64_t nEntries = tree->GetEntries();

    for(Long64_t iEntry=0; iEntry<nEntries; ++iEntry){
        tree->GetEntry(iEntry);
        if(!motherPdg1 || !pOut1 || !pSide1 || !pLong1 || !m1 ||
           !motherPdg2 || !pOut2 || !pSide2 || !pLong2 || !m2 ||
           !oS_dot_oD || !sS_dot_oD || !oS_dot_sD || !sS_dot_sD || !longScale) continue;

        const size_t nPairs = std::min({
            motherPdg1->size(), pOut1->size(), pSide1->size(), pLong1->size(),
            m1->size(), motherPdg2->size(), pOut2->size(), pSide2->size(),
            pLong2->size(), m2->size(), oS_dot_oD->size(), sS_dot_oD->size(),
            oS_dot_sD->size(), sS_dot_sD->size(), longScale->size()
        });

        for(size_t iPair=0; iPair<nPairs; ++iPair){
            source.AddBGT_RR(TVector3(pOut1->at(iPair),
                                      pSide1->at(iPair),
                                      pLong1->at(iPair)),
                             m1->at(iPair),
                             cfg.ResonanceDefaultTauFm,
                             TVector3(pOut2->at(iPair),
                                      pSide2->at(iPair),
                                      pLong2->at(iPair)),
                             m2->at(iPair),
                             cfg.ResonanceDefaultTauFm,
                             oS_dot_oD->at(iPair),
                             sS_dot_oD->at(iPair),
                             oS_dot_sD->at(iPair),
                             sS_dot_sD->at(iPair),
                             longScale->at(iPair));
            ++loaded;
            if(cfg.ResonanceMaxPairs && loaded>=cfg.ResonanceMaxPairs) return loaded;
        }
    }

    return loaded;
}

} // namespace

std::unique_ptr<DLM_CleverReso3D> BuildCATS3DResonanceSource(const FitCATS3DConfig& cfg)
{
    if(!cfg.UseResonanceSource) return nullptr;

    std::unique_ptr<DLM_CleverReso3D> source(new DLM_CleverReso3D());
    source->SetUpReso(0, cfg.ResonanceWeight1);
    source->SetUpReso(1, cfg.ResonanceWeight2);
    source->InitNumMcIter(cfg.ResonanceMcIter);
    source->InitScale(cfg.ResonanceScaleBins,
                      cfg.ResonanceScaleMinFm,
                      cfg.ResonanceScaleMaxFm);
    source->InitRad(cfg.ResonanceRadBins,
                    cfg.ResonanceRadMinFm,
                    cfg.ResonanceRadMaxFm);
    source->InitStability(cfg.ResonanceAlphaBins,
                          cfg.ResonanceAlphaMin,
                          cfg.ResonanceAlphaMax);
    source->SetCoreMode(cfg.ResonanceUseEllipsoidalCore ?
                        DLM_CleverReso3D::cmEllipsoidalLevy :
                        DLM_CleverReso3D::cmSeparableStable);
    source->SetUseLCMS(true);
    source->SetUseModelCore(false);
    source->SetDebugMode(cfg.ResonanceDebug);

    unsigned nPR = 0;
    unsigned nRR = 0;
    const bool needPR = ((1.0-cfg.ResonanceWeight1)*cfg.ResonanceWeight2 > 0.0) ||
                        (cfg.ResonanceWeight1*(1.0-cfg.ResonanceWeight2) > 0.0);
    if(cfg.UseResonancePR && needPR){
        nPR = LoadResonancePR(*source, cfg,
                              cfg.ResonancePRFileName,
                              cfg.ResonancePRTreeName);
        if(!nPR){
            throw std::runtime_error(Form("BuildCATS3DResonanceSource: requested PR/RP resonance source but loaded zero pairs from '%s' :: '%s'",
                                          cfg.ResonancePRFileName.Data(),
                                          cfg.ResonancePRTreeName.Data()));
        }
    }
    if(cfg.UseResonanceRR && cfg.ResonanceWeight1 > 0.0 && cfg.ResonanceWeight2 > 0.0){
        nRR = LoadResonanceRR(*source, cfg,
                              cfg.ResonanceRRFileName,
                              cfg.ResonanceRRTreeName);
        if(!nRR){
            throw std::runtime_error(Form("BuildCATS3DResonanceSource: requested RR resonance source but loaded zero pairs from '%s' :: '%s'",
                                          cfg.ResonanceRRFileName.Data(),
                                          cfg.ResonanceRRTreeName.Data()));
        }
    }

    std::cout << "Loaded CATS3D resonance source kernels: PR=" << nPR
              << " RR=" << nRR
              << " weights=(" << cfg.ResonanceWeight1
              << ", " << cfg.ResonanceWeight2 << ")"
              << " alphaRange=[" << cfg.ResonanceAlphaMin
              << ", " << cfg.ResonanceAlphaMax << "]"
              << " coreMode="
              << (cfg.ResonanceUseEllipsoidalCore ? "ellipsoidal" : "separable")
              << "\n";

    return source;
}
