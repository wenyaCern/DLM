#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TString.h"
#include "TSystem.h"
#include <memory>
#include <stdexcept>

void prepare_ft0a_merged_mt(
 const TString input="File/out_3DCF_mT_cent_qn_phi_tofOnly_LHC23_FT0CqVec_FT0APlane.root",
 const TString lambdaInput="File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root",
 int cent=1,int mtFirst=2,int mtLast=4,
 const TString dataOutput="File/out_3DCF_FT0APlane_cent1_mt0p4to0p7.root",
 const TString lambdaOutput="File/lam/lambda_FT0APlane_cent1_mt0p4to0p7.root"){
 TFile fi(input,"READ");if(fi.IsZombie())throw std::runtime_error("cannot open input");
 TFile fo(dataOutput,"RECREATE");
 for(int p=0;p<6;++p){
  std::unique_ptr<TH3> mix;
  for(int m=mtFirst;m<=mtLast;++m){auto*h=dynamic_cast<TH3*>(fi.Get(Form("H3_mixed_cent%d_mt%d_phi%d",cent,m,p)));if(!h)throw std::runtime_error("missing mixed");if(!mix){mix.reset(dynamic_cast<TH3*>(h->Clone(Form("H3_mixed_cent%d_mtMerged_phi%d",cent,p))));mix->SetDirectory(nullptr);}else mix->Add(h);}
  fo.cd();mix->Write();
  for(int q=0;q<3;++q){std::unique_ptr<TH3> same;
   for(int m=mtFirst;m<=mtLast;++m){auto*h=dynamic_cast<TH3*>(fi.Get(Form("H3_same_cent%d_mt%d_qn%d_phi%d",cent,m,q,p)));if(!h)throw std::runtime_error("missing same");if(!same){same.reset(dynamic_cast<TH3*>(h->Clone(Form("H3_same_cent%d_mtMerged_qn%d_phi%d",cent,q,p))));same->SetDirectory(nullptr);}else same->Add(h);}
   fo.cd();same->Write();
  }
 }
 fo.Close();
 const char* mtTag[3]={"0.40to0.50","0.50to0.60","0.60to0.70"};
 TFile fl(lambdaInput,"READ");if(fl.IsZombie())throw std::runtime_error("cannot open lambda input");
 std::unique_ptr<TH1> eff,wsum;
 for(int i=0;i<3;++i){auto*p=dynamic_cast<TH1*>(fl.Get(Form("hPurityKstar_mT%s_cent10to30",mtTag[i])));auto*w=dynamic_cast<TH1*>(fl.Get(Form("hWeightKstar_mT%s_cent10to30",mtTag[i])));if(!p||!w)throw std::runtime_error("missing lambda purity/weight");
  if(!eff){eff.reset(dynamic_cast<TH1*>(p->Clone("hPurityKstar_mT0.40to0.70_cent10to30")));eff->SetDirectory(nullptr);eff->Reset("ICES");wsum.reset(dynamic_cast<TH1*>(w->Clone("hWeightKstar_mT0.40to0.70_cent10to30")));wsum->SetDirectory(nullptr);wsum->Reset("ICES");}
  for(int b=1;b<=eff->GetNbinsX();++b){eff->SetBinContent(b,eff->GetBinContent(b)+p->GetBinContent(b)*w->GetBinContent(b));wsum->SetBinContent(b,wsum->GetBinContent(b)+w->GetBinContent(b));}
 }
 for(int b=1;b<=eff->GetNbinsX();++b)if(wsum->GetBinContent(b)>0)eff->SetBinContent(b,eff->GetBinContent(b)/wsum->GetBinContent(b));
 TFile flo(lambdaOutput,"RECREATE");eff->Write();wsum->Write();flo.Close();
 printf("wrote %s and %s\n",dataOutput.Data(),lambdaOutput.Data());
}
