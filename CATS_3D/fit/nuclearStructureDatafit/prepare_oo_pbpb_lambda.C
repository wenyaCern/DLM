#include "TFile.h"
#include "TH1.h"
#include <memory>
#include <stdexcept>

void prepare_oo_pbpb_lambda(
 const TString input="../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root",
 const TString output="../inputs/nuclearStructureDatafit/lambda_OOproxy_PbPb_cent10to30_mt0p2to0p7.root"){
 const char*tag[5]={"0.20to0.30","0.30to0.40","0.40to0.50","0.50to0.60","0.60to0.70"};
 TFile f(input,"READ");std::unique_ptr<TH1>e,wSum;
 for(auto t:tag){auto*p=dynamic_cast<TH1*>(f.Get(Form("hPurityKstar_mT%s_cent10to30",t)));auto*w=dynamic_cast<TH1*>(f.Get(Form("hWeightKstar_mT%s_cent10to30",t)));if(!p||!w)throw std::runtime_error("missing purity/weight");if(!e){e.reset(dynamic_cast<TH1*>(p->Clone("hPurityKstar_OOproxy_mt0.244to0.714_cent10to30")));e->SetDirectory(nullptr);e->Reset("ICES");wSum.reset(dynamic_cast<TH1*>(w->Clone("hWeightKstar_OOproxy_mt0.244to0.714_cent10to30")));wSum->SetDirectory(nullptr);wSum->Reset("ICES");}for(int b=1;b<=e->GetNbinsX();++b){e->SetBinContent(b,e->GetBinContent(b)+p->GetBinContent(b)*w->GetBinContent(b));wSum->SetBinContent(b,wSum->GetBinContent(b)+w->GetBinContent(b));}}
 for(int b=1;b<=e->GetNbinsX();++b)if(wSum->GetBinContent(b)>0)e->SetBinContent(b,e->GetBinContent(b)/wSum->GetBinContent(b));TFile o(output,"RECREATE");e->Write();wSum->Write();o.Close();printf("wrote %s\n",output.Data());
}
