#include "TFile.h"
#include "TH2.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace {
struct Mean { double value=0., error=0.; };
Mean selectedMean(TH2* sum, TH2* count, int xlo, int xhi, int ylo, int yhi) {
  double s=0., s2=0., n=0.;
  for(int ix=xlo;ix<=xhi;++ix) for(int iy=ylo;iy<=yhi;++iy) {
    s += sum->GetBinContent(ix,iy);
    const double e=sum->GetBinError(ix,iy); s2 += e*e;
    n += count->GetBinContent(ix,iy);
  }
  Mean out; out.value=s/n;
  const double variance=(s2-s*s/n)/(n-1.);
  out.error=std::sqrt(std::max(0.,variance/n));
  return out;
}
}

void calculate_ft0a_resolution(
 const char* input="/Users/wuwenya/Downloads/AnalysisResults.root",
 const char* csvName="farm_output/ft0a_resolution_vs_q2.csv") {
 TFile f(input,"READ");
 auto*d=dynamic_cast<TDirectory*>(f.Get("femto-dream-producer-task/Event"));
 if(!d) throw std::runtime_error("missing producer Event directory");
 auto*ct=dynamic_cast<TH2*>(d->Get("epReso_FT0CTPC"));
 auto*at=dynamic_cast<TH2*>(d->Get("epReso_FT0ATPC"));
 auto*ca=dynamic_cast<TH2*>(d->Get("epReso_FT0CFT0A"));
 auto*n=dynamic_cast<TH2*>(d->Get("epReso_count"));
 if(!ct||!at||!ca||!n) throw std::runtime_error("missing resolution histograms");
 const int xlo[4]={1,2,4,6},xhi[4]={1,3,5,8};
 const char*cent[4]={"0-10","10-30","30-50","50-80"};
 const int qlo[3]={1,4,8},qhi[3]={3,7,10};
 std::ofstream csv(csvName);csv<<"centrality,q2_bin,N,mean_CTPC,err_CTPC,mean_ATPC,err_ATPC,mean_CA,err_CA,res_FT0A,err_resA_independent,res_FT0C,err_resC_independent\n";
 for(int c=0;c<4;++c){
  printf("CENT %s%%\n",cent[c]);
  for(int q=0;q<3;++q){
   const auto mct=selectedMean(ct,n,xlo[c],xhi[c],qlo[q],qhi[q]);
   const auto mat=selectedMean(at,n,xlo[c],xhi[c],qlo[q],qhi[q]);
   const auto mca=selectedMean(ca,n,xlo[c],xhi[c],qlo[q],qhi[q]);
   const double nn=n->Integral(xlo[c],xhi[c],qlo[q],qhi[q]);
   const double r=std::sqrt(mat.value*mca.value/mct.value);
   const double er=.5*r*std::sqrt(std::pow(mat.error/mat.value,2)+std::pow(mca.error/mca.value,2)+std::pow(mct.error/mct.value,2));
   const double rc=std::sqrt(mct.value*mca.value/mat.value);
   const double erc=.5*rc*std::sqrt(std::pow(mat.error/mat.value,2)+std::pow(mca.error/mca.value,2)+std::pow(mct.error/mct.value,2));
   printf(" q%d N=%.0f  CT=%.6f+/-%.6f AT=%.6f+/-%.6f CA=%.6f+/-%.6f  ResA=%.6f+/-%.6f ResC=%.6f+/-%.6f\n",q,nn,mct.value,mct.error,mat.value,mat.error,mca.value,mca.error,r,er,rc,erc);
   csv<<cent[c]<<','<<q<<','<<nn<<','<<mct.value<<','<<mct.error<<','<<mat.value<<','<<mat.error<<','<<mca.value<<','<<mca.error<<','<<r<<','<<er<<','<<rc<<','<<erc<<'\n';
  }
 }
}
