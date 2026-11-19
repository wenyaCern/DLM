#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <array>
#include <iomanip>
#include <iostream>

namespace {
constexpr std::array<double,3> radius={9.09009,8.08106,9.50016};
constexpr double lambda=0.85*(1.0-0.14)*(1.0-0.14);

void Configure(CATS3D& kitty,double rMax,bool fine){
    kitty.SetMomBinsOut(fine?6:3,0.0,60.0);
    kitty.SetMomBinsSide(fine?2:1,0.0,20.0);
    kitty.SetMomBinsLong(fine?2:1,0.0,20.0);
    kitty.SetGridX(31,-rMax,rMax);
    kitty.SetGridY(31,-rMax,rMax);
    kitty.SetGridZ(31,-rMax,rMax);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);
    CATSparameters pars(CATSparameters::tSource3D,3,true);
    for(unsigned a=0;a<3;++a)pars.SetParameter(a,radius[a]);
    kitty.SetAnaSource(Basics_Source_Gauss3D,pars);
    kitty.SetPdgId(211,211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);
}

void Run(double rMax){
    CATS3D coarseKitty;
    Configure(coarseKitty,rMax,false);
    coarseKitty.KillTheCat();
    CATS3D kitty;
    Configure(kitty,rMax,true);
    kitty.KillTheCat();
    for(unsigned coarse=0;coarse<3;++coarse){
        double c=0;
        for(unsigned dx=0;dx<2;++dx)
            for(unsigned dy=0;dy<2;++dy)
                for(unsigned dz=0;dz<2;++dz)
                    c+=(1.0+lambda*(kitty.GetCorrFun(2*coarse+dx,dy,dz)-1.0))/8.0;
        const double c20=1.0+lambda*(coarseKitty.GetCorrFun(coarse,0,0)-1.0);
        std::cout<<rMax<<','<<(20*coarse+10)<<','<<c20<<','<<c<<','<<c-c20<<'\n';
    }
}
}

int main(){
    std::cout<<std::setprecision(9)<<"rMax_fm,qOut_center_MeV,C20_center,C10_average,delta\n";
    Run(20.0);
    Run(40.0);
}
