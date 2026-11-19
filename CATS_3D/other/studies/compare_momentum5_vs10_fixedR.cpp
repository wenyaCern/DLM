#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <array>
#include <iomanip>
#include <iostream>

namespace {
constexpr std::array<double,3> radius={8.98989,7.99913,9.40678};
constexpr double lambda=0.83*(1.0-0.14)*(1.0-0.14);

void Configure(CATS3D& kitty,unsigned step){
    kitty.SetMomBinsOut(40/step,0.0,40.0);
    kitty.SetMomBinsSide(20/step,0.0,20.0);
    kitty.SetMomBinsLong(20/step,0.0,20.0);
    kitty.SetGridX(31,-40.0,40.0);
    kitty.SetGridY(31,-40.0,40.0);
    kitty.SetGridZ(31,-40.0,40.0);
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

void Run(unsigned step){
    CATS3D kitty;
    Configure(kitty,step);
    kitty.KillTheCat();
    const unsigned n=20/step;
    for(unsigned coarse=0;coarse<2;++coarse){
        double sum=0;
        for(unsigned ix=0;ix<n;++ix)
            for(unsigned iy=0;iy<n;++iy)
                for(unsigned iz=0;iz<n;++iz)
                    sum+=1.0+lambda*(kitty.GetCorrFun(coarse*n+ix,iy,iz)-1.0);
        std::cout<<step<<','<<(10+20*coarse)<<','<<sum/(n*n*n)<<'\n';
    }
}
}

int main(){
    std::cout<<std::setprecision(10)<<"step_MeV,qout_center_MeV,C_model\n";
    Run(10);
    Run(5);
}
