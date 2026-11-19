#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <array>
#include <iomanip>
#include <iostream>

namespace {
constexpr std::array<double,3> radius={8.98989,7.99913,9.40678};
constexpr double lambda=0.83*(1.0-0.14)*(1.0-0.14);

void Run(double rmax,unsigned ngrid){
    CATS3D kitty;
    kitty.SetMomBinsOut(2,0.0,20.0);
    kitty.SetMomBinsSide(2,0.0,20.0);
    kitty.SetMomBinsLong(2,0.0,20.0);
    kitty.SetGridX(ngrid,-rmax,rmax);
    kitty.SetGridY(ngrid,-rmax,rmax);
    kitty.SetGridZ(ngrid,-rmax,rmax);
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
    kitty.KillTheCat();
    double sum=0;
    for(unsigned ix=0;ix<2;++ix)
        for(unsigned iy=0;iy<2;++iy)
            for(unsigned iz=0;iz<2;++iz)
                sum+=1.0+lambda*(kitty.GetCorrFun(ix,iy,iz)-1.0);
    std::cout<<rmax<<','<<ngrid<<','<<sum/8<<std::endl;
}
}

int main(){
    std::cout<<std::setprecision(10)<<"rmax_fm,ngrid,C_center_10MeV\n";
    Run(40,31);
    Run(40,61);
    Run(60,91);
}
