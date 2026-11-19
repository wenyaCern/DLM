#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <array>
#include <iomanip>
#include <iostream>

namespace {
constexpr std::array<double,3> radius0={8.98989,7.99913,9.40678};
constexpr double weakScale=(1.0-0.14)*(1.0-0.14);
constexpr double data10=1.18292, data30=1.07473;
constexpr double baseline10=1.21702/1.218069201;
constexpr double baseline30=1.05322/1.054208303;

void Run(double scale){
    CATS3D kitty;
    kitty.SetMomBinsOut(4,0.0,40.0);
    kitty.SetMomBinsSide(2,0.0,20.0);
    kitty.SetMomBinsLong(2,0.0,20.0);
    kitty.SetGridX(31,-40.0,40.0);
    kitty.SetGridY(31,-40.0,40.0);
    kitty.SetGridZ(31,-40.0,40.0);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);
    CATSparameters pars(CATSparameters::tSource3D,3,true);
    for(unsigned a=0;a<3;++a)pars.SetParameter(a,radius0[a]*scale);
    kitty.SetAnaSource(Basics_Source_Gauss3D,pars);
    kitty.SetPdgId(211,211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);
    kitty.KillTheCat();
    double raw[2]={0,0};
    for(unsigned coarse=0;coarse<2;++coarse){
        for(unsigned ix=0;ix<2;++ix)
            for(unsigned iy=0;iy<2;++iy)
                for(unsigned iz=0;iz<2;++iz)
                    raw[coarse]+=kitty.GetCorrFun(coarse*2+ix,iy,iz)/8.0;
    }
    const double lambdaEff=(data10/baseline10-1)/(raw[0]-1);
    const double pred30=baseline30*(1+lambdaEff*(raw[1]-1));
    std::cout<<scale<<','<<raw[0]<<','<<raw[1]<<','<<lambdaEff/weakScale<<','<<pred30<<','<<pred30-data30<<'\n';
}
}

int main(){
    std::cout<<std::setprecision(9)<<"R_scale,Craw10,Craw30,lambda_gen_for_peak,pred_C30,minus_data30\n";
    for(double scale:{1.0,0.95,0.90,0.85,0.80,0.75})Run(scale);
}
