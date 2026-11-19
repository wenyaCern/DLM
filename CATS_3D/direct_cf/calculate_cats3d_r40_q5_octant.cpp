#include "Basics3D.h"
#include "CATSconstants.h"

#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc,char** argv){
    if(argc<2||argc>4)throw std::runtime_error("Usage: calculate_cats3d_r40_q5_octant OUTPUT.csv [scale] [side|all]");
    const double scale=(argc>=3)?std::stod(argv[2]):1.0;
    const bool all=(argc==4 && std::string(argv[3])=="all");
    if(scale<=0)throw std::runtime_error("Scale must be positive");
    const std::array<double,3> radius=all
        ?std::array<double,3>{8.98989*scale,7.99913*scale,9.40678*scale}
        :std::array<double,3>{8.98989,7.99913*scale,9.40678};
    const double lambda=0.83*(1.0-0.14)*(1.0-0.14);
    CATS3D kitty;
    kitty.SetMomBinsOut(10,0.0,50.0);
    kitty.SetMomBinsSide(10,0.0,50.0);
    kitty.SetMomBinsLong(10,0.0,50.0);
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
    kitty.KillTheCat();
    std::ofstream out(argv[1]);
    if(!out)throw std::runtime_error("Cannot create output CSV");
    out<<"i_out,i_side,i_long,C\n";
    for(unsigned io=0;io<10;++io){
        for(unsigned is=0;is<10;++is){
            for(unsigned il=0;il<10;++il)
                out<<io<<','<<is<<','<<il<<','
                   <<1.0+lambda*(kitty.GetCorrFun(io,is,il)-1.0)<<'\n';
        }
        std::cerr<<"completed qout index "<<io+1<<" / 10\n";
    }
}
