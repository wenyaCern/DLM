#include "AdaptivePipiKernel.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    if(argc!=2) return 2;
    try{
        AdaptivePipiKernel kernel(argv[1]);
        const auto c=kernel.Fold(8.0,8.0,8.0);
        if(c.size()!=1 || !std::isfinite(c[0]) || c[0]<0.9 || c[0]>1.1)
            throw std::runtime_error("Unexpected folded correlation");
        const auto both=kernel.AxisOverlaps(0,-5,5);
        const auto positive=kernel.AxisOverlaps(0,0,5);
        if(both.size()!=1 || positive.size()!=1 ||
           std::abs(both[0].second-10)>1e-9 || std::abs(positive[0].second-5)>1e-9)
            throw std::runtime_error("Positive/negative q bin overlap failed");
        std::cout << "Folded C(q=2.5,2.5,2.5; R=8 fm)=" << c[0]
                  << "; symmetric q-bin overlap passed\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr << e.what() << '\n';
        return 1;
    }
}
