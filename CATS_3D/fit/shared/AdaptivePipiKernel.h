#ifndef ADAPTIVE_PIPI_KERNEL_H
#define ADAPTIVE_PIPI_KERNEL_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

// Positive-q-octant |psi(q,r)|^2, folded with a symmetric Cartesian Gaussian.
class AdaptivePipiKernel {
public:
    explicit AdaptivePipiKernel(const std::string& fileName);
    std::vector<double> Fold(double rout, double rside, double rlong) const;
    std::vector<std::pair<unsigned, double>> AxisOverlaps(unsigned axis,
                                                           double lowMeV,
                                                           double highMeV) const;
    double QCenter(unsigned axis, unsigned i) const;
    size_t QIndex(unsigned ix, unsigned iy, unsigned iz) const;
    unsigned NumQ(unsigned axis) const;

private:
    std::vector<double> qEdges_[3];
    std::vector<double> rEdges_[3];
    std::vector<float> wavefunction2_;
    unsigned nQ_[3] = {0, 0, 0};
    unsigned nR_[3] = {0, 0, 0};
    size_t nR3_ = 0;
};

#endif
