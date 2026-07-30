#include "CATS3D.h"

#include "gsl_sf_bessel.h"
#include "gsl_sf_coulomb.h"
#include "gsl_sf_legendre.h"
#include "gsl_sf_gamma.h"
#include <gsl/gsl_errno.h>

#include <cmath>
#include <cstring>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {
int CATS3DProfileLevel() {
    const char* env = std::getenv("CATS3D_PROFILE");
    return env ? std::atoi(env) : 0;
}

std::chrono::steady_clock::time_point CATS3DNow() {
    return std::chrono::steady_clock::now();
}

double CATS3DSecondsSince(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

// constructor
CATS3D::CATS3D()
    : RedMass(0.0),
      QuantumStatistics(-1),
      IdenticalParticles(false),
      Q1Q2(0),
      CurrentWfMode(wfNone),
      MaxPw(16),

      NumMomBins_out(0),
      NumMomBins_side(0),
      NumMomBins_long(0),
      MomBin_out(nullptr),
      MomBin_side(nullptr),
      MomBin_long(nullptr),
      MomBinCenter_out(nullptr),
      MomBinCenter_side(nullptr),
      MomBinCenter_long(nullptr),

      NumGridBinsX(0),
      NumGridBinsY(0),
      NumGridBinsZ(0),
      GridBinX(nullptr),
      GridBinY(nullptr),
      GridBinZ(nullptr),
      GridBinCenterX(nullptr),
      GridBinCenterY(nullptr),
      GridBinCenterZ(nullptr),
      GridStepX(0.0),
      GridStepY(0.0),
      GridStepZ(0.0),

      AnalyticSource(nullptr),
      ForwardedSource(nullptr),
      SourceContext(nullptr),
      AnaSourcePar(nullptr),
      NumSourcePars(0),

      ShortRangePotential3D(nullptr),
      PotPar3D(nullptr),

      kCorrFun(nullptr),
      kCorrFunErr(nullptr),

      SourceGridReady(false),
      SourceUpdated(false),
      ComputedWaveFunction(false),
      ComputedCorrFunction(false),
      UseFastKernel(false),
      FastKernelReady(false),
      FastKernelSubdivisions(1),
      BaseSourceGrid3D(nullptr),
      kSourceGrid3D(nullptr) {
    pdgID[0] = 0;
    pdgID[1] = 0;
}

// destructor
CATS3D::~CATS3D() {
    DeleteCorrFun();
    if (kSourceGrid3D) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            if (!kSourceGrid3D[uOut]) continue;
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                if (!kSourceGrid3D[uOut][uSide]) continue;
                for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                    if (kSourceGrid3D[uOut][uSide][uLong]) {
                        delete kSourceGrid3D[uOut][uSide][uLong];
                        kSourceGrid3D[uOut][uSide][uLong] = nullptr;
                    }
                }
                delete[] kSourceGrid3D[uOut][uSide];
                kSourceGrid3D[uOut][uSide] = nullptr;
            }
            delete[] kSourceGrid3D[uOut];
            kSourceGrid3D[uOut] = nullptr;
        }
        delete[] kSourceGrid3D;
        kSourceGrid3D = nullptr;
    }

    if (BaseSourceGrid3D) {
        delete BaseSourceGrid3D;
        BaseSourceGrid3D = nullptr;
    }

    DeleteGridBins();
    DeleteAllMomBinDefinitions();

    if (AnaSourcePar) {
        delete AnaSourcePar;
        AnaSourcePar = nullptr;
    }

    if (PotPar3D) {
        delete PotPar3D;
        PotPar3D = nullptr;
    }
}


// delete helpers
void CATS3D::DeleteAllMomBinDefinitions() {
    if (MomBin_out) {
        delete[] MomBin_out;
        MomBin_out = nullptr;
    }
    if (MomBin_side) {
        delete[] MomBin_side;
        MomBin_side = nullptr;
    }
    if (MomBin_long) {
        delete[] MomBin_long;
        MomBin_long = nullptr;
    }

    if (MomBinCenter_out) {
        delete[] MomBinCenter_out;
        MomBinCenter_out = nullptr;
    }
    if (MomBinCenter_side) {
        delete[] MomBinCenter_side;
        MomBinCenter_side = nullptr;
    }
    if (MomBinCenter_long) {
        delete[] MomBinCenter_long;
        MomBinCenter_long = nullptr;
    }

    NumMomBins_out = 0;
    NumMomBins_side = 0;
    NumMomBins_long = 0;
}

void CATS3D::DeleteGridBins() {
    if (GridBinX) {
        delete[] GridBinX;
        GridBinX = nullptr;
    }
    if (GridBinY) {
        delete[] GridBinY;
        GridBinY = nullptr;
    }
    if (GridBinZ) {
        delete[] GridBinZ;
        GridBinZ = nullptr;
    }

    if (GridBinCenterX) {
        delete[] GridBinCenterX;
        GridBinCenterX = nullptr;
    }
    if (GridBinCenterY) {
        delete[] GridBinCenterY;
        GridBinCenterY = nullptr;
    }
    if (GridBinCenterZ) {
        delete[] GridBinCenterZ;
        GridBinCenterZ = nullptr;
    }

    NumGridBinsX = 0;
    NumGridBinsY = 0;
    NumGridBinsZ = 0;

    GridStepX = 0.0;
    GridStepY = 0.0;
    GridStepZ = 0.0;
}

void CATS3D::DeleteCorrFun() {
    if (kCorrFun) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            if (!kCorrFun[uOut]) continue;
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                if (kCorrFun[uOut][uSide]) {
                    delete[] kCorrFun[uOut][uSide];
                    kCorrFun[uOut][uSide] = nullptr;
                }
            }
            delete[] kCorrFun[uOut];
            kCorrFun[uOut] = nullptr;
        }
        delete[] kCorrFun;
        kCorrFun = nullptr;
    }

    if (kCorrFunErr) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            if (!kCorrFunErr[uOut]) continue;
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                if (kCorrFunErr[uOut][uSide]) {
                    delete[] kCorrFunErr[uOut][uSide];
                    kCorrFunErr[uOut][uSide] = nullptr;
                }
            }
            delete[] kCorrFunErr[uOut];
            kCorrFunErr[uOut] = nullptr;
        }
        delete[] kCorrFunErr;
        kCorrFunErr = nullptr;
    }
}

void CATS3D::AllocateCorrFun() {
    DeleteCorrFun();

    kCorrFun = new double**[NumMomBins_out];
    kCorrFunErr = new double**[NumMomBins_out];

    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        kCorrFun[uOut] = new double*[NumMomBins_side];
        kCorrFunErr[uOut] = new double*[NumMomBins_side];

        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            kCorrFun[uOut][uSide] = new double[NumMomBins_long];
            kCorrFunErr[uOut][uSide] = new double[NumMomBins_long];

            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                kCorrFun[uOut][uSide][uLong] = 0.0;
                kCorrFunErr[uOut][uSide][uLong] = 0.0;
            }
        }
    }
}

void CATS3D::DeleteWaveFunctionCache() {
    FastWaveFunction2.clear();
    FastKernelReady = false;
}

unsigned CATS3D::GetFlatMomIndex(const unsigned& uBinOut,
                                 const unsigned& uBinSide,
                                 const unsigned& uBinLong) const {
    return (uBinOut * NumMomBins_side + uBinSide) * NumMomBins_long + uBinLong;
}

unsigned CATS3D::GetFlatGridIndex(const unsigned& uBinX,
                                  const unsigned& uBinY,
                                  const unsigned& uBinZ) const {
    return (uBinX * NumGridBinsY + uBinY) * NumGridBinsZ + uBinZ;
}

// setter & getter
void CATS3D::SetRedMass(const double& redMass) {
    if (RedMass == redMass) return;
    RedMass = redMass;

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

double CATS3D::GetRedMass() const {
    return RedMass;
}

void CATS3D::SetPdgId(const int& id1, const int& id2) {
    if (pdgID[0] == id1 && pdgID[1] == id2) return;

    pdgID[0] = id1;
    pdgID[1] = id2;

    // same logic as original CATS:
    // if QS was manually set to 0/1, do not override IdenticalParticles here
    if (QuantumStatistics != 0 && QuantumStatistics != 1) {
        IdenticalParticles = (id1 == id2);
    }

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false; 
}

void CATS3D::GetPdgId(int& id1, int& id2) const {
    id1 = pdgID[0];
    id2 = pdgID[1];
}

void CATS3D::SetQuantumStatistics(short qs) {
    // keep the CATS convention:
    // 0 = off, 1 = on, otherwise auto
    if (qs != 0 && qs != 1) qs = -1;
    if (QuantumStatistics == qs) return;

    QuantumStatistics = qs;

    bool identical = false;
    switch (QuantumStatistics) {
        case 0:
            identical = false;
            break;
        case 1:
            identical = true;
            break;
        default:
            identical = (pdgID[0] == pdgID[1]);
            break;
    }

    if (IdenticalParticles != identical) {
        IdenticalParticles = identical;
        DeleteWaveFunctionCache();
        ComputedWaveFunction = false;
        ComputedCorrFunction = false;
    }
}

short CATS3D::GetQuantumStatistics() const {
    return QuantumStatistics;
}

void CATS3D::SetQ1Q2(const int& q1q2) {
    if (Q1Q2 == q1q2) return;
    Q1Q2 = q1q2;

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

int CATS3D::GetQ1Q2() const {
    return Q1Q2;
}

void CATS3D::SetMaxPw(const unsigned short& maxPw) {
    const unsigned short safeMaxPw = maxPw ? maxPw : 1;
    if (MaxPw == safeMaxPw) return;

    MaxPw = safeMaxPw;

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

unsigned short CATS3D::GetMaxPw() const {
    return MaxPw;
}

void CATS3D::GetMomentumVector(const unsigned& uBinOut,
                               const unsigned& uBinSide,
                               const unsigned& uBinLong,
                               double& qOut,
                               double& qSide,
                               double& qLong) const {
    qOut = 0.0;
    qSide = 0.0;
    qLong = 0.0;

    if (uBinOut >= NumMomBins_out) return;
    if (uBinSide >= NumMomBins_side) return;
    if (uBinLong >= NumMomBins_long) return;

    qOut  = GetMomentumOut(uBinOut);
    qSide = GetMomentumSide(uBinSide);
    qLong = GetMomentumLong(uBinLong);
}

double CATS3D::GetQNorm(const double& qOut,
                        const double& qSide,
                        const double& qLong) const {
    return std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);
}

double CATS3D::GetRNorm(const double& x,
                        const double& y,
                        const double& z) const {
    return std::sqrt(x*x + y*y + z*z);
}

double CATS3D::GetCosTheta_qr(const double& qOut,
                              const double& qSide,
                              const double& qLong,
                              const double& x,
                              const double& y,
                              const double& z) const {
    const double qNorm = GetQNorm(qOut, qSide, qLong);
    const double rNorm = GetRNorm(x, y, z);

    if (qNorm == 0.0 || rNorm == 0.0) {
        return 1.0;
    }

    const double qDotR = qOut*x + qSide*y + qLong*z;
    double cosTheta = qDotR / (qNorm * rNorm);

    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;

    return cosTheta;
}

double CATS3D::GetRho(const double& q,
                      const double& r) const {
    return q * r * FmToNu;
}

double CATS3D::GetEta(const double& q) const {
    if (q == 0.0) {
        return 0.0;
    }
    return RedMass * double(Q1Q2) * AlphaFS / q;
}

void CATS3D::SetWfMode(const WfMode& mode) {
    if (CurrentWfMode == mode) return;
    CurrentWfMode = mode;

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false; 
}

CATS3D::WfMode CATS3D::GetWfMode() const {
    return CurrentWfMode;
}

bool CATS3D::GetIdenticalParticles() const {
    return IdenticalParticles;
}

void CATS3D::SetFastKernelMode(const bool& useFast) {
    if (UseFastKernel == useFast) return;
    UseFastKernel = useFast;
    DeleteWaveFunctionCache();
    DeleteCorrFun();
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

bool CATS3D::GetFastKernelMode() const {
    return UseFastKernel;
}

void CATS3D::SetFastKernelSubdivisions(const unsigned& subdivisions) {
    const unsigned safeSubdivisions = subdivisions ? subdivisions : 1;
    if (FastKernelSubdivisions == safeSubdivisions) return;
    FastKernelSubdivisions = safeSubdivisions;
    DeleteWaveFunctionCache();
    DeleteCorrFun();
    ComputedCorrFunction = false;
}

unsigned CATS3D::GetFastKernelSubdivisions() const {
    return FastKernelSubdivisions;
}

void CATS3D::BuildFastKernelOnly() {
    BuildFastWaveFunctionKernel();
}

bool CATS3D::IsFastKernelReady() const {
    return FastKernelReady && !FastWaveFunction2.empty();
}

size_t CATS3D::GetFastKernelMomSize() const {
    return size_t(NumMomBins_out) * size_t(NumMomBins_side) * size_t(NumMomBins_long);
}

size_t CATS3D::GetFastKernelGridSize() const {
    const unsigned subdivisions = FastKernelSubdivisions ? FastKernelSubdivisions : 1;
    const size_t nBaseGrid =
        size_t(NumGridBinsX) * size_t(NumGridBinsY) * size_t(NumGridBinsZ);
    const size_t nSubCells =
        size_t(subdivisions) * size_t(subdivisions) * size_t(subdivisions);
    return nBaseGrid * nSubCells;
}

size_t CATS3D::GetFastKernelStorageSize() const {
    return GetFastKernelMomSize() * GetFastKernelGridSize();
}

const std::vector<float>& CATS3D::GetFastWaveFunctionKernel() const {
    return FastWaveFunction2;
}

void CATS3D::SetFastWaveFunctionKernel(const std::vector<float>& kernel,
                                       const unsigned& subdivisions) {
    const unsigned safeSubdivisions = subdivisions ? subdivisions : 1;
    if (FastKernelSubdivisions != safeSubdivisions) {
        FastKernelSubdivisions = safeSubdivisions;
        DeleteCorrFun();
    }

    const size_t expectedSize = GetFastKernelStorageSize();
    if (expectedSize == 0) {
        throw std::runtime_error("CATS3D::SetFastWaveFunctionKernel: momentum/source grids are not configured.");
    }
    if (kernel.size() != expectedSize) {
        throw std::runtime_error("CATS3D::SetFastWaveFunctionKernel: kernel size does not match current q/r grid.");
    }

    FastWaveFunction2 = kernel;
    UseFastKernel = true;
    FastKernelReady = true;
    ComputedWaveFunction = true;
    DeleteCorrFun();
    ComputedCorrFunction = false;
}

// setter & getter of momentums
void CATS3D::SetMomBinsOut(const unsigned& nBins, const double& minMom, const double& maxMom) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetMomBinsOut: nBins cannot be zero.");
    }
    if (minMom > maxMom) {
        throw std::runtime_error("CATS3D::SetMomBinsOut: minMom cannot be larger than maxMom.");
    }
    if (minMom == maxMom && nBins != 1) {
        throw std::runtime_error("CATS3D::SetMomBinsOut: minMom == maxMom is only allowed for nBins == 1.");
    }

    const double binWidth = (maxMom - minMom) / double(nBins);

    bool sameBinning = true;
    if (!MomBin_out || NumMomBins_out != nBins) {
        sameBinning = false;
    } else {
        if (MomBin_out[0] != minMom || MomBin_out[nBins] != maxMom) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (MomBin_out) delete[] MomBin_out;
    if (MomBinCenter_out) delete[] MomBinCenter_out;

    MomBin_out = new double[nBins + 1];
    MomBinCenter_out = new double[nBins];
    NumMomBins_out = nBins;

    for (unsigned uBin = 0; uBin <= NumMomBins_out; ++uBin) {
        MomBin_out[uBin] = minMom + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumMomBins_out; ++uBin) {
        MomBinCenter_out[uBin] = 0.5 * (MomBin_out[uBin] + MomBin_out[uBin + 1]);
    }

    DeleteCorrFun();
    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false; 
}

void CATS3D::SetMomBinsSide(const unsigned& nBins, const double& minMom, const double& maxMom) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetMomBinsSide: nBins cannot be zero.");
    }
    if (minMom > maxMom) {
        throw std::runtime_error("CATS3D::SetMomBinsSide: minMom cannot be larger than maxMom.");
    }
    if (minMom == maxMom && nBins != 1) {
        throw std::runtime_error("CATS3D::SetMomBinsSide: minMom == maxMom is only allowed for nBins == 1.");
    }

    const double binWidth = (maxMom - minMom) / double(nBins);

    bool sameBinning = true;
    if (!MomBin_side || NumMomBins_side != nBins) {
        sameBinning = false;
    } else {
        if (MomBin_side[0] != minMom || MomBin_side[nBins] != maxMom) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (MomBin_side) delete[] MomBin_side;
    if (MomBinCenter_side) delete[] MomBinCenter_side;

    MomBin_side = new double[nBins + 1];
    MomBinCenter_side = new double[nBins];
    NumMomBins_side = nBins;

    for (unsigned uBin = 0; uBin <= NumMomBins_side; ++uBin) {
        MomBin_side[uBin] = minMom + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumMomBins_side; ++uBin) {
        MomBinCenter_side[uBin] = 0.5 * (MomBin_side[uBin] + MomBin_side[uBin + 1]);
    }

    DeleteCorrFun();
    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS3D::SetMomBinsLong(const unsigned& nBins, const double& minMom, const double& maxMom) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetMomBinsLong: nBins cannot be zero.");
    }
    if (minMom > maxMom) {
        throw std::runtime_error("CATS3D::SetMomBinsLong: minMom cannot be larger than maxMom.");
    }
    if (minMom == maxMom && nBins != 1) {
        throw std::runtime_error("CATS3D::SetMomBinsLong: minMom == maxMom is only allowed for nBins == 1.");
    }

    const double binWidth = (maxMom - minMom) / double(nBins);

    bool sameBinning = true;
    if (!MomBin_long || NumMomBins_long != nBins) {
        sameBinning = false;
    } else {
        if (MomBin_long[0] != minMom || MomBin_long[nBins] != maxMom) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (MomBin_long) delete[] MomBin_long;
    if (MomBinCenter_long) delete[] MomBinCenter_long;

    MomBin_long = new double[nBins + 1];
    MomBinCenter_long = new double[nBins];
    NumMomBins_long = nBins;

    for (unsigned uBin = 0; uBin <= NumMomBins_long; ++uBin) {
        MomBin_long[uBin] = minMom + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumMomBins_long; ++uBin) {
        MomBinCenter_long[uBin] = 0.5 * (MomBin_long[uBin] + MomBin_long[uBin + 1]);
    }

    DeleteCorrFun();
    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

unsigned CATS3D::GetNumMomBinsOut() const {
    return NumMomBins_out;
}

unsigned CATS3D::GetNumMomBinsSide() const {
    return NumMomBins_side;
}

unsigned CATS3D::GetNumMomBinsLong() const {
    return NumMomBins_long;
}

double CATS3D::GetMomentumOut(const unsigned& uBin) const {
    if (uBin >= NumMomBins_out || !MomBinCenter_out) return 0.0;
    return MomBinCenter_out[uBin];
}

double CATS3D::GetMomentumSide(const unsigned& uBin) const {
    if (uBin >= NumMomBins_side || !MomBinCenter_side) return 0.0;
    return MomBinCenter_side[uBin];
}

double CATS3D::GetMomentumLong(const unsigned& uBin) const {
    if (uBin >= NumMomBins_long || !MomBinCenter_long) return 0.0;
    return MomBinCenter_long[uBin];
}

double CATS3D::GetMomBinLowEdgeOut(const unsigned& uBin) const {
    if (uBin >= NumMomBins_out || !MomBin_out) return 0.0;
    return MomBin_out[uBin];
}

double CATS3D::GetMomBinUpEdgeOut(const unsigned& uBin) const {
    if (uBin >= NumMomBins_out || !MomBin_out) return 0.0;
    return MomBin_out[uBin + 1];
}

double CATS3D::GetMomBinLowEdgeSide(const unsigned& uBin) const {
    if (uBin >= NumMomBins_side || !MomBin_side) return 0.0;
    return MomBin_side[uBin];
}

double CATS3D::GetMomBinUpEdgeSide(const unsigned& uBin) const {
    if (uBin >= NumMomBins_side || !MomBin_side) return 0.0;
    return MomBin_side[uBin + 1];
}

double CATS3D::GetMomBinLowEdgeLong(const unsigned& uBin) const {
    if (uBin >= NumMomBins_long || !MomBin_long) return 0.0;
    return MomBin_long[uBin];
}

double CATS3D::GetMomBinUpEdgeLong(const unsigned& uBin) const {
    if (uBin >= NumMomBins_long || !MomBin_long) return 0.0;
    return MomBin_long[uBin + 1];
}

double CATS3D::GetTotMomentum3D(const unsigned& uBinOut,
                                const unsigned& uBinSide,
                                const unsigned& uBinLong) const {
    double qOut, qSide, qLong;
    GetMomentumVector(uBinOut, uBinSide, uBinLong, qOut, qSide, qLong);
    return GetQNorm(qOut, qSide, qLong);
}

// setter & getter of 3d source grid
void CATS3D::SetGridX(const unsigned& nBins, const double& minX, const double& maxX) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetGridX: nBins cannot be zero.");
    }
    if (minX > maxX) {
        throw std::runtime_error("CATS3D::SetGridX: minX cannot be larger than maxX.");
    }
    if (minX == maxX && nBins != 1) {
        throw std::runtime_error("CATS3D::SetGridX: minX == maxX is only allowed for nBins == 1.");
    }

    const double binWidth = (maxX - minX) / double(nBins);

    bool sameBinning = true;
    if (!GridBinX || NumGridBinsX != nBins) {
        sameBinning = false;
    } else {
        if (GridBinX[0] != minX || GridBinX[nBins] != maxX) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (GridBinX) delete[] GridBinX;
    if (GridBinCenterX) delete[] GridBinCenterX;

    GridBinX = new double[nBins + 1];
    GridBinCenterX = new double[nBins];
    NumGridBinsX = nBins;
    GridStepX = binWidth;

    for (unsigned uBin = 0; uBin <= NumGridBinsX; ++uBin) {
        GridBinX[uBin] = minX + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumGridBinsX; ++uBin) {
        GridBinCenterX[uBin] = 0.5 * (GridBinX[uBin] + GridBinX[uBin + 1]);
    }

    DeleteWaveFunctionCache();
    DeleteCorrFun();
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

void CATS3D::SetGridY(const unsigned& nBins, const double& minY, const double& maxY) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetGridY: nBins cannot be zero.");
    }
    if (minY > maxY) {
        throw std::runtime_error("CATS3D::SetGridY: minY cannot be larger than maxY.");
    }
    if (minY == maxY && nBins != 1) {
        throw std::runtime_error("CATS3D::SetGridY: minY == maxY is only allowed for nBins == 1.");
    }

    const double binWidth = (maxY - minY) / double(nBins);

    bool sameBinning = true;
    if (!GridBinY || NumGridBinsY != nBins) {
        sameBinning = false;
    } else {
        if (GridBinY[0] != minY || GridBinY[nBins] != maxY) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (GridBinY) delete[] GridBinY;
    if (GridBinCenterY) delete[] GridBinCenterY;

    GridBinY = new double[nBins + 1];
    GridBinCenterY = new double[nBins];
    NumGridBinsY = nBins;
    GridStepY = binWidth;

    for (unsigned uBin = 0; uBin <= NumGridBinsY; ++uBin) {
        GridBinY[uBin] = minY + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumGridBinsY; ++uBin) {
        GridBinCenterY[uBin] = 0.5 * (GridBinY[uBin] + GridBinY[uBin + 1]);
    }

    DeleteWaveFunctionCache();
    DeleteCorrFun();
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

void CATS3D::SetGridZ(const unsigned& nBins, const double& minZ, const double& maxZ) {
    if (!nBins) {
        throw std::runtime_error("CATS3D::SetGridZ: nBins cannot be zero.");
    }
    if (minZ > maxZ) {
        throw std::runtime_error("CATS3D::SetGridZ: minZ cannot be larger than maxZ.");
    }
    if (minZ == maxZ && nBins != 1) {
        throw std::runtime_error("CATS3D::SetGridZ: minZ == maxZ is only allowed for nBins == 1.");
    }

    const double binWidth = (maxZ - minZ) / double(nBins);

    bool sameBinning = true;
    if (!GridBinZ || NumGridBinsZ != nBins) {
        sameBinning = false;
    } else {
        if (GridBinZ[0] != minZ || GridBinZ[nBins] != maxZ) {
            sameBinning = false;
        }
    }
    if (sameBinning) return;

    if (GridBinZ) delete[] GridBinZ;
    if (GridBinCenterZ) delete[] GridBinCenterZ;

    GridBinZ = new double[nBins + 1];
    GridBinCenterZ = new double[nBins];
    NumGridBinsZ = nBins;
    GridStepZ = binWidth;

    for (unsigned uBin = 0; uBin <= NumGridBinsZ; ++uBin) {
        GridBinZ[uBin] = minZ + double(uBin) * binWidth;
    }
    for (unsigned uBin = 0; uBin < NumGridBinsZ; ++uBin) {
        GridBinCenterZ[uBin] = 0.5 * (GridBinZ[uBin] + GridBinZ[uBin + 1]);
    }

    DeleteWaveFunctionCache();
    DeleteCorrFun();
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

unsigned CATS3D::GetNumGridBinsX() const {
    return NumGridBinsX;
}

unsigned CATS3D::GetNumGridBinsY() const {
    return NumGridBinsY;
}

unsigned CATS3D::GetNumGridBinsZ() const {
    return NumGridBinsZ;
}

double CATS3D::GetGridCenterX(const unsigned& uBin) const {
    if (uBin >= NumGridBinsX || !GridBinCenterX) return 0.0;
    return GridBinCenterX[uBin];
}

double CATS3D::GetGridCenterY(const unsigned& uBin) const {
    if (uBin >= NumGridBinsY || !GridBinCenterY) return 0.0;
    return GridBinCenterY[uBin];
}

double CATS3D::GetGridCenterZ(const unsigned& uBin) const {
    if (uBin >= NumGridBinsZ || !GridBinCenterZ) return 0.0;
    return GridBinCenterZ[uBin];
}

double CATS3D::GetGridStepX() const {
    return GridStepX;
}

double CATS3D::GetGridStepY() const {
    return GridStepY;
}

double CATS3D::GetGridStepZ() const {
    return GridStepZ;
} 

// setter & getter of analytic source evaluation
void CATS3D::SetAnaSource(double (*AS)(double*), CATSparameters& Pars) {
    if (!AS) {
        throw std::runtime_error("CATS3D::SetAnaSource: source function pointer is null.");
    }

    AnalyticSource = AS;
    ForwardedSource = nullptr;
    SourceContext = nullptr;

    if (AnaSourcePar) {
        delete AnaSourcePar;
        AnaSourcePar = nullptr;
    }

    AnaSourcePar = new CATSparameters(Pars);
    NumSourcePars = Pars.GetNumPars();

    DeleteWaveFunctionCache();
    DeleteCorrFun();

    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

void CATS3D::SetAnaSource(double (*FS)(void*, double*),
                          void* context,
                          const unsigned& numparameters) {
    if (!FS) {
        throw std::runtime_error("CATS3D::SetAnaSource: forwarded source function pointer is null.");
    }
    if (!context) {
        throw std::runtime_error("CATS3D::SetAnaSource: forwarded source context is null.");
    }

    if (ForwardedSource == FS &&
        SourceContext == context &&
        AnaSourcePar &&
        NumSourcePars == numparameters) {
        return;
    }

    AnalyticSource = nullptr;
    ForwardedSource = FS;
    SourceContext = context;

    if (AnaSourcePar) {
        delete AnaSourcePar;
        AnaSourcePar = nullptr;
    }

    AnaSourcePar = new CATSparameters(CATSparameters::tSource3D,
                                      numparameters,
                                      true);
    NumSourcePars = numparameters;

    DeleteWaveFunctionCache();
    DeleteCorrFun();

    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

void CATS3D::SetAnaSource(const unsigned& WhichPar, const double& Value, const bool& SmallChange) {
    if (!AnaSourcePar) {
        throw std::runtime_error("CATS3D::SetAnaSource(par,val): analytic source parameters are not set.");
    }

    if (WhichPar >= NumSourcePars) {
        throw std::runtime_error("CATS3D::SetAnaSource(par,val): parameter index out of range.");
    }

    if (AnaSourcePar->GetParameter(WhichPar) == Value) return;

    AnaSourcePar->SetParameter(WhichPar, Value, false);

    if (!SmallChange) {
        DeleteWaveFunctionCache();
        DeleteCorrFun();
    }
    else {
        DeleteCorrFun();
    }

    SourceGridReady = SourceGridReady ? SmallChange : false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}

double CATS3D::GetAnaSourcePar(const unsigned& WhichPar) const {
    if (!AnaSourcePar) {
        return 0.0;
    }
    if (WhichPar >= NumSourcePars) {
        return 0.0;
    }
    return AnaSourcePar->GetParameter(WhichPar);
}

double CATS3D::EvaluateTheSource(CATSparameters* Pars) const{
    if((!AnalyticSource && !ForwardedSource) || !Pars){
        printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported a crash! This is most likely a bug!\n");
        printf("         Please contact the developers, reporting that source or Pars=NULL\n");
        return 0;
    }
    if(ForwardedSource){
        return ForwardedSource(SourceContext, Pars->GetParameters());
    }
    if(AnalyticSource){
        return AnalyticSource(Pars->GetParameters());
    }

    return 0;
}

double CATS3D::EvaluateTheSource(const double& qOut,
                                 const double& qSide,
                                 const double& qLong,
                                 const double& rOut,
                                 const double& rSide,
                                 const double& rLong) const {
    if ((!AnalyticSource && !ForwardedSource) || !AnaSourcePar) {
        return 0.0;
    }

    AnaSourcePar->SetVariable(0, qOut, false);
    AnaSourcePar->SetVariable(1, qSide, false);
    AnaSourcePar->SetVariable(2, qLong, false);
    AnaSourcePar->SetVariable(3, rOut, false);
    AnaSourcePar->SetVariable(4, rSide, false);
    AnaSourcePar->SetVariable(5, rLong, false);

    double* Pars = AnaSourcePar->GetParameters();
    if (ForwardedSource) {
        return ForwardedSource(SourceContext, Pars);
    }
    return AnalyticSource(Pars);
}

// status query
bool CATS3D::SourceStatus() const {
    return (AnalyticSource != nullptr || ForwardedSource != nullptr) &&
           AnaSourcePar != nullptr;
}

bool CATS3D::PotentialStatus() const {
    switch (CurrentWfMode) {
        case wfNone:
            return false;
        case wfPlaneWave:
            return true;
        case wfCoulomb:
            return true;
        case wfExternal3D:
            return false;
        case wfSchrodinger3D:
            return HasShortRangePotential3D();
        default:
            return false;
    }
}


bool CATS3D::CkStatus() const {
    return ComputedCorrFunction && kCorrFun != nullptr;
}

double CATS3D::CoulombPotential(const double& Radius) const {
    // Define of AlphaFS is in CATSconstants
    // Radius is expected in natural units (1/MeV)
    return Q1Q2 * AlphaFS / (std::fabs(Radius) + 1e-64); 
}

double CATS3D::GetCorrFun(const unsigned& uBinOut,
                          const unsigned& uBinSide,
                          const unsigned& uBinLong) const {
    if (!kCorrFun) return 0.0;
    if (uBinOut >= NumMomBins_out) return 0.0;
    if (uBinSide >= NumMomBins_side) return 0.0;
    if (uBinLong >= NumMomBins_long) return 0.0;
    return kCorrFun[uBinOut][uBinSide][uBinLong];
}

double CATS3D::GetCorrFunErr(const unsigned& uBinOut,
                             const unsigned& uBinSide,
                             const unsigned& uBinLong) const {
    if (!kCorrFunErr) return 0.0;
    if (uBinOut >= NumMomBins_out) return 0.0;
    if (uBinSide >= NumMomBins_side) return 0.0;
    if (uBinLong >= NumMomBins_long) return 0.0;
    return kCorrFunErr[uBinOut][uBinSide][uBinLong];
}

// Wave Function !
double CATS3D::PlanePartialWave3D(const double& r,
                                  const double& q,
                                  const unsigned short& l) const {
    const double rho = GetRho(q, r);

    if (rho == 0.0) {
        return 0.0;
    }

    // keep the same convention as 1D CATS:
    // u_l(r) = r * j_l(qr)
    if (rho > 0.0) {
        return r * FmToNu * gsl_sf_bessel_jl(l, rho);
    } else {
        return std::pow(-1.0, int(l)) * r * FmToNu * gsl_sf_bessel_jl(l, -rho);
    }
}

double CATS3D::CoulombPartialWave3D(const double& r,
                                    const double& q,
                                    const unsigned short& l) const {
    if (q == 0.0) {
        return 0.0;
    }

    const double eta = GetEta(q);
    const double rho = GetRho(q, r);

    if (rho == 0.0) {
        return 0.0;
    }

    double result = 0.0;
    double overflow = 0.0;

    gsl_sf_coulomb_wave_F_array(
        l,
        0,
        eta,
        std::fabs(rho),
        &result,
        &overflow
    );

    result /= q;

    // same parity handling as in 1D CATS
    if (rho < 0.0 && (l % 2 == 1)) {
        result = -result;
    }

    return result;
}

double CATS3D::CoulombPhaseShift3D(const unsigned short& l,
                                   const double& eta) const {
    // For Q1Q2 = 0 or eta = 0, Coulomb phase shift vanishes
    if (eta == 0.0) {
        return 0.0;
    }

    gsl_sf_result lnr;
    gsl_sf_result arg;

    const int status = gsl_sf_lngamma_complex_e(
        double(l) + 1.0,
        eta,
        &lnr,
        &arg
    );

    if (status != GSL_SUCCESS) {
        return 0.0;
    }

    return arg.val;
}

double CATS3D::ReferencePartialWave3D(const double& r,
                                      const double& q,
                                      const unsigned short& l) const {
    return Q1Q2 ? CoulombPartialWave3D(r, q, l)
                : PlanePartialWave3D(r, q, l);
}

double CATS3D::LegendrePl3D(const unsigned short& l,
                            const double& cosTheta) const {
    double x = cosTheta;
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;

    return gsl_sf_legendre_Pl(l, x);
}
    
// Could be 1+cos(2qDotr) but we return e^{iqDotr} here for uniform
std::complex<double> CATS3D::EvaluatePlaneWaveAmplitude(
    const double& qOut,
    const double& qSide,
    const double& qLong,
    const double& x,
    const double& y,
    const double& z) const {

    const double qDotR =
        qOut  * x * FmToNu +
        qSide * y * FmToNu +
        qLong * z * FmToNu;

    return std::complex<double>(std::cos(qDotR/2.), std::sin(qDotR/2.));
}

// phase-shift from Coulomb (CoulombPhaseShift3D()) included
std::complex<double> CATS3D::EvaluateCoulombPWAmplitude(
    const double& qOut,
    const double& qSide,
    const double& qLong,
    const double& x,
    const double& y,
    const double& z) const {

    const double q = 1./2.*GetQNorm(qOut, qSide, qLong);
    const double r = GetRNorm(x, y, z);
    const double eta = GetEta(q);

    // avoid singular points for now
    if (q == 0.0 || r == 0.0) {
        return std::complex<double>(1.0, 0.0);
    }

    const double cosTheta = GetCosTheta_qr(qOut, qSide, qLong, x, y, z);

    std::complex<double> psi(0.0, 0.0);

    for (unsigned short l = 0; l <= MaxPw; ++l) {
        // identical spin-0 bosons -> only even l
        if (IdenticalParticles && (l % 2 == 1)) continue;

        const double ul = ReferencePartialWave3D(r, q, l);
        const double Pl = LegendrePl3D(l, cosTheta);

        // radial part R_l = u_l / r_natural
        const double rNu = r * FmToNu;
        if (rNu == 0.0) continue;

        const double Rl = ul / rNu;

        // i^l
        std::complex<double> il(1.0, 0.0);
        switch (l % 4) {
            case 0: il = std::complex<double>(1.0, 0.0); break;
            case 1: il = std::complex<double>(0.0, 1.0); break;
            case 2: il = std::complex<double>(-1.0, 0.0); break;
            case 3: il = std::complex<double>(0.0, -1.0); break;
        }

        psi += double(2 * l + 1) * il * Rl * Pl;

        // const double sigmaL = CoulombPhaseShift3D(l, eta);
        // const std::complex<double> phaseSigma(std::cos(sigmaL), std::sin(sigmaL));
        // TODO: Try to understand the physics of the phase-shift !!!
        // Now we don't have phase transition
         

    }

    // if (IdenticalParticles) psi *= std::sqrt(2.0);

    return psi;
}

double CATS3D::EvaluateWaveFunction2(const unsigned& uBinOut,
                                     const unsigned& uBinSide,
                                     const unsigned& uBinLong,
                                     const double& x,
                                     const double& y,
                                     const double& z) const {
    if (uBinOut >= NumMomBins_out) return 0.0;
    if (uBinSide >= NumMomBins_side) return 0.0;
    if (uBinLong >= NumMomBins_long) return 0.0;

    double qOut = 0.0;
    double qSide = 0.0;
    double qLong = 0.0;
    GetMomentumVector(uBinOut, uBinSide, uBinLong, qOut, qSide, qLong);

    std::complex<double> psi(0.0, 0.0);

    switch (CurrentWfMode) {
        case wfPlaneWave: {
            psi = EvaluatePlaneWaveAmplitude(qOut, qSide, qLong, x, y, z);

            if (IdenticalParticles) {
                const std::complex<double> psiEx =
                    EvaluatePlaneWaveAmplitude(qOut, qSide, qLong, -x, -y, -z);
                psi = (psi + psiEx) / std::sqrt(2.0);
            }

            return std::norm(psi);
        }

        case wfCoulomb: {
            // For Coulomb PW, QS for identical spin-0 bosons
            // is already implemented in the PW sum by keeping only even l.
            psi = EvaluateCoulombPWAmplitude(qOut, qSide, qLong, x, y, z);

            if (IdenticalParticles) {
                const std::complex<double> psiEx =
                    EvaluateCoulombPWAmplitude(qOut, qSide, qLong, -x, -y, -z);
                psi = (psi + psiEx) / std::sqrt(2.0);
            }
            return std::norm(psi);
        }

        case wfSchrodinger3D: {
            psi = EvaluateShortRangePWAmplitude(qOut, qSide, qLong, x, y, z);
            return std::norm(psi);
        }

        case wfExternal3D:
        case wfNone:
        default:
            return 0.0;
    }
}

// short-range potential and wave-function propogation
std::complex<double> CATS3D::EvaluateShortRangePWAmplitude(
    const double& qOut,
    const double& qSide,
    const double& qLong,
    const double& x,
    const double& y,
    const double& z) const {
    (void)qOut;
    (void)qSide;
    (void)qLong;
    (void)x;
    (void)y;
    (void)z;

    // Placeholder for future short-range + Coulomb/plane-wave propagation.
    // For now, no propagated 3D wave function is available.
    return std::complex<double>(0.0, 0.0);
}

void CATS3D::SetShortRangePotential3D(double (*pot)(double*), CATSparameters& Pars) {
    ShortRangePotential3D = pot;

    if (PotPar3D) {
        delete PotPar3D;
        PotPar3D = nullptr;
    }

    if (pot) {
        PotPar3D = new CATSparameters(Pars);
    }

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS3D::RemoveShortRangePotential3D() {
    ShortRangePotential3D = nullptr;

    if (PotPar3D) {
        delete PotPar3D;
        PotPar3D = nullptr;
    }

    DeleteWaveFunctionCache();
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

bool CATS3D::HasShortRangePotential3D() const {
    return ShortRangePotential3D != nullptr;
}

// At the moment the propogating is follow the 1-D!
// To be impplemented to 3-D with 3-D (anisotropic) potential!
void CATS3D::PropagatingFunction3D(double& Basic,
                                   double& Full,
                                   const double& Radius,
                                   const double& Momentum,
                                   const unsigned short& l) const {
    const double rNu = Radius * FmToNu;

    // Basic = Coulomb + centrifugal - q^2
    Basic =
        2.0 * RedMass * CoulombPotential(rNu) +
        double(l) * (double(l) + 1.0) / (rNu * rNu + 1e-64) -
        Momentum * Momentum;

    Full = Basic;

    if (!ShortRangePotential3D || !PotPar3D) {
        return;
    }

    // Convention for future 3D short-range potential parameters:
    // variable 0 -> Radius [fm]
    // variable 1 -> Momentum [MeV]
    PotPar3D->SetVariable(0, Radius, true);
    PotPar3D->SetVariable(1, Momentum, true);

    double* Parameters = PotPar3D->GetParameters();
    Full = Basic + 2.0 * RedMass * ShortRangePotential3D(Parameters);
}


void CATS3D::ComputeWaveFunction() {
    if (!PotentialStatus()) {
        ComputedWaveFunction = false;
        return;
    }

    switch (CurrentWfMode) {
        case wfPlaneWave:
            ComputedWaveFunction = true;
            return;

        case wfCoulomb:
            // Workflow-compatible placeholder.
            // A full Coulomb wave function can be implemented later.
            ComputedWaveFunction = true;
            return;

        case wfExternal3D:
        case wfSchrodinger3D:
            // placeholder for future propagated 3D wave function
            ComputedWaveFunction = HasShortRangePotential3D();
            return;
        case wfNone:
        default:
            ComputedWaveFunction = false;
            return;
    }
}

void CATS3D::ComputeTotWaveFunction(const bool& ReallocateTotWaveFun) {
    (void)ReallocateTotWaveFun;

    if (!ComputedWaveFunction) {
        ComputeWaveFunction();
    }

    // First minimal CATS3D version:
    // no dedicated total-WF cache is stored yet.
    // |psi|^2 is evaluated on the fly in EvaluateWaveFunction2(...)
}

// Kill the cat !
void CATS3D::KillTheCat() {
    const int profileLevel = CATS3DProfileLevel();
    const auto tKill0 = CATS3DNow();
    if (!SourceStatus()) {
        throw std::runtime_error("CATS3D::KillTheCat: analytic source is not configured.");
    }

    if (!PotentialStatus()) {
        throw std::runtime_error("CATS3D::KillTheCat: wave-function mode is not ready.");
    }

    if (UseFastKernel) {
        if (profileLevel) {
            std::cout << "CATS3D::KillTheCat profile: fastKernel=1"
                      << " fastSubdiv=" << FastKernelSubdivisions
                      << " fastReady=" << FastKernelReady
                      << " wfComputed=" << ComputedWaveFunction
                      << " sourceGridReady=" << SourceGridReady
                      << " sourceUpdated=" << SourceUpdated << "\n";
        }
        if (!ComputedWaveFunction) {
            const auto tWf0 = CATS3DNow();
            ComputeWaveFunction();
            if (profileLevel) {
                std::cout << "  ComputeWaveFunction: "
                          << CATS3DSecondsSince(tWf0) << " s\n";
            }
        }
        if (!FastKernelReady) {
            const auto tFast0 = CATS3DNow();
            BuildFastWaveFunctionKernel();
            if (profileLevel) {
                std::cout << "  BuildFastWaveFunctionKernel total: "
                          << CATS3DSecondsSince(tFast0) << " s\n";
            }
        }
        const auto tFold0 = CATS3DNow();
        FoldFastSourceAndWF();
        if (profileLevel) {
            std::cout << "  FoldFastSourceAndWF total: "
                      << CATS3DSecondsSince(tFold0) << " s\n"
                      << "CATS3D::KillTheCat profile total: "
                      << CATS3DSecondsSince(tKill0) << " s\n";
        }
        return;
    }

    if (!SourceGridReady) {
        const auto tSetup0 = CATS3DNow();
        SetUpSourceGrid();
        if (profileLevel) {
            std::cout << "  SetUpSourceGrid: "
                      << CATS3DSecondsSince(tSetup0) << " s\n";
        }
    }

    if (!SourceUpdated) {
        const auto tUpdate0 = CATS3DNow();
        UpdateSourceGrid();
        if (profileLevel) {
            std::cout << "  UpdateSourceGrid: "
                      << CATS3DSecondsSince(tUpdate0) << " s\n";
        }
    }

    if (!ComputedWaveFunction) {
        const auto tWf0 = CATS3DNow();
        ComputeWaveFunction();
        if (profileLevel) {
            std::cout << "  ComputeWaveFunction: "
                      << CATS3DSecondsSince(tWf0) << " s\n";
        }
    }

    const auto tTotWf0 = CATS3DNow();
    ComputeTotWaveFunction(false);
    if (profileLevel) {
        std::cout << "  ComputeTotWaveFunction: "
                  << CATS3DSecondsSince(tTotWf0) << " s\n";
    }
    const auto tFold0 = CATS3DNow();
    FoldSourceAndWF();
    if (profileLevel) {
        std::cout << "  FoldSourceAndWF: "
                  << CATS3DSecondsSince(tFold0) << " s\n"
                  << "CATS3D::KillTheCat profile total: "
                  << CATS3DSecondsSince(tKill0) << " s\n";
    }
}

// source grid setup using CATSelder
void CATS3D::SetUpSourceGrid() {
    if (!SourceStatus()) {
        throw std::runtime_error("CATS3D::SetUpSourceGrid: analytic source is not configured.");
    }

    if (!GridBinX || !GridBinY || !GridBinZ) {
        throw std::runtime_error("CATS3D::SetUpSourceGrid: source grid bins are not configured.");
    }

    if (!MomBinCenter_out || !MomBinCenter_side || !MomBinCenter_long) {
        throw std::runtime_error("CATS3D::SetUpSourceGrid: momentum bins are not configured.");
    }

    // -----------------------------
    // Clean previous source-grid objects
    // -----------------------------
    if (kSourceGrid3D) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            if (!kSourceGrid3D[uOut]) continue;
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                if (!kSourceGrid3D[uOut][uSide]) continue;
                for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                    if (kSourceGrid3D[uOut][uSide][uLong]) {
                        delete kSourceGrid3D[uOut][uSide][uLong];
                        kSourceGrid3D[uOut][uSide][uLong] = nullptr;
                    }
                }
                delete[] kSourceGrid3D[uOut][uSide];
                kSourceGrid3D[uOut][uSide] = nullptr;
            }
            delete[] kSourceGrid3D[uOut];
            kSourceGrid3D[uOut] = nullptr;
        }
        delete[] kSourceGrid3D;
        kSourceGrid3D = nullptr;
    }

    if (BaseSourceGrid3D) {
        delete BaseSourceGrid3D;
        BaseSourceGrid3D = nullptr;
    }

    // -----------------------------
    // Settings for CATSelder
    // -----------------------------
    const short DIM = 3;

    // Use the already-defined OSL/source coordinate ranges
    double* MEAN = new double[DIM];
    double* LENGTH = new double[DIM];

    const double minX = GridBinX[0];
    const double maxX = GridBinX[NumGridBinsX];
    const double minY = GridBinY[0];
    const double maxY = GridBinY[NumGridBinsY];
    const double minZ = GridBinZ[0];
    const double maxZ = GridBinZ[NumGridBinsZ];

    MEAN[0]   = 0.5 * (minX + maxX);
    MEAN[1]   = 0.5 * (minY + maxY);
    MEAN[2]   = 0.5 * (minZ + maxZ);

    LENGTH[0] = (maxX - minX);
    LENGTH[1] = (maxY - minY);
    LENGTH[2] = (maxZ - minZ);

    // Same philosophy as original CATS:
    // if user has not supplied custom values, use safe defaults.
    // double LIMIT = 1.0 / 1024.0;
    // short MAXDEPTH = 14;
    // const short MINDEPTH = 7;
    // less end nodes to decrease the running time !!
    double LIMIT = 1.0 / 128.0;
    short MAXDEPTH = 9;
    const short MINDEPTH = 5;

    // Set the source variables to the average momentum point
    // Convention:
    //   Pars[0] = qOut
    //   Pars[1] = qSide
    //   Pars[2] = qLong
    //   Pars[3] = rOut
    //   Pars[4] = rSide
    //   Pars[5] = rLong
    // TODO: need to think about it!!!!!!
    const double qOutMean =
        0.5 * (MomBin_out[0] + MomBin_out[NumMomBins_out]);
    const double qSideMean =
        0.5 * (MomBin_side[0] + MomBin_side[NumMomBins_side]);
    const double qLongMean =
        0.5 * (MomBin_long[0] + MomBin_long[NumMomBins_long]);

    AnaSourcePar->SetVariable(0, qOutMean, false);
    AnaSourcePar->SetVariable(1, qSideMean, false);
    AnaSourcePar->SetVariable(2, qLongMean, false);

    // Base source grid: adaptive grid in (rOut, rSide, rLong)
    BaseSourceGrid3D = new CATSelder(
        DIM,
        MINDEPTH,
        MAXDEPTH,
        LIMIT,
        MEAN,
        LENGTH,
        this,
        AnaSourcePar,
        nullptr,
        0,
        true,
        true
    );

    // -----------------------------
    // Allocate kSourceGrid3D[uOut][uSide][uLong]
    // -----------------------------
    kSourceGrid3D = new CATSelder***[NumMomBins_out];
    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        kSourceGrid3D[uOut] = new CATSelder**[NumMomBins_side];
        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            kSourceGrid3D[uOut][uSide] = new CATSelder*[NumMomBins_long];
            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                kSourceGrid3D[uOut][uSide][uLong] = nullptr;
            }
        }
    }

    // Create one source-grid instance per 3D momentum bin
    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                AnaSourcePar->SetVariable(0, MomBinCenter_out[uOut], false);
                AnaSourcePar->SetVariable(1, MomBinCenter_side[uSide], false);
                AnaSourcePar->SetVariable(2, MomBinCenter_long[uLong], false);

                kSourceGrid3D[uOut][uSide][uLong] = new CATSelder(
                    BaseSourceGrid3D,
                    this,
                    AnaSourcePar,
                    nullptr,
                    0,
                    true,
                    true
                );
            }
        }
    }

    delete[] MEAN;
    delete[] LENGTH;

    SourceGridReady = true;
    SourceUpdated = true;
}

void CATS3D::UpdateSourceGrid() {
    if (!BaseSourceGrid3D) {
        SourceUpdated = false;
        return;
    }

    BaseSourceGrid3D->Update3D();

    if (kSourceGrid3D) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            if (!kSourceGrid3D[uOut]) continue;
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                if (!kSourceGrid3D[uOut][uSide]) continue;
                for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                    if (kSourceGrid3D[uOut][uSide][uLong]) {
                        kSourceGrid3D[uOut][uSide][uLong]->Update3D();
                    }
                }
            }
        }
    }

    SourceUpdated = true;
}

void CATS3D::BuildFastWaveFunctionKernel() {
    const int profileLevel = CATS3DProfileLevel();
    const auto tBuild0 = CATS3DNow();
    if (!GridBinCenterX || !GridBinCenterY || !GridBinCenterZ) {
        throw std::runtime_error("CATS3D::BuildFastWaveFunctionKernel: source grid bins are not configured.");
    }
    if (!MomBinCenter_out || !MomBinCenter_side || !MomBinCenter_long) {
        throw std::runtime_error("CATS3D::BuildFastWaveFunctionKernel: momentum bins are not configured.");
    }
    if (!PotentialStatus()) {
        throw std::runtime_error("CATS3D::BuildFastWaveFunctionKernel: wave-function mode is not ready.");
    }

    const unsigned subdivisions = FastKernelSubdivisions ? FastKernelSubdivisions : 1;
    const unsigned nMom = NumMomBins_out * NumMomBins_side * NumMomBins_long;
    const unsigned nBaseGrid = NumGridBinsX * NumGridBinsY * NumGridBinsZ;
    const size_t nSubCells = size_t(subdivisions) * size_t(subdivisions) * size_t(subdivisions);
    const size_t nGrid = size_t(nBaseGrid) * nSubCells;
    if (profileLevel) {
        std::cout << "CATS3D::BuildFastWaveFunctionKernel profile: nMom="
                  << nMom << " nBaseGrid=" << nBaseGrid
                  << " subdivisions=" << subdivisions
                  << " nKernelGrid=" << nGrid
                  << " totalEval=" << size_t(nMom) * size_t(nGrid) << "\n";
    }
    FastWaveFunction2.assign(size_t(nMom) * size_t(nGrid), 0.0f);

    auto fastGridIndex = [this, subdivisions](const unsigned uX,
                                              const unsigned uY,
                                              const unsigned uZ,
                                              const unsigned sX,
                                              const unsigned sY,
                                              const unsigned sZ) -> size_t {
        const size_t baseIndex = GetFlatGridIndex(uX, uY, uZ);
        if (subdivisions <= 1) return baseIndex;
        return (((baseIndex * subdivisions + sX) * subdivisions + sY) * subdivisions + sZ);
    };

    unsigned qCounter = 0;
    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                const auto tQ0 = CATS3DNow();
                const unsigned qIndex = GetFlatMomIndex(uOut, uSide, uLong);

                for (unsigned uX = 0; uX < NumGridBinsX; ++uX) {
                    for (unsigned uY = 0; uY < NumGridBinsY; ++uY) {
                        for (unsigned uZ = 0; uZ < NumGridBinsZ; ++uZ) {
                            for (unsigned sX = 0; sX < subdivisions; ++sX) {
                                const double rOut = subdivisions <= 1
                                    ? GridBinCenterX[uX]
                                    : GridBinX[uX] + (double(sX) + 0.5) * GridStepX / double(subdivisions);
                                for (unsigned sY = 0; sY < subdivisions; ++sY) {
                                    const double rSide = subdivisions <= 1
                                        ? GridBinCenterY[uY]
                                        : GridBinY[uY] + (double(sY) + 0.5) * GridStepY / double(subdivisions);
                                    for (unsigned sZ = 0; sZ < subdivisions; ++sZ) {
                                        const double rLong = subdivisions <= 1
                                            ? GridBinCenterZ[uZ]
                                            : GridBinZ[uZ] + (double(sZ) + 0.5) * GridStepZ / double(subdivisions);
                                        const size_t rIndex = fastGridIndex(uX, uY, uZ, sX, sY, sZ);
                                        const double wf2 =
                                            EvaluateWaveFunction2(uOut, uSide, uLong,
                                                                  rOut, rSide, rLong);
                                        FastWaveFunction2[size_t(qIndex) * nGrid + rIndex] =
                                            float(wf2);
                                    }
                                }
                            }
                        }
                    }
                }
                ++qCounter;
                if (profileLevel > 1) {
                    std::cout << "  fast WF q-bin " << qCounter << "/"
                              << nMom << " time="
                              << CATS3DSecondsSince(tQ0) << " s\n";
                }
            }
        }
    }

    FastKernelReady = true;
    ComputedWaveFunction = true;
    if (profileLevel) {
        std::cout << "CATS3D::BuildFastWaveFunctionKernel profile total: "
                  << CATS3DSecondsSince(tBuild0) << " s\n";
    }
}

void CATS3D::FoldFastSourceAndWF() {
    const int profileLevel = CATS3DProfileLevel();
    const auto tFoldFast0 = CATS3DNow();
    const unsigned subdivisions = FastKernelSubdivisions ? FastKernelSubdivisions : 1;
    if (!FastKernelReady || FastWaveFunction2.empty()) {
        BuildFastWaveFunctionKernel();
    }
    if (!SourceStatus()) {
        throw std::runtime_error("CATS3D::FoldFastSourceAndWF: analytic source is not configured.");
    }

    AllocateCorrFun();

    const unsigned nBaseGrid = NumGridBinsX * NumGridBinsY * NumGridBinsZ;
    const size_t nSubCells = size_t(subdivisions) * size_t(subdivisions) * size_t(subdivisions);
    const size_t nGrid = size_t(nBaseGrid) * nSubCells;
    std::vector<double> sourceWeight(nGrid, 0.0);
    auto fastGridIndex = [this, subdivisions](const unsigned uX,
                                              const unsigned uY,
                                              const unsigned uZ,
                                              const unsigned sX,
                                              const unsigned sY,
                                              const unsigned sZ) -> size_t {
        const size_t baseIndex = GetFlatGridIndex(uX, uY, uZ);
        if (subdivisions <= 1) return baseIndex;
        return (((baseIndex * subdivisions + sX) * subdivisions + sY) * subdivisions + sZ);
    };
    if (profileLevel) {
        const unsigned nMom = NumMomBins_out * NumMomBins_side * NumMomBins_long;
        std::cout << "CATS3D::FoldFastSourceAndWF profile: nMom="
                  << nMom << " nBaseGrid=" << nBaseGrid
                  << " subdivisions=" << subdivisions
                  << " nKernelGrid=" << nGrid
                  << " grid=(" << NumGridBinsX << ","
                  << NumGridBinsY << "," << NumGridBinsZ << ")\n";
    }

    double sourceInt = 0.0;
    const double cellVolume = GridStepX * GridStepY * GridStepZ;
    const double subCellVolume = cellVolume /
        double(subdivisions * subdivisions * subdivisions);

    // Fast kernel mode assumes the analytic source depends on r and source
    // parameters, not explicitly on q. This is true for the current Gaussian
    // and DLM_CleverReso3D sources used by the 3-D fit macro.
    const double qOutRef = 0.0;
    const double qSideRef = 0.0;
    const double qLongRef = 0.0;

    const auto tSource0 = CATS3DNow();
    if (subdivisions <= 1) {
        for (unsigned uX = 0; uX < NumGridBinsX; ++uX) {
            const auto tPlane0 = CATS3DNow();
            const double rOut = GridBinCenterX[uX];
            for (unsigned uY = 0; uY < NumGridBinsY; ++uY) {
                const double rSide = GridBinCenterY[uY];
                for (unsigned uZ = 0; uZ < NumGridBinsZ; ++uZ) {
                    const double rLong = GridBinCenterZ[uZ];
                    const unsigned rIndex = GetFlatGridIndex(uX, uY, uZ);
                    const double src =
                        EvaluateTheSource(qOutRef, qSideRef, qLongRef,
                                          rOut, rSide, rLong);
                    const double weight = src * cellVolume;
                    sourceWeight[rIndex] += weight;
                    sourceInt += weight;
                }
            }
            if (profileLevel > 1) {
                std::cout << "  source x-plane " << (uX+1) << "/"
                          << NumGridBinsX << " time="
                          << CATS3DSecondsSince(tPlane0)
                          << " s cumulative="
                          << CATS3DSecondsSince(tSource0) << " s\n";
            }
        }
    }
    else {
        for (unsigned uX = 0; uX < NumGridBinsX; ++uX) {
            const auto tPlane0 = CATS3DNow();
            for (unsigned uY = 0; uY < NumGridBinsY; ++uY) {
                for (unsigned uZ = 0; uZ < NumGridBinsZ; ++uZ) {
                    for (unsigned sX = 0; sX < subdivisions; ++sX) {
                        const double rOut =
                            GridBinX[uX] + (double(sX) + 0.5) * GridStepX / double(subdivisions);
                        for (unsigned sY = 0; sY < subdivisions; ++sY) {
                            const double rSide =
                                GridBinY[uY] + (double(sY) + 0.5) * GridStepY / double(subdivisions);
                            for (unsigned sZ = 0; sZ < subdivisions; ++sZ) {
                                const double rLong =
                                    GridBinZ[uZ] + (double(sZ) + 0.5) * GridStepZ / double(subdivisions);
                                const double src =
                                    EvaluateTheSource(qOutRef, qSideRef, qLongRef,
                                                      rOut, rSide, rLong);
                                const double weight = src * subCellVolume;
                                const size_t rIndex = fastGridIndex(uX, uY, uZ, sX, sY, sZ);
                                sourceWeight[rIndex] += weight;
                                sourceInt += weight;
                            }
                        }
                    }
                }
            }
            if (profileLevel > 1) {
                std::cout << "  source-subcell x-plane " << (uX+1) << "/"
                          << NumGridBinsX << " time="
                          << CATS3DSecondsSince(tPlane0)
                          << " s cumulative="
                          << CATS3DSecondsSince(tSource0) << " s\n";
            }
        }
    }
    const double sourceSec = CATS3DSecondsSince(tSource0);
    if (profileLevel) {
        std::cout << "  sourceWeight loop: " << sourceSec
                  << " s sourceInt=" << sourceInt << "\n";
    }

    if (sourceInt <= 0.0) {
        for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
            for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
                for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                    kCorrFun[uOut][uSide][uLong] = 0.0;
                    kCorrFunErr[uOut][uSide][uLong] = 1e6;
                }
            }
        }
        ComputedCorrFunction = true;
        SourceUpdated = true;
        if (profileLevel) {
            std::cout << "CATS3D::FoldFastSourceAndWF profile total: "
                      << CATS3DSecondsSince(tFoldFast0)
                      << " s (empty source)\n";
        }
        return;
    }

    const auto tCorr0 = CATS3DNow();
    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {
                double corr = 0.0;
                const unsigned qIndex = GetFlatMomIndex(uOut, uSide, uLong);
                const size_t offset = size_t(qIndex) * nGrid;
                for (size_t rIndex = 0; rIndex < nGrid; ++rIndex) {
                    corr += sourceWeight[rIndex] *
                            double(FastWaveFunction2[offset + rIndex]);
                }

                kCorrFun[uOut][uSide][uLong] = corr / sourceInt;
                kCorrFunErr[uOut][uSide][uLong] = 0.0;
            }
        }
    }
    if (profileLevel) {
        std::cout << "  correlation fold loop: "
                  << CATS3DSecondsSince(tCorr0) << " s\n"
                  << "CATS3D::FoldFastSourceAndWF profile total: "
                  << CATS3DSecondsSince(tFoldFast0) << " s\n";
    }

    SourceGridReady = true;
    SourceUpdated = true;
    ComputedCorrFunction = true;
}

void CATS3D::FoldSourceAndWF() {
    if (!kSourceGrid3D) {
        throw std::runtime_error("CATS3D::FoldSourceAndWF: source grid is not set up.");
    }

    AllocateCorrFun();

    if (!SourceUpdated) {
        UpdateSourceGrid();
    }

    for (unsigned uOut = 0; uOut < NumMomBins_out; ++uOut) {
        for (unsigned uSide = 0; uSide < NumMomBins_side; ++uSide) {
            for (unsigned uLong = 0; uLong < NumMomBins_long; ++uLong) {

                CATSelder* grid = kSourceGrid3D[uOut][uSide][uLong];
                if (!grid) {
                    kCorrFunErr[uOut][uSide][uLong] = 1e6;
                    continue;
                }

                const unsigned NumGridPts = grid->GetNumEndNodes();

                double corr = 0.0;
                double corrErr2 = 0.0;
                double sourceInt = 0.0;

                for (unsigned uGrid = 0; uGrid < NumGridPts; ++uGrid) {
                    const double rOut  = grid->GetParValue(uGrid, 0);
                    const double rSide = grid->GetParValue(uGrid, 1);
                    const double rLong = grid->GetParValue(uGrid, 2);

                    const double SourceVal = grid->GetGridValue(uGrid);
                    const double SourceErr = grid->GetGridError(uGrid);

                    if (!SourceVal) continue;

                    sourceInt += SourceVal;

                    const double WaveFunVal =
                        EvaluateWaveFunction2(uOut, uSide, uLong, rOut, rSide, rLong);

                    const double Integrand = SourceVal * WaveFunVal;
                    corr += Integrand;
                    corrErr2 += std::pow(Integrand * SourceErr, 2.0);
                }

                // For a properly normalized source grid, sourceInt should be ~1.
                // If not exactly 1 because of grid/adaptive precision, renormalize.
                if (sourceInt > 0.0) {
                    corr /= sourceInt;
                    corrErr2 /= (sourceInt * sourceInt);
                }

                kCorrFun[uOut][uSide][uLong] = corr;
                kCorrFunErr[uOut][uSide][uLong] = std::sqrt(corrErr2);
            }
        }
    }

    ComputedCorrFunction = true;
}
