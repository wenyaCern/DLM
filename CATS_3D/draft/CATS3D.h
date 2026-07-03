#ifndef CATS3D_H
#define CATS3D_H

#include <complex>
#include "CATStools.h"
#include "CATSconstants.h"
#include "CATS.h"

class CATS3D;
class CATSelder;
class CATSnode;
class CATSparameters;

class CATS3D {
friend class CATSelder;
public:
    CATS3D();
    ~CATS3D();

    // ---------------------------
    // basic physics setup
    // ---------------------------
    void SetRedMass(const double& redMass);
    double GetRedMass() const;

    void SetPdgId(const int& id1, const int& id2);
    void GetPdgId(int& id1, int& id2) const;

    // 0 = off; 1 = on; else auto from PDG IDs
    void SetQuantumStatistics(short qs);
    short GetQuantumStatistics() const;

    void SetQ1Q2(const int& q1q2);
    int GetQ1Q2() const;

    // ---------------------------
    // wave-function backend mode
    // ---------------------------
    enum WfMode {
        wfNone = 0,
        wfPlaneWave = 1,
        wfCoulomb = 2,
        wfExternal3D = 3,
        wfSchrodinger3D = 4
    };

    void SetWfMode(const WfMode& mode);
    WfMode GetWfMode() const;
    bool GetIdenticalParticles() const;

    // ---------------------------
    // momentum binning in 3D
    // ---------------------------
    void SetMomBinsOut(const unsigned& nBins, const double& minMom, const double& maxMom);
    void SetMomBinsSide(const unsigned& nBins, const double& minMom, const double& maxMom);
    void SetMomBinsLong(const unsigned& nBins, const double& minMom, const double& maxMom);

    unsigned GetNumMomBinsOut() const;
    unsigned GetNumMomBinsSide() const;
    unsigned GetNumMomBinsLong() const;

    double GetMomentumOut(const unsigned& uBin) const;
    double GetMomentumSide(const unsigned& uBin) const;
    double GetMomentumLong(const unsigned& uBin) const;

    double GetMomBinLowEdgeOut(const unsigned& uBin) const;
    double GetMomBinUpEdgeOut(const unsigned& uBin) const;

    double GetMomBinLowEdgeSide(const unsigned& uBin) const;
    double GetMomBinUpEdgeSide(const unsigned& uBin) const;

    double GetMomBinLowEdgeLong(const unsigned& uBin) const;
    double GetMomBinUpEdgeLong(const unsigned& uBin) const;

    // total |k| from 3D bin centers
    double GetTotMomentum3D(const unsigned& uBinOut,
                            const unsigned& uBinSide,
                            const unsigned& uBinLong) const;

    // ---------------------------
    // source grid in Cartesian 3D
    // ---------------------------
    void SetGridX(const unsigned& nBins, const double& minX, const double& maxX);
    void SetGridY(const unsigned& nBins, const double& minY, const double& maxY);
    void SetGridZ(const unsigned& nBins, const double& minZ, const double& maxZ);

    unsigned GetNumGridBinsX() const;
    unsigned GetNumGridBinsY() const;
    unsigned GetNumGridBinsZ() const;

    double GetGridCenterX(const unsigned& uBin) const;
    double GetGridCenterY(const unsigned& uBin) const;
    double GetGridCenterZ(const unsigned& uBin) const;

    double GetGridStepX() const;
    double GetGridStepY() const;
    double GetGridStepZ() const;

    // ---------------------------
    // analytic source interface
    // ---------------------------
    void SetAnaSource(double (*AS)(double*), CATSparameters& Pars);
    void SetAnaSource(const unsigned& WhichPar, const double& Value, const bool& SmallChange=false);
    double GetAnaSourcePar(const unsigned& WhichPar) const;

    // ---------------------------
    // main workflow
    // ---------------------------
    void KillTheCat();

    bool SourceStatus() const;
    bool PotentialStatus() const;
    bool CkStatus() const;

    void ComputeWaveFunction();
    void ComputeTotWaveFunction(const bool& ReallocateTotWaveFun = true);
    void FoldSourceAndWF();

    // ---------------------------
    // output
    // ---------------------------
    double GetCorrFun(const unsigned& uBinOut,
                      const unsigned& uBinSide,
                      const unsigned& uBinLong) const;

    double GetCorrFunErr(const unsigned& uBinOut,
                         const unsigned& uBinSide,
                         const unsigned& uBinLong) const;

protected:
    // ---------------------------
    // internal helpers
    // ---------------------------
    void SetUpSourceGrid();
    void UpdateSourceGrid();

    double CoulombPotential(const double& Radius) const;

    double EvaluateTheSource(CATSparameters* Pars) const;
    double EvaluateTheSource(const double& qOut,
                             const double& qSide,
                             const double& qLong,
                             const double& rOut,
                             const double& rSide,
                             const double& rLong) const;

    double EvaluateWaveFunction2(const unsigned& uBinOut,
                                 const unsigned& uBinSide,
                                 const unsigned& uBinLong,
                                 const double& x,
                                 const double& y,
                                 const double& z) const;

    // utility allocation helpers
    void DeleteAllMomBinDefinitions();
    void DeleteGridBins();
    void DeleteCorrFun();
    void DeleteWaveFunctionCache();

    CATSelder* BaseSourceGrid3D;
    CATSelder**** kSourceGrid3D;

private:
    // ---------------------------
    // basic physical parameters
    // ---------------------------
    double RedMass;
    int pdgID[2];
    short QuantumStatistics;
    bool IdenticalParticles;
    int Q1Q2;

    WfMode CurrentWfMode;

    // ---------------------------
    // momentum bins
    // ---------------------------
    unsigned NumMomBins_out;
    unsigned NumMomBins_side;
    unsigned NumMomBins_long;

    double* MomBin_out;
    double* MomBin_side;
    double* MomBin_long;

    double* MomBinCenter_out;
    double* MomBinCenter_side;
    double* MomBinCenter_long;

    // ---------------------------
    // source grid bins in x,y,z
    // ---------------------------
    unsigned NumGridBinsX;
    unsigned NumGridBinsY;
    unsigned NumGridBinsZ;

    double* GridBinX;
    double* GridBinY;
    double* GridBinZ;

    double* GridBinCenterX;
    double* GridBinCenterY;
    double* GridBinCenterZ;

    double GridStepX;
    double GridStepY;
    double GridStepZ;

    // ---------------------------
    // analytic source definition
    // ---------------------------
    double (*AnalyticSource)(double*);
    CATSparameters* AnaSourcePar;
    unsigned NumSourcePars;

    // ---------------------------
    // final correlation function
    // ---------------------------
    double*** kCorrFun;
    double*** kCorrFunErr;

    // ---------------------------
    // status flags
    // ---------------------------
    bool SourceGridReady;
    bool SourceUpdated;
    bool ComputedWaveFunction;
    bool ComputedCorrFunction;
};

#endif // CATS3D_H
