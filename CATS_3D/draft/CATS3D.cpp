#include "CATS3D.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

// constructor
CATS3D::CATS3D()
    : RedMass(0.0),
      QuantumStatistics(-1),
      IdenticalParticles(false),
      Q1Q2(0),
      CurrentWfMode(wfNone),

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
      AnaSourcePar(nullptr),
      NumSourcePars(0),

      kCorrFun(nullptr),
      kCorrFunErr(nullptr),

      SourceGridReady(false),
      SourceUpdated(false),
      ComputedWaveFunction(false),
      ComputedCorrFunction(false),
      BaseSourceGrid3D(nullptr),
      kSourceGrid3D(nullptr) {
    pdgID[0] = 0;
    pdgID[1] = 0;
}

// destructor
CATS3D::~CATS3D() {
    DeleteCorrFun();
    DeleteGridBins();
    DeleteAllMomBinDefinitions();

    if (AnaSourcePar) {
        delete AnaSourcePar;
        AnaSourcePar = nullptr;
    }

    if (BaseSourceGrid3D) {
        delete BaseSourceGrid3D;
        BaseSourceGrid3D = nullptr;
    }

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

void CATS3D::DeleteWaveFunctionCache() {
    // First minimal CATS3D version:
    // no dedicated wave-function cache is stored yet.
    // This hook is kept for future extensions.
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
        if (uBin < NumMomBins_out) {
            MomBinCenter_out[uBin] = 0.5 * (MomBin_out[uBin] + MomBin_out[uBin + 1]);
        }
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
        if (uBin < NumMomBins_side) {
            MomBinCenter_side[uBin] = 0.5 * (MomBin_side[uBin] + MomBin_side[uBin + 1]);
        }
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
        if (uBin < NumMomBins_long) {
            MomBinCenter_long[uBin] = 0.5 * (MomBin_long[uBin] + MomBin_long[uBin + 1]);
        }
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
    if (uBinOut >= NumMomBins_out || uBinSide >= NumMomBins_side || uBinLong >= NumMomBins_long) {
        return 0.0;
    }

    const double qOut = MomBinCenter_out ? MomBinCenter_out[uBinOut] : 0.0;
    const double qSide = MomBinCenter_side ? MomBinCenter_side[uBinSide] : 0.0;
    const double qLong = MomBinCenter_long ? MomBinCenter_long[uBinLong] : 0.0;

    return std::sqrt(qOut * qOut + qSide * qSide + qLong * qLong);
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
        if (uBin < NumGridBinsX) {
            GridBinCenterX[uBin] = 0.5 * (GridBinX[uBin] + GridBinX[uBin + 1]);
        }
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
        if (uBin < NumGridBinsY) {
            GridBinCenterY[uBin] = 0.5 * (GridBinY[uBin] + GridBinY[uBin + 1]);
        }
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
        if (uBin < NumGridBinsZ) {
            GridBinCenterZ[uBin] = 0.5 * (GridBinZ[uBin] + GridBinZ[uBin + 1]);
        }
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
    if(!AnalyticSource || !Pars){
        printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported a crash! This is most likely a bug!\n");
        printf("         Please contact the developers, reporting that Pars=NULL\n");    
        return 0;
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
    if (!AnalyticSource || !AnaSourcePar) {
        return 0.0;
    }

    AnaSourcePar->SetVariable(0, qOut, false);
    AnaSourcePar->SetVariable(1, qSide, false);
    AnaSourcePar->SetVariable(2, qLong, false);
    AnaSourcePar->SetVariable(3, rOut, false);
    AnaSourcePar->SetVariable(4, rSide, false);
    AnaSourcePar->SetVariable(5, rLong, false);

    double* Pars = AnaSourcePar->GetParameters();
    return AnalyticSource(Pars);
}

// status query
bool CATS3D::SourceStatus() const {
    return AnalyticSource != nullptr && AnaSourcePar != nullptr;
}

bool CATS3D::PotentialStatus() const {
    // In the first minimal CATS3D version:
    // - plane wave needs no extra potential setup
    // - Coulomb only needs Q1Q2 and RedMass at evaluation time
    // - external3D / Schrodinger3D are not implemented yet
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
            return false;
        default:
            return false;
    }
}

bool CATS3D::CkStatus() const {
    return ComputedCorrFunction && kCorrFun != nullptr;
}

double CATS3D::CoulombPotential(const double& Radius) const {
    // Define of AlphaFS is in CATSconstants
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
double CATS3D::EvaluateWaveFunction2(const unsigned& uBinOut,
                                     const unsigned& uBinSide,
                                     const unsigned& uBinLong,
                                     const double& x,
                                     const double& y,
                                     const double& z) const {
    if (uBinOut >= NumMomBins_out) return 0.0;
    if (uBinSide >= NumMomBins_side) return 0.0;
    if (uBinLong >= NumMomBins_long) return 0.0;

    const double qOut  = GetMomentumOut(uBinOut);
    const double qSide = GetMomentumSide(uBinSide);
    const double qLong = GetMomentumLong(uBinLong);

    // q·r in natural units:
    // q [MeV], r [fm], so multiply r by FmToNu [1/MeV/fm]
    const double qDotR =
        qOut  * x * FmToNu +
        qSide * y * FmToNu +
        qLong * z * FmToNu;

    switch (CurrentWfMode) {
        case wfPlaneWave:
            if (IdenticalParticles) {
                return 1.0 + std::cos(2.0 * qDotR);
            } else {
                return 1.0;
            }

        case wfCoulomb:
            // First minimal version:
            // keep the workflow slot, but do not yet claim a full Coulomb solution.
            // For now we fall back to the plane-wave expression.
            if (IdenticalParticles) {
                return 1.0 + std::cos(2.0 * qDotR);
            } else {
                return 1.0;
            }

        case wfExternal3D:
        case wfSchrodinger3D:
        case wfNone:
        default:
            return 0.0;
    }
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
    if (!SourceStatus()) {
        throw std::runtime_error("CATS3D::KillTheCat: analytic source is not configured.");
    }

    if (!PotentialStatus()) {
        throw std::runtime_error("CATS3D::KillTheCat: wave-function mode is not ready.");
    }

    if (!SourceGridReady) {
        SetUpSourceGrid();
    }

    if (!SourceUpdated) {
        UpdateSourceGrid();
    }

    if (!ComputedWaveFunction) {
        ComputeWaveFunction();
    }

    ComputeTotWaveFunction(false);
    FoldSourceAndWF();
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
    double LIMIT = 1.0 / 1024.0;
    short MAXDEPTH = 14;
    const short MINDEPTH = 7;

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

void CATS3D::FoldSourceAndWF() {
    if (!kSourceGrid3D) {
        throw std::runtime_error("CATS3D::FoldSourceAndWF: source grid is not set up.");
    }

    // Reallocate correlation arrays every time for now.
    // This is simpler and safer during development.
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
