#include "CATS3D.h"
#include "CATSconstants.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<double> qEdgesMeV() {
    std::vector<double> e{0.0};
    for (int q = 5; q <= 40; q += 5) e.push_back(q);
    for (int q = 50; q <= 100; q += 10) e.push_back(q);
    for (int q = 120; q <= 200; q += 20) e.push_back(q);
    return e;
}

// Resolve the rapidly varying femtoscopic peak first, then coarsen toward
// the high-q baseline. The same axis is used for out, side, and long.
std::vector<double> farmQEdgesMeV() {
    std::vector<double> e{0.0};
    for (int i = 1; i <= 12; ++i) e.push_back(2.5 * i);       // 0..30
    for (int i = 1; i <= 10; ++i) e.push_back(30.0 + 5.0 * i); // 30..80
    for (int i = 1; i <= 4; ++i) e.push_back(80.0 + 10.0 * i); // 80..120
    for (int i = 1; i <= 4; ++i) e.push_back(120.0 + 20.0 * i);// 120..200
    return e;
}

std::vector<double> rEdgesFm() {
    std::vector<double> positive{0.0};
    for (int i = 1; i <= 10; ++i) positive.push_back(1.5 * i);
    for (int i = 1; i <= 5; ++i) positive.push_back(15.0 + 3.0 * i);
    positive.push_back(35.0);
    positive.push_back(40.0);
    std::vector<double> both;
    both.reserve(2 * positive.size() - 1);
    for (auto it = positive.rbegin(); it != positive.rend(); ++it) {
        if (*it > 0.0) both.push_back(-*it);
    }
    both.insert(both.end(), positive.begin(), positive.end());
    return both;
}


// Finer near the Coulomb-sensitive origin; retain the full +/-40 fm range.
std::vector<double> farmREdgesFm() {
    std::vector<double> positive{0.0};
    for (int i = 1; i <= 16; ++i) positive.push_back(0.75 * i);      // 0..12
    for (int i = 1; i <= 6; ++i) positive.push_back(12.0 + 2.0 * i); // 12..24
    for (int i = 1; i <= 5; ++i) positive.push_back(24.0 + 3.2 * i); // 24..40
    std::vector<double> both;
    both.reserve(2 * positive.size() - 1);
    for (auto it = positive.rbegin(); it != positive.rend(); ++it) {
        if (*it > 0.0) both.push_back(-*it);
    }
    both.insert(both.end(), positive.begin(), positive.end());
    return both;
}

std::vector<double> centers(const std::vector<double>& edges) {
    std::vector<double> result;
    result.reserve(edges.size() - 1);
    for (size_t i = 0; i + 1 < edges.size(); ++i) result.push_back(0.5 * (edges[i] + edges[i + 1]));
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: build_adaptive_pipi_coulomb_kernel OUTPUT.root [--farm] [--smoke] [--estimate]\n";
        return 2;
    }
    bool smoke = false, farm = false, estimate = false;
    for (int i = 2; i < argc; ++i) {
        const std::string flag(argv[i]);
        if (flag == "--smoke") smoke = true;
        else if (flag == "--farm") farm = true;
        else if (flag == "--estimate") estimate = true;
        else {
            std::cerr << "Unknown option: " << flag << '\n';
            return 2;
        }
    }

    auto qEdges = smoke ? std::vector<double>{0.0,5.0}
                        : (farm ? farmQEdgesMeV() : qEdgesMeV());
    auto rEdges = farm ? farmREdgesFm() : rEdgesFm();
    auto qCenters = centers(qEdges);
    auto rCenters = centers(rEdges);
    const unsigned nQ = unsigned(qCenters.size());
    const unsigned nR = unsigned(rCenters.size());
    const size_t nR3 = size_t(nR) * nR * nR;
    const size_t total = size_t(nQ) * nQ * nQ * nR3;
    std::cout << "Adaptive pi+pi+ Coulomb + QS wavefunction kernel\n"
              << "q: " << nQ << "^3 positive-octant cells, edges 0.."
              << qEdges.back() << " MeV/c\n"
              << "r: " << nR << "^3 Cartesian cells, edges -40..40 fm\n"
              << "raw float payload: " << total << " values, " << double(total * sizeof(float)) / 1e9 << " GB\n"
              << "mode: " << (farm ? "farm" : "local")
              << (smoke ? " smoke (one q cell)" : " full") << std::endl;
    if (estimate) return 0;

    CATS3D kitty;
    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5 * Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);

    TFile out(argv[1], "RECREATE");
    if (out.IsZombie()) throw std::runtime_error("Cannot create output ROOT file");
    out.SetCompressionLevel(1);

    unsigned formatVersion = 1;
    int pdg1 = 211, pdg2 = 211;
    int q1q2 = 1;
    double reducedMassMeV = 0.5 * Mass_Pic;
    TTree metadata("KernelMetadata", "Nonuniform Cartesian wavefunction kernel metadata");
    metadata.Branch("formatVersion", &formatVersion);
    metadata.Branch("pdg1", &pdg1);
    metadata.Branch("pdg2", &pdg2);
    metadata.Branch("q1q2", &q1q2);
    metadata.Branch("reducedMassMeV", &reducedMassMeV);
    metadata.Branch("qOutEdgesMeV", &qEdges);
    metadata.Branch("qSideEdgesMeV", &qEdges);
    metadata.Branch("qLongEdgesMeV", &qEdges);
    metadata.Branch("rOutEdgesFm", &rEdges);
    metadata.Branch("rSideEdgesFm", &rEdges);
    metadata.Branch("rLongEdgesFm", &rEdges);
    metadata.Fill();
    metadata.Write();

    unsigned iqOut = 0, iqSide = 0, iqLong = 0;
    double qOutMeV = 0, qSideMeV = 0, qLongMeV = 0;
    std::vector<float> wavefunction2(nR3);
    TTree kernel("KernelByQ", "One spatial |psi|^2 array per positive q cell");
    kernel.Branch("iqOut", &iqOut);
    kernel.Branch("iqSide", &iqSide);
    kernel.Branch("iqLong", &iqLong);
    kernel.Branch("qOutMeV", &qOutMeV);
    kernel.Branch("qSideMeV", &qSideMeV);
    kernel.Branch("qLongMeV", &qLongMeV);
    kernel.Branch("wavefunction2", &wavefunction2);
    kernel.SetAutoFlush(10);
    kernel.SetAutoSave(100'000'000);

    const auto started = std::chrono::steady_clock::now();
    for (iqOut = 0; iqOut < nQ; ++iqOut) {
        qOutMeV = qCenters[iqOut];
        kitty.SetMomBinsOut(1, qOutMeV, qOutMeV);
        for (iqSide = 0; iqSide < nQ; ++iqSide) {
            qSideMeV = qCenters[iqSide];
            kitty.SetMomBinsSide(1, qSideMeV, qSideMeV);
            for (iqLong = 0; iqLong < nQ; ++iqLong) {
                qLongMeV = qCenters[iqLong];
                kitty.SetMomBinsLong(1, qLongMeV, qLongMeV);
                size_t index = 0;
                for (double x : rCenters) {
                    for (double y : rCenters) {
                        for (double z : rCenters) {
                            const double value = kitty.EvaluateWaveFunction2(0, 0, 0, x, y, z);
                            if (!std::isfinite(value) || value < 0.0) {
                                throw std::runtime_error("Nonfinite or negative wavefunction value");
                            }
                            wavefunction2[index++] = float(value);
                        }
                    }
                }
                kernel.Fill();
                if (smoke) break;
            }
            if (smoke) break;
        }
        kernel.AutoSave("SaveSelf");
        out.Flush();
        const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        std::cout << "completed qOut plane " << iqOut + 1 << "/" << nQ
                  << " entries=" << kernel.GetEntries() << " elapsed=" << elapsed << " s" << std::endl;
        if (smoke) break;
    }
    out.cd();
    kernel.Write("", TObject::kOverwrite);
    out.Close();
    return 0;
}
