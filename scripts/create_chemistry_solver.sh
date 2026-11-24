#!/bin/bash
# RFAM Chemistry Integration - Step 1: Create modified solver with 1-step chemistry

cat > /d/realFluidAllMach/solver_3d_tetra_chemistry.cpp << 'EOF'
// RFAM tetrahedral solver WITH CHEMISTRY
// Based on apps/solver_3d_tetra.cpp + GlobalReaction1Step integration

// [Copy all includes from original apps/solver_3d_tetra.cpp lines 1-50]
#include "TetraMesh.H"
#include "TetraFiniteVolume.H"
#include "TetraLUSGS.H"
#include "CompleteCellState.H"
#include "EquilibriumThermoModel.H"
#include "CoolPropThermoModel.H"
#include "TableAcoustic.H"
#include "ThermoCounters.H"
#include "StateAcceptance.H"
#include "TetraStateRecovery.H"
#include "TetraPressureOutlet.H"
#include "TetraInletAllocation.H"
#include "TetraBlueprintWall.H"
#include "TetraMolecularAdmission.H"
#include "TetraSteadyMonitor.H"

// NEW: Add chemistry headers
#include "src/combustion/FiniteRate/GlobalReaction1Step.H"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>
#include <iterator>
#include <map>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef RFAM_USE_OPENMP
#include <omp.h>
#endif

using namespace std;

// Chemistry activation temperature threshold
constexpr double ACTIVATION_T_MIN = 750.0;  // K, below this no reaction

// Simple 1-step chemistry integration
// This is a MINIMAL implementation for proof-of-concept
class SimpleChemistryIntegrator {
private:
    GlobalReaction1Step chemistry_;
    bool enabled_;
    double activationTemperature_;
    
public:
    SimpleChemistryIntegrator(bool enable = true, double Tmin = 750.0)
        : chemistry_(8.6e11, 0.0, 125520.0, 1.0, 1.5)
        , enabled_(enable)
        , activationTemperature_(Tmin)
    {}
    
    // Compute chemistry source terms for one cell
    // Returns energy source qdot [W/m³] and species sources omega [kg/m³/s]
    void computeCell(
        double T, double p, double rho,
        double Y_O2, double Y_C3H8,
        double& qdot,
        double& omega_O2,
        double& omega_C3H8
    ) {
        qdot = 0.0;
        omega_O2 = 0.0;
        omega_C3H8 = 0.0;
        
        if (!enabled_) return;
        if (T < activationTemperature_) return;  // Too cold to ignite
        if (Y_O2 < 1e-6 || Y_C3H8 < 1e-6) return;  // Need both reactants
        
        // Call chemistry model
        std::vector<double> Y = {Y_O2, Y_C3H8};
        std::vector<double> omega(2);
        
        try {
            chemistry_.computeSourceTerms(T, p, rho, Y, omega, qdot);
            omega_O2 = omega[0];
            omega_C3H8 = omega[1];
        } catch (...) {
            // Chemistry failed, return zero source
            qdot = 0.0;
            omega_O2 = 0.0;
            omega_C3H8 = 0.0;
        }
    }
    
    void setEnabled(bool enable) { enabled_ = enable; }
    bool isEnabled() const { return enabled_; }
};

// [Rest of solver code would go here - for now this is a template]

int main(int argc, char** argv) {
    cout << "RFAM Tetra Solver WITH CHEMISTRY (Proof of Concept)" << endl;
    cout << "Chemistry: 1-step C3H8 + O2 → products" << endl;
    cout << "Activation temperature: " << ACTIVATION_T_MIN << " K" << endl;
    cout << "" << endl;
    
    // Parse command line for --enable-chemistry flag
    bool enableChemistry = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--enable-chemistry") == 0) {
            enableChemistry = true;
            cout << "Chemistry ENABLED by command line flag" << endl;
        }
    }
    
    if (!enableChemistry) {
        cout << "Chemistry DISABLED (use --enable-chemistry to enable)" << endl;
        cout << "Running cold flow / mixing mode only" << endl;
    }
    
    // Create chemistry integrator
    SimpleChemistryIntegrator chemistry(enableChemistry, ACTIVATION_T_MIN);
    
    // [Full solver implementation would follow here]
    // Key integration points:
    //
    // 1. In LUSGS loop after flux calculation:
    //    for (size_t c = 0; c < mesh.cellCount(); ++c) {
    //        double T = temperature[c];
    //        double p = pressure[c];
    //        double rho = density[c];
    //        double Y_O2 = rhoY_O2[c] / rho;
    //        double Y_C3H8 = rhoY_C3H8[c] / rho;
    //        
    //        double qdot, omega_O2, omega_C3H8;
    //        chemistry.computeCell(T, p, rho, Y_O2, Y_C3H8, 
    //                              qdot, omega_O2, omega_C3H8);
    //        
    //        // Apply source terms (explicit for now)
    //        double dt = timestep[c];
    //        double vol = volume[c];
    //        
    //        rhoE[c] += qdot * dt;
    //        rhoY_O2[c] += omega_O2 * dt;
    //        rhoY_C3H8[c] += omega_C3H8 * dt;
    //    }
    //
    // 2. Monitor chemistry statistics each step:
    //    int nReacting = 0;
    //    double totalHeatRelease = 0.0;
    //    for (each cell) {
    //        if (qdot[c] > 0) {
    //            nReacting++;
    //            totalHeatRelease += qdot[c] * volume[c];
    //        }
    //    }
    //    cout << "Reacting cells: " << nReacting << " / " << totalCells << endl;
    //    cout << "Total heat release: " << totalHeatRelease / 1e6 << " MW" << endl;
    
    cout << "TODO: Full solver implementation" << endl;
    cout << "This is a TEMPLATE showing integration points" << endl;
    
    return 0;
}
EOF

echo "Created solver_3d_tetra_chemistry.cpp template"
echo "Next: Copy full solver code and insert chemistry calls"
