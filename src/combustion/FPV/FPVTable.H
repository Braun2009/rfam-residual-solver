#ifndef FPVTable_H
#define FPVTable_H

#include <vector>
#include <string>
#include <map>

namespace Foam
{

class FPVTable
{
private:
    // Table dimensions
    std::vector<double> Z_grid_;   // Mixture fraction
    std::vector<double> C_grid_;   // Progress variable
    std::vector<double> h_grid_;   // Enthalpy
    
    int nZ_, nC_, nh_;
    
    // Table data (stored as flat arrays for efficiency)
    std::vector<double> T_;           // Temperature
    std::vector<double> rho_;         // Density
    std::vector<double> omega_dot_;   // Source term
    
    // Species mass fractions (13 species)
    std::vector<std::vector<double>> Y_;  // [species][grid_point]
    std::vector<std::string> species_names_;
    
    // Helper: convert (i,j,k) to flat index
    int index(int i, int j, int k) const {
        return i * nC_ * nh_ + j * nh_ + k;
    }
    
    // Trilinear interpolation
    double interpolate3D(
        const std::vector<double>& data,
        double Z, double C, double h
    ) const;
    
    // Find bracket indices
    void findBracket(
        const std::vector<double>& grid,
        double value,
        int& i0, int& i1,
        double& weight
    ) const;

public:
    FPVTable();
    
    // Load table from CSV file
    bool loadFromCSV(const std::string& filename);

    // DESIGN INVARIANT (2026-07-05): FPVTable has exactly ONE output method.
    // getTemperature/getDensity/getSpeciesMassFractions were REMOVED because they
    // violated the core RFAM principle "thermodynamics from RealFluidState only,
    // never from tables". If you need T/rho/Y, use RealFluidState::updatePTY or
    // updateRhoEY. The FPV table is a combustion source-term model, not a
    // thermodynamic state model. [[realfluidallmach-dev: FPV quarantine 2026-07-05]]
    double getOmegaDot(double Z, double C, double h) const;
    
    // Getters
    int nSpecies() const { return species_names_.size(); }
    const std::vector<std::string>& speciesNames() const { return species_names_; }
    
    // Bounds
    double Z_min() const { return Z_grid_.front(); }
    double Z_max() const { return Z_grid_.back(); }
    double C_min() const { return C_grid_.front(); }
    double C_max() const { return C_grid_.back(); }
    double h_min() const { return h_grid_.front(); }
    double h_max() const { return h_grid_.back(); }
    
    void print() const;
};

} // End namespace Foam

#endif
