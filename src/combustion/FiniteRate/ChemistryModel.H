#ifndef ChemistryModel_H
#define ChemistryModel_H

#include <vector>
#include <string>
#include <iostream>

namespace Foam
{

/*---------------------------------------------------------------------------*\
    ChemistryModel - Abstract base class for finite-rate chemistry
    
    Design principle:
    - Chemistry ONLY provides source terms (omega_i, qdot)
    - Thermodynamics (T, rho, p) is handled by RealFluidState
    - Energy equation updates temperature via conservation
\*---------------------------------------------------------------------------*/

class ChemistryModel
{
protected:
    int nSpecies_;
    int nReactions_;
    std::vector<std::string> speciesNames_;
    std::vector<double> molecularWeights_;  // kg/kmol
    
    static constexpr double R_UNIVERSAL = 8314.0;  // J/(kmol·K)

public:
    ChemistryModel() : nSpecies_(0), nReactions_(0) {}
    virtual ~ChemistryModel() = default;
    
    // Core interface: compute source terms
    // omega_i [kg/m³/s] - species production rate
    // qdot [J/m³/s] - volumetric heat release rate
    virtual void computeSourceTerms(
        double T,                          // Temperature [K]
        double p,                          // Pressure [Pa]
        double rho,                        // Density [kg/m³]
        const std::vector<double>& Y,     // Mass fractions [-]
        std::vector<double>& omega_i,     // OUT: species rates
        double& qdot                       // OUT: heat release
    ) const = 0;
    
    // Auxiliary interfaces
    virtual int nSpecies() const { return nSpecies_; }
    virtual int nReactions() const { return nReactions_; }
    
    virtual const std::vector<std::string>& speciesNames() const { 
        return speciesNames_; 
    }
    
    virtual const std::vector<double>& molecularWeights() const {
        return molecularWeights_;
    }
    
    // Get species index by name
    virtual int speciesIndex(const std::string& name) const {
        for (int i = 0; i < nSpecies_; i++) {
            if (speciesNames_[i] == name) return i;
        }
        return -1;
    }
    
    // Print mechanism info
    virtual void print() const {
        std::cout << "Chemistry Mechanism:" << std::endl;
        std::cout << "  Species: " << nSpecies_ << std::endl;
        std::cout << "  Reactions: " << nReactions_ << std::endl;
    }
};

} // End namespace Foam

#endif
