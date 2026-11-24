#ifndef SpeciesMapper_H
#define SpeciesMapper_H

/*---------------------------------------------------------------------------*\
    SpeciesMapper - bridge between a chemistry model's local species space
    and the global SpeciesRegistry ordering.

    Stage 1.5 (Species/Energy Consistency). A ChemistryModel works in a small
    "chem-space" (e.g. TwoStepCO: {C3H8,O2,CO,CO2,H2O,N2}, 6 species). The CFD
    state / EOS work in "global-space" (SpeciesRegistry, 9 species). CVODE
    integrates the small chem-space ODE; this mapper moves mass fractions and
    source terms between the two spaces WITHOUT any implicit ordering assumption.

        chemToGlobal_[j] = global index of chem species j   (built once)

    Heat release stays a scalar qdot (sensible-energy convention) and is added
    to rhoE explicitly by the caller; the mapper only routes species arrays.
\*---------------------------------------------------------------------------*/

#include "SpeciesRegistry.H"
#include <string>
#include <vector>
#include <stdexcept>

namespace Foam
{

class SpeciesMapper
{
private:
    const SpeciesRegistry& reg_;
    std::vector<int> chemToGlobal_;   // chem index -> global index
    int nChem_, nGlobal_;

public:
    SpeciesMapper(const SpeciesRegistry& reg,
                  const std::vector<std::string>& chemNames)
        : reg_(reg), nChem_((int)chemNames.size()), nGlobal_(reg.size())
    {
        chemToGlobal_.resize(nChem_);
        for (int j = 0; j < nChem_; ++j)
        {
            int g = reg_.index(chemNames[j]);
            if (g < 0)
                throw std::runtime_error(
                    "SpeciesMapper: chem species '" + chemNames[j] +
                    "' not in SpeciesRegistry");
            chemToGlobal_[j] = g;
        }
    }

    int nChem()   const { return nChem_; }
    int nGlobal() const { return nGlobal_; }
    int globalIndex(int chemJ) const { return chemToGlobal_[chemJ]; }

    // chem-space mass fractions -> global-space (unlisted globals = 0)
    void mapChemToGlobal(const std::vector<double>& Ychem,
                         std::vector<double>& Yglobal) const
    {
        Yglobal.assign(nGlobal_, 0.0);
        for (int j = 0; j < nChem_; ++j) Yglobal[chemToGlobal_[j]] = Ychem[j];
    }

    // global-space mass fractions -> chem-space (gather)
    void mapGlobalToChem(const std::vector<double>& Yglobal,
                         std::vector<double>& Ychem) const
    {
        Ychem.assign(nChem_, 0.0);
        for (int j = 0; j < nChem_; ++j) Ychem[j] = Yglobal[chemToGlobal_[j]];
    }

    // chem-space source terms -> global-space (scatter; same routing as Y)
    void mapOmega(const std::vector<double>& omegaChem,
                  std::vector<double>& omegaGlobal) const
    {
        omegaGlobal.assign(nGlobal_, 0.0);
        for (int j = 0; j < nChem_; ++j)
            omegaGlobal[chemToGlobal_[j]] = omegaChem[j];
    }
};

} // End namespace Foam

#endif
