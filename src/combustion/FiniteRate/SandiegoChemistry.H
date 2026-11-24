/* SandiegoChemistry - 59-species San Diego mechanism as an RFAM
   ChemistryModel. Source terms and heat release are evaluated from the
   generated SandiegoKineticsData.H / SandiegoKineticsRate.H tables, which
   are cross-checked against Cantera 3.1 for all 271 reactions. */

#ifndef SandiegoChemistry_H
#define SandiegoChemistry_H

#include "ChemistryModel.H"
#include "SandiegoKineticsRate.H"
#include "SandiegoNasa7Data.H"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

class SandiegoChemistry : public ChemistryModel
{
public:
    SandiegoChemistry()
    {
        nSpecies_ = SandiegoNasa7::nSpecies;
        nReactions_ = SandiegoKinetics::nReactions;
        speciesNames_.resize(nSpecies_);
        molecularWeights_.resize(nSpecies_);
        for (int i = 0; i < nSpecies_; ++i)
        {
            speciesNames_[i] = SandiegoNasa7::species[i].name;
            // ChemistryModel stores kg/kmol; San Diego table is kg/mol.
            molecularWeights_[i] =
                SandiegoNasa7::species[i].mw * 1000.0;
        }
    }

    void computeSourceTerms(
        double T, double p, double rho,
        const std::vector<double>& Y,
        std::vector<double>& omega_i,
        double& qdot) const override
    {
        if (static_cast<int>(Y.size()) != nSpecies_)
            throw std::invalid_argument(
                "SandiegoChemistry: composition dimension mismatch");
        omega_i.assign(nSpecies_, 0.0);
        double C[59];
        for (int i = 0; i < nSpecies_; ++i)
        {
            const double mw = SandiegoNasa7::species[i].mw;  // kg/mol
            C[i] = (Y[i] > 0.0) ? rho * Y[i] / mw : 0.0;     // mol/m3
        }
        qdot = 0.0;
        for (int j = 0; j < nReactions_; ++j)
        {
            const SandiegoKinetics::SandiegoReactionData& r =
                SandiegoKinetics::reactions[j];
            const double net = SandiegoKinetics::netRateOfProgress(r, T, C);
            if (net == 0.0) continue;
            for (int k = 0; k < r.nReactants; ++k)
            {
                const int s = r.reactants[k].species;
                omega_i[s] -= r.reactants[k].coeff * net
                    * SandiegoNasa7::species[s].mw;
            }
            for (int k = 0; k < r.nProducts; ++k)
            {
                const int s = r.products[k].species;
                omega_i[s] += r.products[k].coeff * net
                    * SandiegoNasa7::species[s].mw;
            }
            // Heat release: -net * DH, DH = sum nu_i hf_i, hf = h(298.15).
            double dh = 0.0;
            for (int k = 0; k < r.nReactants; ++k)
                dh -= r.reactants[k].coeff
                    * SandiegoNasa7::hMole(r.reactants[k].species, 298.15);
            for (int k = 0; k < r.nProducts; ++k)
                dh += r.products[k].coeff
                    * SandiegoNasa7::hMole(r.products[k].species, 298.15);
            qdot += -net * dh;   // J/m3/s
        }
    }

    void print() const override
    {
        std::cout << "SandiegoChemistry: " << nSpecies_
                  << " species, " << nReactions_
                  << " reactions (San Diego mechanism)\n";
    }
};

} // End namespace Foam

#endif
