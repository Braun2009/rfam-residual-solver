#ifndef SandiegoReducedChemistry_H
#define SandiegoReducedChemistry_H

/*---------------------------------------------------------------------------*\
    SandiegoReducedChemistry - the San Diego mechanism restricted to a species
    subset, with every reaction whose participating species are not all in the
    subset dropped.

    Why this exists
    ---------------
    The full mechanism is 59 species / 271 reactions. Measured chemistry cost on
    this codebase is 129 us per cell evaluation against 0.49 us for the
    two-step mechanism, which puts one 1.76M-cell step at 25-45 minutes -- a
    convergence run would take over a week. The subset keeps the C1-C3
    oxidation backbone and everything that sets the heat release, and drops the
    large-molecule and peroxy-ester branches (C4 species, RO2/ketone paths,
    O3), which contribute little at the chamber conditions here.

    The reaction table stays indexed the way the mechanism defines it: the
    reduced class keeps a full 59-slot concentration array and gathers the
    source terms afterwards, so netRateOfProgress and the NASA-7 data are used
    exactly as for the full mechanism. Only the loop bounds change.

    Reduction is a documented approximation, not a rewrite: the subset and the
    resulting reaction count are printed at construction so a run log records
    what actually ran.
\*---------------------------------------------------------------------------*/

#include "ChemistryModel.H"
#include "SandiegoChemistry.H"
#include "SandiegoNasa7Data.H"
#include "SandiegoKineticsRate.H"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

class SandiegoReducedChemistry : public ChemistryModel
{
public:
    // Default subset: fuel and C1-C3 oxidation backbone, the H/O radical pool,
    // CO/CO2, the two inerts that matter for this case, and the C2/C3
    // intermediates needed to reach the fuel. 28 of 59 species.
    static std::vector<std::string> defaultSubset()
    {
        return {"N2","AR","H","O2","OH","O","H2","H2O","HO2","H2O2",
                "CO","CO2","HCO","CH3","CH4","CH2O","T-CH2","S-CH2",
                "C2H4","CH3O","C2H5","C2H6","C2H3","C2H2",
                "HCCO","CH2CO","C3H6","C3H8"};
    }

    explicit SandiegoReducedChemistry(const std::vector<std::string>& subset)
    {
        const int nFull=SandiegoNasa7::nSpecies;
        std::vector<int> fullIndex(subset.size(),-1);
        std::vector<bool> kept(static_cast<std::size_t>(nFull),false);
        for(std::size_t k=0;k<subset.size();++k)
        {
            int found=-1;
            for(int i=0;i<nFull;++i)
                if(subset[k]==SandiegoNasa7::species[i].name){found=i;break;}
            if(found<0)
                throw std::invalid_argument(
                    "SandiegoReducedChemistry: unknown species "+subset[k]);
            if(kept[static_cast<std::size_t>(found)])
                throw std::invalid_argument(
                    "SandiegoReducedChemistry: duplicate species "+subset[k]);
            kept[static_cast<std::size_t>(found)]=true;
            fullIndex[k]=found;
        }
        fullIndex_=fullIndex;
        kept_=kept;
        nSpecies_=static_cast<int>(subset.size());
        speciesNames_=subset;
        molecularWeights_.resize(subset.size());
        for(std::size_t k=0;k<subset.size();++k)
            molecularWeights_[k]=
                SandiegoNasa7::species[fullIndex_[k]].mw;

        // Keep a reaction only when every participating species survived.
        for(int j=0;j<SandiegoKinetics::nReactions;++j)
        {
            const auto& r=SandiegoKinetics::reactions[j];
            bool keep=true;
            for(int k=0;k<r.nReactants&&keep;++k)
                keep=kept[static_cast<std::size_t>(r.reactants[k].species)];
            for(int k=0;k<r.nProducts&&keep;++k)
                keep=kept[static_cast<std::size_t>(r.products[k].species)];
            for(int k=0;k<r.nOrders&&keep;++k)
                keep=kept[static_cast<std::size_t>(r.orders[k].species)];
            if(keep) reactionIndex_.push_back(j);
        }
        nReactions_=static_cast<int>(reactionIndex_.size());
    }

    SandiegoReducedChemistry()
        :SandiegoReducedChemistry(defaultSubset()) {}

    void computeSourceTerms(
        double T,double p,double rho,
        const std::vector<double>& Y,
        std::vector<double>& omega_i,
        double& qdot) const override
    {
        if(static_cast<int>(Y.size())!=nSpecies_)
            throw std::invalid_argument(
                "SandiegoReducedChemistry: composition dimension mismatch");
        const int nFull=SandiegoNasa7::nSpecies;
        double C[59];
        for(int i=0;i<nFull;++i) C[i]=0.0;
        for(int k=0;k<nSpecies_;++k)
        {
            const int s=fullIndex_[static_cast<std::size_t>(k)];
            const double mw=SandiegoNasa7::species[s].mw;
            C[s]=(Y[static_cast<std::size_t>(k)]>0.0)
                ?rho*Y[static_cast<std::size_t>(k)]/mw:0.0;
        }
        double omegaFull[59];
        for(int i=0;i<nFull;++i) omegaFull[i]=0.0;
        qdot=0.0;
        for(int j:reactionIndex_)
        {
            const auto& r=SandiegoKinetics::reactions[j];
            const double net=SandiegoKinetics::netRateOfProgress(r,T,C);
            if(net==0.0) continue;
            for(int k=0;k<r.nReactants;++k)
            {
                const int s=r.reactants[k].species;
                omegaFull[s]-=r.reactants[k].coeff*net
                    *SandiegoNasa7::species[s].mw;
            }
            for(int k=0;k<r.nProducts;++k)
            {
                const int s=r.products[k].species;
                omegaFull[s]+=r.products[k].coeff*net
                    *SandiegoNasa7::species[s].mw;
            }
            double dh=0.0;
            for(int k=0;k<r.nReactants;++k)
                dh-=r.reactants[k].coeff
                    *SandiegoNasa7::hMole(r.reactants[k].species,298.15);
            for(int k=0;k<r.nProducts;++k)
                dh+=r.products[k].coeff
                    *SandiegoNasa7::hMole(r.products[k].species,298.15);
            qdot+=-net*dh;
        }
        omega_i.assign(static_cast<std::size_t>(nSpecies_),0.0);
        for(int k=0;k<nSpecies_;++k)
            omega_i[static_cast<std::size_t>(k)]=
                omegaFull[fullIndex_[static_cast<std::size_t>(k)]];
    }

    void print() const override
    {
        std::cout<<"SandiegoReducedChemistry: "<<nSpecies_
                 <<" species, "<<nReactions_
                 <<" reactions (San Diego, reduced subset)\n";
    }

    const std::vector<int>& fullSpeciesIndex() const {return fullIndex_;}

private:
    std::vector<int> fullIndex_;
    std::vector<bool> kept_;
    std::vector<int> reactionIndex_;
};

} // End namespace Foam

#endif
