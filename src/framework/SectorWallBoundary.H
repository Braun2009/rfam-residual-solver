#ifndef SectorWallBoundary_H
#define SectorWallBoundary_H

#include "Domain.H"
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorWallBoundary
{
    double temperature_;

    void applyImpl(Domain& d,bool inStepTransaction) const
    {
        if(!d.stateReady()||d.tripped()
            ||d.stepTransactionActive()!=inStepTransaction)
            throw std::logic_error(
                inStepTransaction
                    ?"SectorWallBoundary inner update requires an active transaction"
                    :"SectorWallBoundary requires a ready Domain outside a transaction");
        std::vector<Domain::CompleteCellReplacement> replacements;
        replacements.reserve(static_cast<std::size_t>(d.Nx()*d.Nz()));
        for(int i=0;i<d.Nx();++i)
        for(int k=0;k<d.Nz();++k)
        {
            PrimitiveVars ghost=d.W(i,d.Ny()-2,k);
            ghost.T=temperature_;
            for(double& velocity:ghost.U) velocity=-velocity;
            const thermo::EquilibriumResult initialized=
                d.thermoModel().initializePTY(
                    ghost.p,ghost.T,ghost.Y);
            if(!initialized.success())
                throw std::runtime_error(
                    "SectorWallBoundary PTY recovery failed");
            replacements.push_back(Domain::CompleteCellReplacement{
                i,d.Ny()-1,k,ghost,*initialized.state});
        }
        const bool published=inStepTransaction
            ?d.replaceCompleteCellsInStepTransaction(replacements)
            :d.replaceCompleteCells(replacements);
        if(!published)
            throw std::runtime_error(
                "SectorWallBoundary publication failed: "+d.failure().detail);
    }

public:
    explicit SectorWallBoundary(double temperature):temperature_(temperature) {}

    void apply(Domain& d) const
    { applyImpl(d,false); }

    void applyInStepTransaction(Domain& d) const
    { applyImpl(d,true); }
};

} // End namespace Foam
#endif
