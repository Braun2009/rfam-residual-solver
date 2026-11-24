#ifndef FixedConservedSourceOperator_H
#define FixedConservedSourceOperator_H

/*---------------------------------------------------------------------------*\
    FixedConservedSourceOperator - immutable conservative rate field.

    Used by the coarse FAS equation to add the fixed source S_H=-tau_H after
    the production transport/diffusion operators have assembled their update.
    The source owns no evolving state and advertises no artificial Jacobian or
    stability bound.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include <stdexcept>
#include <utility>
#include <vector>

namespace Foam
{

class FixedConservedSourceOperator : public Operator
{
    std::vector<ConservedVars> rate_;
    std::vector<unsigned char> active_;
    std::string name_;

    static void addScaled(
        ConservedVars& target,const ConservedVars& source,double scale)
    {
        target.rho+=scale*source.rho;
        target.rhoE+=scale*source.rhoE;
        for(std::size_t component=0;component<target.rhoU.size();++component)
            target.rhoU[component]+=scale*source.rhoU[component];
        for(std::size_t component=0;component<target.rhoY.size();++component)
            target.rhoY[component]+=scale*source.rhoY[component];
        target.rhoZ+=scale*source.rhoZ;
        target.rhoC+=scale*source.rhoC;
    }

public:
    explicit FixedConservedSourceOperator(
        std::vector<ConservedVars> rate,
        std::vector<unsigned char> active={},
        std::string name="fixed conserved FAS source")
        :rate_(std::move(rate)),active_(std::move(active)),
         name_(std::move(name))
    {
        if(!active_.empty()&&active_.size()!=rate_.size())
            throw std::invalid_argument(
                "fixed conserved source active mask size mismatch");
        for(unsigned char value:active_)
            if(value>1U)
                throw std::invalid_argument(
                    "fixed conserved source active mask is invalid");
    }

    void apply(Domain& domain,double dt) const override
    {
        if(rate_.size()!=static_cast<std::size_t>(domain.N()))
            throw std::runtime_error(
                "fixed conserved source dimensions do not match Domain");
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active_.empty()
                &&active_[static_cast<std::size_t>(cell)]==0U) continue;
            const ConservedVars& source=rate_[static_cast<std::size_t>(cell)];
            if(source.rhoY.size()!=static_cast<std::size_t>(domain.nSpecies()))
                throw std::runtime_error(
                    "fixed conserved source species layout mismatch");
            addScaled(domain.Unew(cell),source,dt);
        }
    }

    // ---- residual path ----
    // rate_ is dU/dt (what apply() adds). The residual is the left-hand side,
    //   dU/dt + R = 0,  so a source contributes R -= rate_.  See Residual.H.
    bool supportsResidual() const override { return true; }

    void residual(const Domain& domain,ResidualField& R) const override
    {
        if(rate_.size()!=static_cast<std::size_t>(domain.N()))
            throw std::runtime_error(
                "fixed conserved source dimensions do not match Domain");
        if(R.size()!=rate_.size())
            throw std::runtime_error(
                "fixed conserved source residual field size mismatch");
        for(int cell=0;cell<domain.N();++cell)
        {
            if(!active_.empty()
                &&active_[static_cast<std::size_t>(cell)]==0U) continue;
            const ConservedVars& source=rate_[static_cast<std::size_t>(cell)];
            if(source.rhoY.size()!=static_cast<std::size_t>(domain.nSpecies()))
                throw std::runtime_error(
                    "fixed conserved source species layout mismatch");
            addScaledRate(R[static_cast<std::size_t>(cell)],source,-1.0);
        }
    }

    double maxStableDt(const Domain&,double) const override
    { return 1.0e30; }

    void localStableDt(
        const Domain& domain,double,std::vector<double>& dt) const override
    {
        dt.assign(static_cast<std::size_t>(domain.N()),1.0e30);
    }

    const std::vector<ConservedVars>& rate() const { return rate_; }
    const std::vector<unsigned char>& active() const { return active_; }
    std::string name() const override { return name_; }
};

} // End namespace Foam

#endif
