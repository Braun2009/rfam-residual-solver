#ifndef FluxScheme_H
#define FluxScheme_H

/*---------------------------------------------------------------------------*\
    FluxScheme - abstract flux interface for the framework.

    AUSMFluxAdapter wraps the existing AUSMPlusUpFlux (UNCHANGED) so it plugs
    into the operator framework without modifying the verified module.
\*---------------------------------------------------------------------------*/

#include "ConservedVars.H"
#include "AcousticStateAdmission.H"
#include "EquilibriumThermoModel.H"
#include "AUSMPlusUpFlux.H"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

namespace Foam
{

class FluxScheme
{
public:
    virtual ~FluxScheme() = default;
    virtual void computeFlux(
        const PrimitiveVars& left,
        const PrimitiveVars& right,
        const std::vector<double>& normal,
        ConservedVars& flux
    ) const = 0;
    virtual std::string name() const = 0;
    // Optional interface-consistency hook (AUSM+up implements it; the default
    // is a no-op so other schemes are unaffected).
    virtual void setConsistentDensityHook(
        std::function<bool(double, double, double, double&)>) {}
    virtual std::size_t consistentDensityFallbacks() const { return 0; }
    virtual void setSpeciesNames(const std::vector<std::string>& speciesNames)
    { (void)speciesNames; }
    virtual void computeFluxCached(
        const PrimitiveVars& left,const PrimitiveVars& right,
        const std::vector<double>& normal,ConservedVars& flux,
        double rhoLeft,double soundLeft,double enthalpyLeft,
        double rhoRight,double soundRight,double enthalpyRight) const = 0;

    void computeFluxFromEquilibrium(
        const PrimitiveVars& left,
        const thermo::EquilibriumState& leftThermo,
        const PrimitiveVars& right,
        const thermo::EquilibriumState& rightThermo,
        const std::vector<double>& normal,
        ConservedVars& flux) const
    {
        if (!solverConsumerAdmitsAcousticState(leftThermo)
            || !solverConsumerAdmitsAcousticState(rightThermo))
            throw std::runtime_error(
                "FluxScheme: acoustic state is not admitted by the "
                "production solver");

        if (!std::isfinite(leftThermo.rho) || leftThermo.rho <= 0.0
            || !std::isfinite(rightThermo.rho) || rightThermo.rho <= 0.0
            || !std::isfinite(leftThermo.h)
            || !std::isfinite(rightThermo.h))
            throw std::runtime_error(
                "FluxScheme: invalid equilibrium flux payload");

        computeFluxCached(
            left,right,normal,flux,
            leftThermo.rho,leftThermo.acoustic.soundSpeed(),leftThermo.h,
            rightThermo.rho,rightThermo.acoustic.soundSpeed(),rightThermo.h);
    }
};

class AUSMFluxAdapter : public FluxScheme
{
private:
    AUSMPlusUpFlux ausm_;
public:
    AUSMFluxAdapter()=default;
    explicit AUSMFluxAdapter(const std::vector<std::string>& speciesNames)
        :ausm_(speciesNames) {}
    void setSpeciesNames(const std::vector<std::string>& speciesNames) override
    { ausm_.setSpeciesNames(speciesNames); }
    void setConsistentDensityHook(
        std::function<bool(double, double, double, double&)> hook) override
    { ausm_.setConsistentDensityHook(std::move(hook)); }
    std::size_t consistentDensityFallbacks() const override
    { return ausm_.consistentDensityFallbacks(); }
    void computeFlux(
        const PrimitiveVars& left,
        const PrimitiveVars& right,
        const std::vector<double>& normal,
        ConservedVars& flux
    ) const override
    {
        ausm_.computeFlux(left, right, normal, flux);
    }
    std::string name() const override { return "AUSM+up"; }
    void computeFluxCached(
        const PrimitiveVars& left,const PrimitiveVars& right,
        const std::vector<double>& normal,ConservedVars& flux,
        double rhoLeft,double soundLeft,double enthalpyLeft,
        double rhoRight,double soundRight,double enthalpyRight) const override
    {
        ausm_.computeFluxCached(left,right,normal,flux,rhoLeft,soundLeft,enthalpyLeft,
                                rhoRight,soundRight,enthalpyRight);
    }
};

} // End namespace Foam

#endif
