#ifndef ChemistrySourceOperator_H
#define ChemistrySourceOperator_H

/*---------------------------------------------------------------------------*\
    ChemistrySourceOperator - finite-rate chemistry as a framework Operator.

    Phase 7. Thin wrapper over FiniteRateReactor. Per cell:
      - read global (rho, p, T, Y) from the Domain (t^n state)
      - 0-D chem-space CVODE sub-step (via FiniteRateReactor)
      - ADDITIVE increment of the chem species' global rhoY into Unew
        (only chem slots change; convection's rhoY update is preserved)
      - EXPLICIT chemical heat: Unew.rhoE += dE   (sensible-energy convention)
    T/p are NOT set by chemistry; they come from the EOS reback in
    Domain.commit (化学不越权). Operates on ALL cells (chemistry is local).

    STEP 2 (full 1D reacting flow) PREREQUISITE: the EOS used by
    StateConversion/Domain.commit must be species-aware over the global
    SpeciesRegistry set (RealFluidState::setSpecies). Until then this operator
    is validated only via the standalone 0-D reactor (test_phase7_chem0d).
\*---------------------------------------------------------------------------*/

#include "Operator.H"              // framework (-Isrc/framework)
#include "FiniteRateReactor.H"     // (-Isrc/combustion)
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

class ChemistrySourceOperator : public Operator
{
private:
    ChemistryModel* chemistry_ = nullptr;
    std::vector<std::unique_ptr<FiniteRateReactor>> reactors_;
    std::shared_ptr<const thermo::EquilibriumThermoModel> boundThermoOwner_;
    std::string boundThermoModelId_;
    std::vector<std::string> boundThermoSpeciesNames_;
    EnergyConvention boundEnergyConvention_ =
        EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;

    static int reactorCount()
    {
        int count = 1;
#ifdef _OPENMP
        count = omp_get_max_threads();
#endif
        return count;
    }

    void bindThermoOwner(
        const SpeciesRegistry& registry,
        std::shared_ptr<const thermo::EquilibriumThermoModel> owner)
    {
        if (!owner)
            throw std::invalid_argument(
                "ChemistrySourceOperator: null thermo owner");
        try
        {
            boundThermoModelId_ = owner->modelId();
            boundThermoSpeciesNames_ = owner->speciesNames();
            boundEnergyConvention_ = owner->energyConvention();
        }
        catch (...)
        {
            throw std::invalid_argument(
                "ChemistrySourceOperator: unreadable thermo owner identity");
        }
        if (boundThermoModelId_.empty()
            || boundThermoModelId_.find_first_not_of(" \t\r\n")
                == std::string::npos
            || boundThermoSpeciesNames_ != registry.names()
            || boundEnergyConvention_ != EnergyConvention::
                SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
            throw std::invalid_argument(
                "ChemistrySourceOperator: incompatible thermo owner identity");
        boundThermoOwner_ = std::move(owner);
    }

    bool domainUsesBoundThermoOwner(const Domain& domain) const noexcept
    {
        if (!boundThermoOwner_) return true;
        try
        {
            const auto domainOwner = domain.thermoModelOwner();
            if (!domainOwner
                || domainOwner.get() != boundThermoOwner_.get()
                || boundThermoOwner_->modelId() != boundThermoModelId_
                || boundThermoOwner_->speciesNames()
                    != boundThermoSpeciesNames_
                || boundThermoOwner_->energyConvention()
                    != boundEnergyConvention_
                || domain.thermoModelId() != boundThermoModelId_
                || domain.thermoModel().speciesNames()
                    != boundThermoSpeciesNames_
                || domain.speciesNames() != boundThermoSpeciesNames_
                || domain.energyConvention() != boundEnergyConvention_)
                return false;
            for (const auto& reactor : reactors_)
                if (!reactor || !reactor->hasInjectedThermoModel()
                    || reactor->boundThermoModelId()
                        != boundThermoModelId_)
                    return false;
            return true;
        }
        catch (...)
        {
            return false;
        }
    }

    static double estimateFromCommittedThermo(
        ChemistryModel* chemistry,
        const SpeciesMapper& mapper,
        double density,
        double pressure,
        double temperature,
        double internalEnergy,
        const std::vector<double>& globalMassFractions,
        double safetyFactor)
    {
        if (!chemistry
            || !std::isfinite(density) || density <= 0.0
            || !std::isfinite(pressure) || pressure <= 0.0
            || !std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(internalEnergy)
            || !std::isfinite(safetyFactor) || safetyFactor <= 0.0
            || globalMassFractions.empty())
            return 0.0;
        for (double massFraction : globalMassFractions)
            if (!std::isfinite(massFraction) || massFraction < 0.0)
                return 0.0;

        std::vector<double> chemistryY(
            static_cast<std::size_t>(mapper.nChem()), 0.0);
        std::vector<double> omega;
        mapper.mapGlobalToChem(globalMassFractions, chemistryY);
        double heatRelease = 0.0;
        try
        {
            chemistry->computeSourceTerms(
                temperature, pressure, density,
                chemistryY, omega, heatRelease);
        }
        catch (...)
        {
            return 0.0;
        }
        if (!std::isfinite(heatRelease)
            || omega.size() != chemistryY.size())
            return 0.0;

        double timeScale = 1.0e30;
        if (std::abs(heatRelease) > 1.0)
            timeScale = std::min(
                timeScale,
                density*std::max(std::abs(internalEnergy), 1.0)
                    /std::abs(heatRelease));
        for (std::size_t species = 0;
             species < chemistryY.size(); ++species)
        {
            if (!std::isfinite(omega[species])) return 0.0;
            if (omega[species] < 0.0 && chemistryY[species] > 0.0)
                timeScale = std::min(
                    timeScale,
                    density*chemistryY[species]/(-omega[species]));
        }
        return timeScale >= 1.0e29
            ? 1.0e30 : safetyFactor*timeScale;
    }

    int threadIndex() const
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    FiniteRateReactor& reactorForThread() const
    {
        return *reactors_[threadIndex()];
    }

public:
    ChemistrySourceOperator(ChemistryModel* chem, const SpeciesRegistry& reg)
        : chemistry_(chem)
    {
        const int count = reactorCount();
        reactors_.reserve(count);
        for (int i = 0; i < count; ++i)
            reactors_.emplace_back(new FiniteRateReactor(chem, reg));
    }

    ChemistrySourceOperator(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoOwner)
        : chemistry_(chemistry)
    {
        bindThermoOwner(registry, std::move(thermoOwner));
        const int count = reactorCount();
        reactors_.reserve(count);
        for (int i = 0; i < count; ++i)
            reactors_.emplace_back(new FiniteRateReactor(
                chemistry, registry, boundThermoOwner_));
    }

    bool hasBoundThermoOwner() const
    { return static_cast<bool>(boundThermoOwner_); }
    const std::string& boundThermoModelId() const
    { return boundThermoModelId_; }

    void apply(Domain& d, double dt) const override
    {
        if (!domainUsesBoundThermoOwner(d))
            throw std::runtime_error(
                "ChemistrySourceOperator: Domain thermo owner mismatch");
        const SpeciesMapper& m = reactors_.front()->mapper();
        int failures = 0;
        #pragma omp parallel for reduction(+:failures) schedule(guided)
        for (int i = 0; i < d.N(); ++i)
        {
            FiniteRateReactor& reactor = reactorForThread();
            const PrimitiveVars& w = d.W(i);
            double rho = d.U(i).rho;
            double T = w.T, p = w.p;

            std::vector<double> Yg = w.Y;     // global Y (copy of t^n)
            std::vector<double> Yold = Yg;
            if (reactor.estimateConservativeMaxDt(
                rho,p,T,Yg,1.0) >= 1.0e29) continue;

            double kineticEnergy = 0.0;
            for (double velocity : w.U)
                kineticEnergy += 0.5*velocity*velocity;
            const double specificEnergy =
                d.U(i).rhoE/rho - kineticEnergy;
            double dE = 0.0;
            double finalPressure = p;
            const bool chemistrySuccess = reactor.reactConservative(
                rho,specificEnergy,Yg,dt,dE,T,finalPressure);
            if (!chemistrySuccess)
            {
                ++failures;
                continue;
            }

            ConservedVars& un = d.Unew(i);
            double couplingScale = 1.0;
            std::vector<double> increments(m.nChem(),0.0);
            for (int j = 0; j < m.nChem(); ++j)
            {
                int g = m.globalIndex(j);
                increments[j] = rho*(Yg[g]-Yold[g]);
                if (increments[j] < 0.0 && un.rhoY[g]+increments[j] < 0.0)
                    couplingScale = std::min(
                        couplingScale,un.rhoY[g]/(-increments[j]));
            }
            couplingScale = std::clamp(couplingScale,0.0,1.0);
            if (couplingScale < 1.0) couplingScale *= 1.0-1.0e-6;
            for (int j = 0; j < m.nChem(); ++j)
                un.rhoY[m.globalIndex(j)] += couplingScale*increments[j];
            un.rhoE += couplingScale*dE;
        }
        if (failures > 0)
            throw std::runtime_error(
                "ChemistrySourceOperator: conservative chemistry integration failed");
    }

    // ⑤ chemistry dt limiter: scan all cells, return the minimum chemistry
    // timescale scaled by safetyFactor. SolverBase takes the min over all
    // operators (convective CFL and chemistry), so a stiff igniting cell
    // automatically tightens the global dt without touching the framework.
    //
    // The limiter uses the same real-fluid e/Cv state as reactConservative.
    double maxStableDt(const Domain& d, double /*CFL*/) const override
    {
        if (!domainUsesBoundThermoOwner(d))
            throw std::runtime_error(
                "ChemistrySourceOperator: Domain thermo owner mismatch");
        const double safetyFactor = 0.1;
        double dtMin = 1e30;
        #pragma omp parallel for reduction(min:dtMin) schedule(static)
        for (int i = 0; i < d.N(); ++i)
        {
            FiniteRateReactor& reactor = reactorForThread();
            const PrimitiveVars& w = d.W(i);
            double dt_i = reactor.estimateConservativeMaxDt(
                d.U(i).rho, w.p, w.T, w.Y, safetyFactor);
            const thermo::EquilibriumState& equilibrium = d.Q(i);
            if (boundThermoOwner_ && dt_i <= 0.0
                && !equilibrium.thermodynamics.has_value()
                && (equilibrium.region
                        == thermo::EquilibriumRegion::SINGLE_VAPOR
                    || equilibrium.region
                        == thermo::EquilibriumRegion::SUPERCRITICAL))
                dt_i = estimateFromCommittedThermo(
                    chemistry_, reactor.mapper(),
                    d.U(i).rho, equilibrium.p, equilibrium.T,
                    equilibrium.e, equilibrium.overallMassFractions,
                    safetyFactor);
            if (dt_i < dtMin) dtMin = dt_i;
        }
        return dtMin;
    }

    std::string name() const override { return "ChemistrySource (finite-rate CVODE)"; }
};

} // End namespace Foam

#endif
