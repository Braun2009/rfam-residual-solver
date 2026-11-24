#ifndef RealFluidCVODEWrapper_H
#define RealFluidCVODEWrapper_H

#include "ChemistryModel.H"
#include "EquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRPhaseThermoOwnerAdapter.H"
#include "RealFluidState.H"
#include "RhsCompositionPreparation.H"
#include "SpeciesMapper.H"
#include "SpeciesRegistry.H"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_types.h>

namespace Foam
{

enum class RealFluidCVODEFailureClass
{
    NONE,
    RETRYABLE_NUMERICAL,
    PERMANENT_ADMISSION
};

struct RealFluidCVODEPhaseState
{
    double rho = std::numeric_limits<double>::quiet_NaN();
    double specificEnergy = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> massFractions;
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double pressure = std::numeric_limits<double>::quiet_NaN();
};

class RealFluidCVODEWrapper
{
private:
    ChemistryModel* chemistry_;
    const SpeciesRegistry& registry_;
    const SpeciesMapper& mapper_;
    int nChem_;
    int nEq_;

    double rho_;
    double rtol_;
    double speciesAtol_;
    double energyAtol_;
    double temperatureSeed_;
    std::vector<double> globalTemplate_;
    std::vector<double> rhsChemY_;
    std::vector<double> rhsGlobalY_;
    std::vector<double> rhsOmega_;
    RealFluidState thermo_;
    std::shared_ptr<const thermo::EquilibriumThermoModel> equilibriumThermo_;
    std::string boundThermoModelId_;
    std::vector<std::string> boundThermoSpeciesNames_;
    EnergyConvention boundEnergyConvention_ =
        EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    std::unique_ptr<phaseEquilibrium::PRPhaseThermoOwnerAdapter>
        phaseThermoOwner_;
    std::optional<phaseEquilibrium::PRPhaseThermoParentView> phaseParent_;
    std::optional<phaseEquilibrium::PRPhaseThermoLiveToken> phaseToken_;
    std::optional<RealFluidCVODEPhaseState> phaseInitialState_;
    double initializedTime_ = 0.0;
    bool readyForAdvance_ = false;
    bool phaseModeActive_ = false;
    RealFluidCVODEFailureClass lastFailureClass_ =
        RealFluidCVODEFailureClass::NONE;
    phaseEquilibrium::PRPhaseThermoOwnerStatus lastPhaseOwnerStatus_ =
        phaseEquilibrium::PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION;
    phaseEquilibrium::PRPhaseDisplacementStatus lastPhaseDisplacementStatus_ =
        phaseEquilibrium::PRPhaseDisplacementStatus::INVALID_CONFIGURATION;
    std::optional<phaseEquilibrium::PRPhaseDisplacementDiagnostics>
        lastPhaseDisplacementDiagnostics_;
    std::uint64_t phaseProbeCount_ = 0U;
    std::uint64_t phaseAcceptCount_ = 0U;

    SUNContext sunctx_;
    void* cvodeMem_;
    N_Vector y_;
    N_Vector absoluteTolerances_;
    SUNMatrix matrix_;
    SUNLinearSolver linearSolver_;
    bool initialized_;

    static bool admittedChemistryRegion(thermo::EquilibriumRegion region)
    {
        return region == thermo::EquilibriumRegion::SINGLE_VAPOR
            || region == thermo::EquilibriumRegion::SUPERCRITICAL;
    }

    bool injectedThermoIdentityMatches() const noexcept
    {
        if (!equilibriumThermo_) return false;
        try
        {
            return !boundThermoModelId_.empty()
                && !boundThermoSpeciesNames_.empty()
                && equilibriumThermo_->modelId() == boundThermoModelId_
                && equilibriumThermo_->speciesNames()
                    == boundThermoSpeciesNames_
                && equilibriumThermo_->energyConvention()
                    == boundEnergyConvention_
                && boundEnergyConvention_ == EnergyConvention::
                    SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
        }
        catch (...)
        {
            return false;
        }
    }

    void closePhaseTransaction() noexcept
    {
        if (phaseThermoOwner_) phaseThermoOwner_->abort();
        phaseParent_.reset();
        phaseToken_.reset();
        phaseInitialState_.reset();
        phaseModeActive_ = false;
    }

    bool permanentPhaseFailure() noexcept
    {
        lastFailureClass_ =
            RealFluidCVODEFailureClass::PERMANENT_ADMISSION;
        readyForAdvance_ = false;
        closePhaseTransaction();
        return false;
    }

    static bool permanentPhaseOwnerStatus(
        phaseEquilibrium::PRPhaseThermoOwnerStatus status) noexcept
    {
        using Status = phaseEquilibrium::PRPhaseThermoOwnerStatus;
        switch (status)
        {
            case Status::INVALID_CONFIGURATION:
            case Status::INVALID_PARENT_ADMISSION:
            case Status::INVALID_PARENT_STATE_CONTRACT:
            case Status::INVALID_PARENT_PHASE_ADMISSION:
            case Status::TRANSACTION_ALREADY_ACTIVE:
            case Status::NO_ACTIVE_TRANSACTION:
            case Status::PARENT_BINDING_DRIFT:
            case Status::STALE_TOKEN:
            case Status::DISPLACEMENT_NOT_AUTHORIZED:
            case Status::SEQUENCE_OVERFLOW:
                return true;
            case Status::SUCCESS:
            case Status::INVALID_REQUEST:
            case Status::DISPLACEMENT_AUTHORITY_FAILED:
            case Status::KERNEL_FAILURE:
            case Status::RECOVERY_CONTRACT_FAILED:
                return false;
        }
        return true;
    }

    bool probePhaseThermoState(
        double density,
        double internalEnergy,
        const std::vector<double>& massFractions,
        double temperatureGuess,
        double& temperature,
        double& pressure)
    {
        if (!phaseModeActive_ || !phaseThermoOwner_ || !phaseParent_
            || !phaseToken_ || !phaseThermoOwner_->transactionActive())
            return permanentPhaseFailure();

        phaseEquilibrium::PRPhaseThermoOwnerRequest request;
        request.phaseDensity = density;
        request.phaseSpecificInternalEnergy = internalEnergy;
        request.phaseMassFractions = massFractions;
        if (phaseParent_->committed)
            request.pressureSeed = phaseParent_->committed->Q.p;
        request.temperatureSeed = temperatureGuess;
        const phaseEquilibrium::PRPhaseThermoProbeResult recovered =
            phaseThermoOwner_->probe(
                *phaseParent_, request, *phaseToken_);
        ++phaseProbeCount_;
        lastPhaseOwnerStatus_ = recovered.status;
        lastPhaseDisplacementStatus_ = recovered.diagnostics.displacementStatus;
        lastPhaseDisplacementDiagnostics_ =
            recovered.diagnostics.displacementDiagnostics;
        if (!recovered.success() || !recovered.phaseState)
        {
            if (permanentPhaseOwnerStatus(recovered.status))
                return permanentPhaseFailure();
            lastFailureClass_ =
                RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            return false;
        }

        const double candidateTemperature = recovered.phaseState->temperature;
        const double candidatePressure = recovered.phaseState->pressure;
        if (!std::isfinite(candidateTemperature)
            || candidateTemperature <= 0.0
            || !std::isfinite(candidatePressure) || candidatePressure <= 0.0)
        {
            lastFailureClass_ =
                RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            return false;
        }

        temperature = candidateTemperature;
        pressure = candidatePressure;
        return true;
    }

    bool acceptPhaseThermoState(
        double density,
        double internalEnergy,
        const std::vector<double>& massFractions,
        double temperatureGuess,
        double& temperature,
        double& pressure)
    {
        if (!phaseModeActive_ || !phaseThermoOwner_ || !phaseParent_
            || !phaseToken_ || !phaseThermoOwner_->transactionActive())
            return permanentPhaseFailure();

        phaseEquilibrium::PRPhaseThermoOwnerRequest request;
        request.phaseDensity = density;
        request.phaseSpecificInternalEnergy = internalEnergy;
        request.phaseMassFractions = massFractions;
        if (phaseParent_->committed)
            request.pressureSeed = phaseParent_->committed->Q.p;
        request.temperatureSeed = temperatureGuess;
        const phaseEquilibrium::PRPhaseThermoOwnerResult recovered =
            phaseThermoOwner_->recover(
                *phaseParent_, request, *phaseToken_);
        ++phaseAcceptCount_;
        lastPhaseOwnerStatus_ = recovered.status;
        lastPhaseDisplacementStatus_ = recovered.diagnostics.displacementStatus;
        lastPhaseDisplacementDiagnostics_ =
            recovered.diagnostics.displacementDiagnostics;
        if (!recovered.success() || !recovered.phaseState
            || !recovered.acceptedToken)
        {
            const bool permanent = permanentPhaseOwnerStatus(
                recovered.status);
            if (permanent) return permanentPhaseFailure();
            lastFailureClass_ =
                RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }

        const double candidateTemperature = recovered.phaseState->temperature;
        const double candidatePressure = recovered.phaseState->pressure;
        if (!std::isfinite(candidateTemperature)
            || candidateTemperature <= 0.0
            || !std::isfinite(candidatePressure) || candidatePressure <= 0.0)
        {
            lastFailureClass_ =
                RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }

        phaseToken_ = *recovered.acceptedToken;
        temperature = candidateTemperature;
        pressure = candidatePressure;
        return true;
    }

    bool admittedInjectedResult(
        const thermo::EquilibriumResult& result,
        const std::vector<double>& overallMassFractions) const
    {
        if (!injectedThermoIdentityMatches()
            || result.diagnostics.modelId != boundThermoModelId_
            || !result.success()
            || !result.state.has_value()
            || !admittedChemistryRegion(result.state->region)
            || result.state->overallMassFractions.size()
                != overallMassFractions.size()
            || !equilibriumThermo_->stateValidForModel(*result.state))
            return false;

        constexpr double compositionTolerance = 1.0e-8;
        for (std::size_t species = 0;
             species < overallMassFractions.size(); ++species)
        {
            const double expected = overallMassFractions[species];
            const double recovered =
                result.state->overallMassFractions[species];
            if (!std::isfinite(expected) || !std::isfinite(recovered)
                || std::abs(recovered - expected) > compositionTolerance)
                return false;
        }
        return true;
    }

    void bindInjectedThermoModel(
        std::shared_ptr<const thermo::EquilibriumThermoModel> model)
    {
        if (!model)
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: null equilibrium thermo model");

        try
        {
            boundThermoModelId_ = model->modelId();
            boundThermoSpeciesNames_ = model->speciesNames();
            boundEnergyConvention_ = model->energyConvention();
        }
        catch (...)
        {
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: unreadable thermo model identity");
        }

        if (boundThermoModelId_.empty()
            || boundThermoModelId_.find_first_not_of(" \t\r\n")
                == std::string::npos
            || boundThermoSpeciesNames_.empty()
            || boundThermoSpeciesNames_ != registry_.names()
            || boundEnergyConvention_ != EnergyConvention::
                SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: incompatible thermo model identity");

        equilibriumThermo_ = std::move(model);
        if (!injectedThermoIdentityMatches())
        {
            equilibriumThermo_.reset();
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: unstable thermo model identity");
        }
    }

    bool recoverThermoState(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        double temperatureGuess,
        double& temperature,
        double& pressure)
    {
        if (phaseModeActive_)
            return probePhaseThermoState(
                density, internalEnergy, overallMassFractions,
                temperatureGuess, temperature, pressure);

        if (equilibriumThermo_)
        {
            if (!injectedThermoIdentityMatches()) return false;
            try
            {
                const thermo::EquilibriumResult recovered =
                    equilibriumThermo_->recoverRhoE(
                        density, internalEnergy, overallMassFractions);
                if (!admittedInjectedResult(
                        recovered, overallMassFractions))
                    return false;
                temperature = recovered.state->T;
                pressure = recovered.state->p;
                return std::isfinite(temperature) && temperature > 0.0
                    && std::isfinite(pressure) && pressure > 0.0;
            }
            catch (...)
            {
                return false;
            }
        }

        const auto status = thermo_.updateRhoEY(
            density, internalEnergy, overallMassFractions,
            temperatureGuess);
        if (status != RealFluidState::Status::SUCCESS) return false;
        temperature = thermo_.T();
        pressure = thermo_.p();
        return true;
    }

    static int rhsFailure(
        RealFluidCVODEWrapper* self, bool permanent = false)
    {
        if (permanent || self->lastFailureClass_
                == RealFluidCVODEFailureClass::PERMANENT_ADMISSION)
        {
            self->permanentPhaseFailure();
            return -1;
        }
        self->lastFailureClass_ =
            RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
        return 1;
    }

    static int rhsFn(
        sunrealtype /*time*/, N_Vector y, N_Vector ydot, void* userData)
    {
        auto* self = static_cast<RealFluidCVODEWrapper*>(userData);
        if (self->equilibriumThermo_
            && !self->injectedThermoIdentityMatches())
            return rhsFailure(self, true);
        const sunrealtype* values = N_VGetArrayPointer(y);
        sunrealtype* derivatives = N_VGetArrayPointer(ydot);

        const std::vector<double>& molecularWeights =
            self->chemistry_->molecularWeights();
        for (int c = 0; c < self->nChem_; ++c)
            self->rhsChemY_[c] = values[c];
        if (!RhsCompositionPreparation::prepare(
                self->rho_, self->registry_, self->mapper_, molecularWeights,
                self->globalTemplate_, self->rhsChemY_, self->rhsGlobalY_,
                nullptr))
            return rhsFailure(self);

        const double specificEnergy = values[self->nChem_];
        if (!std::isfinite(specificEnergy)) return rhsFailure(self);

        double temperature = std::numeric_limits<double>::quiet_NaN();
        double pressure = std::numeric_limits<double>::quiet_NaN();
        if (!self->recoverThermoState(
                self->rho_, specificEnergy, self->rhsGlobalY_,
                self->temperatureSeed_, temperature, pressure))
            return rhsFailure(self);

        double heatRelease = 0.0;
        self->chemistry_->computeSourceTerms(
            temperature, pressure, self->rho_,
            self->rhsChemY_, self->rhsOmega_, heatRelease);

        if (!std::isfinite(heatRelease)) return rhsFailure(self);
        for (int c = 0; c < self->nChem_; ++c)
        {
            derivatives[c] = self->rhsOmega_[c] / self->rho_;
            if (!std::isfinite(derivatives[c])) return rhsFailure(self);
        }
        derivatives[self->nChem_] = heatRelease / self->rho_;
        return std::isfinite(derivatives[self->nChem_])
            ? 0 : rhsFailure(self);
    }

    void initializeIntegrator(
        double time,
        double rho,
        const std::vector<double>& globalY,
        double specificEnergy,
        double temperatureGuess)
    {
        rho_ = rho;
        globalTemplate_ = globalY;
        temperatureSeed_ = temperatureGuess;
        initializedTime_ = time;
        readyForAdvance_ = false;
        lastFailureClass_ = RealFluidCVODEFailureClass::NONE;

        if (initialized_)
        {
            sunrealtype* values = N_VGetArrayPointer(y_);
            for (int c = 0; c < nChem_; ++c)
                values[c] = globalY[mapper_.globalIndex(c)];
            values[nChem_] = specificEnergy;
            if (CVodeReInit(cvodeMem_, time, y_))
                throw std::runtime_error(
                    "RealFluidCVODEWrapper: CVodeReInit failed");
            sunrealtype* tolerances =
                N_VGetArrayPointer(absoluteTolerances_);
            for (int c = 0; c < nChem_; ++c)
                tolerances[c] = speciesAtol_;
            tolerances[nChem_] = energyAtol_;
            CVodeSVtolerances(cvodeMem_, rtol_, absoluteTolerances_);
            CVodeSetUserData(cvodeMem_, this);
            CVodeSetMaxNumSteps(cvodeMem_, 1000000);
            readyForAdvance_ = true;
            return;
        }

        if (SUNContext_Create(SUN_COMM_NULL, &sunctx_))
            throw std::runtime_error(
                "RealFluidCVODEWrapper: SUNContext_Create failed");
        y_ = N_VNew_Serial(nEq_, sunctx_);
        if (!y_)
            throw std::runtime_error(
                "RealFluidCVODEWrapper: N_Vector allocation failed");
        absoluteTolerances_ = N_VNew_Serial(nEq_, sunctx_);
        if (!absoluteTolerances_)
            throw std::runtime_error(
                "RealFluidCVODEWrapper: tolerance allocation failed");

        sunrealtype* values = N_VGetArrayPointer(y_);
        for (int c = 0; c < nChem_; ++c)
            values[c] = globalY[mapper_.globalIndex(c)];
        values[nChem_] = specificEnergy;
        sunrealtype* tolerances =
            N_VGetArrayPointer(absoluteTolerances_);
        for (int c = 0; c < nChem_; ++c)
            tolerances[c] = speciesAtol_;
        tolerances[nChem_] = energyAtol_;

        cvodeMem_ = CVodeCreate(CV_BDF, sunctx_);
        if (!cvodeMem_)
            throw std::runtime_error(
                "RealFluidCVODEWrapper: CVodeCreate failed");
        if (CVodeInit(cvodeMem_, &RealFluidCVODEWrapper::rhsFn, time, y_))
            throw std::runtime_error(
                "RealFluidCVODEWrapper: CVodeInit failed");

        CVodeSVtolerances(cvodeMem_, rtol_, absoluteTolerances_);
        CVodeSetUserData(cvodeMem_, this);
        matrix_ = SUNDenseMatrix(nEq_, nEq_, sunctx_);
        linearSolver_ = SUNLinSol_Dense(y_, matrix_, sunctx_);
        if (CVodeSetLinearSolver(cvodeMem_, linearSolver_, matrix_))
            throw std::runtime_error(
                "RealFluidCVODEWrapper: CVodeSetLinearSolver failed");
        CVodeSetMaxNumSteps(cvodeMem_, 1000000);
        initialized_ = true;
        readyForAdvance_ = true;
    }

    void freeAll()
    {
        if (cvodeMem_) { CVodeFree(&cvodeMem_); cvodeMem_ = nullptr; }
        if (linearSolver_)
        { SUNLinSolFree(linearSolver_); linearSolver_ = nullptr; }
        if (matrix_) { SUNMatDestroy(matrix_); matrix_ = nullptr; }
        if (absoluteTolerances_)
        { N_VDestroy(absoluteTolerances_); absoluteTolerances_ = nullptr; }
        if (y_) { N_VDestroy(y_); y_ = nullptr; }
        if (sunctx_) { SUNContext_Free(&sunctx_); sunctx_ = nullptr; }
        initialized_ = false;
        readyForAdvance_ = false;
        closePhaseTransaction();
    }

public:
    RealFluidCVODEWrapper(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        const SpeciesMapper& mapper)
    :
        chemistry_(chemistry), registry_(registry), mapper_(mapper),
        nChem_(mapper.nChem()), nEq_(mapper.nChem() + 1),
        rho_(0.0), rtol_(1.0e-6), speciesAtol_(1.0e-11),
        energyAtol_(1.0e-2),
        temperatureSeed_(300.0),
        globalTemplate_(registry.size(), 0.0),
        rhsChemY_(nChem_, 0.0),
        rhsGlobalY_(registry.size(), 0.0),
        rhsOmega_(nChem_, 0.0),
        equilibriumThermo_(nullptr),
        sunctx_(nullptr), cvodeMem_(nullptr), y_(nullptr),
        absoluteTolerances_(nullptr), matrix_(nullptr),
        linearSolver_(nullptr), initialized_(false)
    {
        thermo_.setSpecies(registry_.names());
    }

    RealFluidCVODEWrapper(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        const SpeciesMapper& mapper,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel)
        : RealFluidCVODEWrapper(chemistry, registry, mapper)
    {
        bindInjectedThermoModel(std::move(thermoModel));
    }

    RealFluidCVODEWrapper(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        const SpeciesMapper& mapper,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel,
        const phaseEquilibrium::PRPhaseRhoERecovery& phaseRecovery)
        : RealFluidCVODEWrapper(
            chemistry, registry, mapper, std::move(thermoModel))
    {
        phaseThermoOwner_ = std::make_unique<
            phaseEquilibrium::PRPhaseThermoOwnerAdapter>(phaseRecovery);
    }

    RealFluidCVODEWrapper(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        const SpeciesMapper& mapper,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel,
        const phaseEquilibrium::PRPhaseRhoERecovery& phaseRecovery,
        const phaseEquilibrium::PRPhaseDisplacementAuthority&
            phaseDisplacement)
        : RealFluidCVODEWrapper(
            chemistry, registry, mapper, std::move(thermoModel))
    {
        phaseThermoOwner_ = std::make_unique<
            phaseEquilibrium::PRPhaseThermoOwnerAdapter>(
                phaseRecovery, phaseDisplacement, *equilibriumThermo_);
    }

    ~RealFluidCVODEWrapper() { freeAll(); }

    bool hasInjectedThermoModel() const
    {
        return static_cast<bool>(equilibriumThermo_);
    }

    bool injectedThermoIdentityValid() const
    {
        return equilibriumThermo_ && injectedThermoIdentityMatches();
    }

    bool hasPhaseThermoOwner() const noexcept
    {
        return static_cast<bool>(phaseThermoOwner_);
    }

    bool phaseTransactionActive() const noexcept
    {
        return phaseModeActive_ && phaseThermoOwner_
            && phaseThermoOwner_->transactionActive();
    }

    RealFluidCVODEFailureClass lastFailureClass() const noexcept
    {
        return lastFailureClass_;
    }

    phaseEquilibrium::PRPhaseThermoOwnerStatus
    lastPhaseOwnerStatus() const noexcept
    {
        return lastPhaseOwnerStatus_;
    }

    phaseEquilibrium::PRPhaseDisplacementStatus
    lastPhaseDisplacementStatus() const noexcept
    {
        return lastPhaseDisplacementStatus_;
    }

    const std::optional<phaseEquilibrium::PRPhaseDisplacementDiagnostics>&
    lastPhaseDisplacementDiagnostics() const noexcept
    {
        return lastPhaseDisplacementDiagnostics_;
    }

    std::uint64_t phaseProbeCount() const noexcept
    {
        return phaseProbeCount_;
    }

    std::uint64_t phaseAcceptCount() const noexcept
    {
        return phaseAcceptCount_;
    }

    bool lastFailureIsPermanent() const noexcept
    {
        return lastFailureClass_
            == RealFluidCVODEFailureClass::PERMANENT_ADMISSION;
    }

    void abortPhase() noexcept
    {
        readyForAdvance_ = false;
        closePhaseTransaction();
    }

    const std::string& boundThermoModelId() const
    {
        return boundThermoModelId_;
    }

    const std::vector<std::string>& boundThermoSpeciesNames() const
    {
        return boundThermoSpeciesNames_;
    }

    EnergyConvention boundEnergyConvention() const
    {
        return boundEnergyConvention_;
    }

    bool recoverAdmittedRhoEY(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        double temperatureGuess,
        double& temperature,
        double& pressure)
    {
        return recoverThermoState(
            density, internalEnergy, overallMassFractions,
            temperatureGuess, temperature, pressure);
    }

    bool initializeAdmittedPTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions,
        double& density,
        double& internalEnergy,
        double& heatCapacityCv)
    {
        if (equilibriumThermo_)
        {
            if (!injectedThermoIdentityMatches()) return false;
            try
            {
                const thermo::EquilibriumResult initialized =
                    equilibriumThermo_->initializePTY(
                        pressure, temperature, overallMassFractions);
                if (!admittedInjectedResult(
                        initialized, overallMassFractions)
                    || !initialized.state->thermodynamics.has_value())
                    return false;
                density = initialized.state->rho;
                internalEnergy = initialized.state->e;
                heatCapacityCv = initialized.state->thermodynamics->Cv;
                return std::isfinite(density) && density > 0.0
                    && std::isfinite(internalEnergy)
                    && std::isfinite(heatCapacityCv)
                    && heatCapacityCv > 0.0;
            }
            catch (...)
            {
                return false;
            }
        }

        if (thermo_.updatePTY(
                pressure, temperature, overallMassFractions)
            != RealFluidState::Status::SUCCESS)
            return false;
        density = thermo_.rho();
        internalEnergy = thermo_.e();
        heatCapacityCv = thermo_.Cv();
        return true;
    }

    void setTolerances(double relative, double speciesAbsolute)
    {
        if (!std::isfinite(relative) || relative <= 0.0 ||
            !std::isfinite(speciesAbsolute) || speciesAbsolute <= 0.0)
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: invalid tolerances");
        rtol_ = relative;
        speciesAtol_ = speciesAbsolute;
    }

    void setTolerances(
        double relative,
        double speciesAbsolute,
        double energyAbsolute)
    {
        setTolerances(relative,speciesAbsolute);
        if (!std::isfinite(energyAbsolute) || energyAbsolute <= 0.0)
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: invalid energy tolerance");
        energyAtol_ = energyAbsolute;
    }

    bool initPhase(
        double time,
        const phaseEquilibrium::PRPhaseThermoParentView& parent,
        double temperatureSeed,
        RealFluidCVODEPhaseState& initialPhase)
    {
        readyForAdvance_ = false;
        lastFailureClass_ = RealFluidCVODEFailureClass::NONE;
        lastPhaseOwnerStatus_ =
            phaseEquilibrium::PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION;
        lastPhaseDisplacementStatus_ = phaseEquilibrium::
            PRPhaseDisplacementStatus::INVALID_CONFIGURATION;
        lastPhaseDisplacementDiagnostics_.reset();
        phaseProbeCount_ = 0U;
        phaseAcceptCount_ = 0U;
        closePhaseTransaction();
        if (!phaseThermoOwner_ || !equilibriumThermo_
            || !injectedThermoIdentityMatches()
            || parent.model != equilibriumThermo_.get()
            || !std::isfinite(time)
            || !std::isfinite(temperatureSeed) || temperatureSeed <= 0.0)
            return permanentPhaseFailure();

        const phaseEquilibrium::PRPhaseThermoOwnerResult opened =
            phaseThermoOwner_->open(parent, temperatureSeed);
        const thermo::EquilibriumPhaseProperties* parentPhase = nullptr;
        if (parent.committed)
        {
            if (parent.admission.phaseIdentity
                    == phaseEquilibrium::PRPhaseIdentity::VAPOR
                && parent.committed->Q.vapor)
                parentPhase = &*parent.committed->Q.vapor;
            else if (parent.admission.phaseIdentity
                    == phaseEquilibrium::PRPhaseIdentity::SUPERCRITICAL
                && parent.committed->Q.singlePhase)
                parentPhase = &*parent.committed->Q.singlePhase;
        }
        if (!opened.success() || !opened.phaseState
            || !opened.acceptedToken
            || !parentPhase
            || parentPhase->massFractions.size()
                != static_cast<std::size_t>(registry_.size()))
            return permanentPhaseFailure();

        phaseParent_ = parent;
        phaseToken_ = *opened.acceptedToken;
        RealFluidCVODEPhaseState authoritative;
        authoritative.rho = parentPhase->rho;
        authoritative.specificEnergy = parentPhase->e;
        authoritative.massFractions = parentPhase->massFractions;
        authoritative.temperature = parent.committed->Q.T;
        authoritative.pressure = parent.committed->Q.p;
        phaseInitialState_ = authoritative;
        phaseModeActive_ = true;
        try
        {
            initializeIntegrator(
                time, authoritative.rho,
                authoritative.massFractions,
                authoritative.specificEnergy,
                authoritative.temperature);
        }
        catch (...)
        {
            permanentPhaseFailure();
            throw;
        }
        initialPhase = std::move(authoritative);
        return true;
    }

    void init(
        double time, double rho,
        const std::vector<double>& globalY,
        double specificEnergy, double temperatureGuess)
    {
        readyForAdvance_ = false;
        lastFailureClass_ = RealFluidCVODEFailureClass::NONE;
        closePhaseTransaction();
        if (!std::isfinite(rho) || rho <= 0.0 ||
            globalY.size() != static_cast<size_t>(registry_.size()) ||
            !std::isfinite(specificEnergy) ||
            !std::isfinite(temperatureGuess) || temperatureGuess <= 0.0)
            throw std::invalid_argument(
                "RealFluidCVODEWrapper: invalid initial state");

        if (equilibriumThermo_)
        {
            double admittedTemperature =
                std::numeric_limits<double>::quiet_NaN();
            double admittedPressure =
                std::numeric_limits<double>::quiet_NaN();
            if (!recoverThermoState(
                    rho, specificEnergy, globalY, temperatureGuess,
                    admittedTemperature, admittedPressure))
                throw std::invalid_argument(
                    "RealFluidCVODEWrapper: initial thermo state is not "
                    "an admitted single-vapor or supercritical state");
            temperatureGuess = admittedTemperature;
        }
        initializeIntegrator(
            time, rho, globalY, specificEnergy, temperatureGuess);
    }

    bool advanceTo(
        double outputTime,
        std::vector<double>& globalY,
        double& specificEnergy,
        double& temperature,
        double& pressure)
    {
        if (!initialized_ || !readyForAdvance_)
            throw std::runtime_error(
                "RealFluidCVODEWrapper: advanceTo before init");
        if (equilibriumThermo_ && !injectedThermoIdentityMatches())
        {
            lastFailureClass_ =
                RealFluidCVODEFailureClass::PERMANENT_ADMISSION;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }
        if (!std::isfinite(outputTime) || outputTime < initializedTime_)
        {
            lastFailureClass_ =
                RealFluidCVODEFailureClass::PERMANENT_ADMISSION;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }

        if (phaseModeActive_ && outputTime == initializedTime_)
        {
            if (!phaseInitialState_)
                return permanentPhaseFailure();
            double validatedTemperature =
                std::numeric_limits<double>::quiet_NaN();
            double validatedPressure =
                std::numeric_limits<double>::quiet_NaN();
            if (!acceptPhaseThermoState(
                    phaseInitialState_->rho,
                    phaseInitialState_->specificEnergy,
                    phaseInitialState_->massFractions,
                    phaseInitialState_->temperature,
                    validatedTemperature,
                    validatedPressure))
                return false;
            std::vector<double> candidateY =
                phaseInitialState_->massFractions;
            const double candidateEnergy =
                phaseInitialState_->specificEnergy;
            const double candidateTemperature = phaseInitialState_->temperature;
            const double candidatePressure = phaseInitialState_->pressure;
            closePhaseTransaction();
            readyForAdvance_ = false;
            lastFailureClass_ = RealFluidCVODEFailureClass::NONE;
            globalY = std::move(candidateY);
            specificEnergy = candidateEnergy;
            temperature = candidateTemperature;
            pressure = candidatePressure;
            return true;
        }

        sunrealtype reachedTime = 0.0;
        const int flag = CVode(
            cvodeMem_, outputTime, y_, &reachedTime, CV_NORMAL);
        const double timeTolerance = 16.0*std::numeric_limits<double>::epsilon()
            *std::max(std::abs(outputTime),1.0e-30);
        if (flag != CV_SUCCESS ||
            std::abs(static_cast<double>(reachedTime)-outputTime)>timeTolerance)
        {
            if (lastFailureClass_ == RealFluidCVODEFailureClass::NONE)
                lastFailureClass_ =
                    RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }

        const sunrealtype* values = N_VGetArrayPointer(y_);
        std::vector<double> candidateY = globalTemplate_;
        for (int c = 0; c < nChem_; ++c)
            candidateY[mapper_.globalIndex(c)] = values[c];
        const double candidateEnergy = values[nChem_];
        double candidateTemperature =
            std::numeric_limits<double>::quiet_NaN();
        double candidatePressure =
            std::numeric_limits<double>::quiet_NaN();
        const bool recovered = std::isfinite(candidateEnergy)
            && (phaseModeActive_
                ? acceptPhaseThermoState(
                    rho_, candidateEnergy, candidateY, temperatureSeed_,
                    candidateTemperature, candidatePressure)
                : recoverThermoState(
                    rho_, candidateEnergy, candidateY, temperatureSeed_,
                    candidateTemperature, candidatePressure));
        if (!recovered)
        {
            if (lastFailureClass_ == RealFluidCVODEFailureClass::NONE)
                lastFailureClass_ =
                    RealFluidCVODEFailureClass::RETRYABLE_NUMERICAL;
            readyForAdvance_ = false;
            closePhaseTransaction();
            return false;
        }
        if (phaseModeActive_)
        {
            closePhaseTransaction();
            readyForAdvance_ = false;
        }
        lastFailureClass_ = RealFluidCVODEFailureClass::NONE;
        globalY = std::move(candidateY);
        specificEnergy = candidateEnergy;
        temperature = candidateTemperature;
        pressure = candidatePressure;
        return true;
    }

    long numSteps() const
    {
        long count = 0;
        if (cvodeMem_) CVodeGetNumSteps(cvodeMem_, &count);
        return count;
    }
};

} // End namespace Foam

#endif
