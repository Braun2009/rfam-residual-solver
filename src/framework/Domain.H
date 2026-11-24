#ifndef Domain_H
#define Domain_H

/*---------------------------------------------------------------------------*\
    Domain - unified mesh + state container (Stage 0/4 framework)

    Supports 1D and 2D structured grids with a SINGLE linear storage layout.
    All existing 1D API is preserved UNCHANGED so the MVP/precond solvers
    remain byte-for-byte identical.

    1D (Ny==1): U(i) / W(i) / N() / dx() / L() / x(i)  — legacy path.
    2D (Ny> 1): U(i,j) / W(i,j) / Nx() / Ny() / dx() / dy() / Ly()
                x(i,j) / y(i,j).

    Linear storage: index = j * Nx_ + i.  The 1D path just uses Ny_=1, Ly_=0.

    Additive-update workflow (unchanged):
        beginStep(): Unew <- U
        operators:   read U/W (t^n), write/accumulate into Unew
        commit():    NaN/neg check; U <- Unew; refresh primitives W
\*---------------------------------------------------------------------------*/

#include "CompleteCellState.H"
#include "StateContract.H"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Foam
{

class SectorRestartOwnerCodec;

class Domain
{
public:
    enum class FailureCode
    {
        NONE,
        INVALID_CONSERVED_STATE,
        STATE_RECOVERY_FAILED,
        PRIMITIVE_CONVERSION_FAILED,
        OPERATOR_EXCEPTION,
        TRANSACTION_PROTOCOL
    };

    struct FailureInfo
    {
        FailureCode code = FailureCode::NONE;
        int cell = -1;
        std::string component;
        std::string detail;
        thermo::EquilibriumStatus equilibriumStatus =
            thermo::EquilibriumStatus::SUCCESS;
        thermo::EquilibriumDiagnostics thermoDiagnostics;
    };

    struct CompleteCellReplacement
    {
        int i = -1;
        int j = -1;
        int k = -1;
        PrimitiveVars primitive;
        thermo::EquilibriumState equilibrium;
    };

    static const char* failureCodeName(FailureCode code)
    {
        switch (code)
        {
            case FailureCode::NONE: return "NONE";
            case FailureCode::INVALID_CONSERVED_STATE:
                return "INVALID_CONSERVED_STATE";
            case FailureCode::STATE_RECOVERY_FAILED:
                return "STATE_RECOVERY_FAILED";
            case FailureCode::PRIMITIVE_CONVERSION_FAILED:
                return "PRIMITIVE_CONVERSION_FAILED";
            case FailureCode::OPERATOR_EXCEPTION: return "OPERATOR_EXCEPTION";
            case FailureCode::TRANSACTION_PROTOCOL:
                return "TRANSACTION_PROTOCOL";
        }
        return "UNKNOWN";
    }

private:
    friend class SectorRestartOwnerCodec;

    int    Nx_, Ny_, Nz_, nSpecies_;
    double Lx_, Ly_, Lz_, dx_, dy_, dz_;
    std::vector<CompleteCellState> cells_;
    std::vector<ConservedVars> Unew_;
    std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel_;
    bool externalThermoModel_ = false;
    bool stateReady_ = false;
    bool   tripped_;
    FailureInfo failure_;
    ConservativeStateContract stateContract_;

    // A complete solver step is a transaction over both the authoritative
    // conserved field and its derived primitive cache.  RK stages may commit
    // temporary states, but a rejected stage restores these exact snapshots.
    std::vector<CompleteCellState> cellsStepStart_;
    std::vector<ConservedVars> UnewStepStart_;
    bool stepTransactionActive_ = false;
    std::uint64_t acceptedSteps_ = 0;
    std::uint64_t rejectedSteps_ = 0;

    // RK2 temporary storage: holds U^n during predictor-corrector
    mutable std::vector<ConservedVars> Usaved_;

    // Explicit initialization/import staging. Committed W and Q never expose a
    // mutable reference, so a failed bulk conversion cannot desynchronize them.
    std::vector<PrimitiveVars> primitiveStage_;
    bool primitiveStageActive_ = true;

    // Linear index: k*Nx*Ny + j*Nx + i.  1D: Ny_=Nz_=1 → index=i.
    int idx(int i, int j) const { return j * Nx_ + i; }
    int idx(int i, int j, int k) const { return k * Nx_ * Ny_ + j * Nx_ + i; }

    void setFailure(
        FailureCode code,int cell,
        const std::string& component,const std::string& detail,
        thermo::EquilibriumStatus equilibriumStatus =
            thermo::EquilibriumStatus::SUCCESS,
        const thermo::EquilibriumDiagnostics& thermoDiagnostics = {})
    {
        // Boundary/state conversion can report failures from an OpenMP loop.
        // Keep the first failure deterministic and avoid racing std::string.
        #pragma omp critical(rfam_domain_failure)
        {
            tripped_ = true;
            if (failure_.code == FailureCode::NONE)
                failure_ = FailureInfo{
                    code,cell,component,detail,
                    equilibriumStatus,thermoDiagnostics};
        }
    }

    void configureThermoModel(
        std::shared_ptr<const thermo::EquilibriumThermoModel> model)
    {
        externalThermoModel_ = static_cast<bool>(model);
        thermoModel_ = model ? std::move(model)
            : thermo::makeLegacySinglePhaseThermoModel();
        if (!thermoModel_)
            throw std::invalid_argument("Domain requires an equilibrium model");
        if (externalThermoModel_ && thermoModel_->speciesNames().empty())
            throw std::invalid_argument(
                "an external equilibrium model must declare its species layout");
        if (thermoModel_->energyConvention()
            != stateContract_.energyConvention())
            throw std::invalid_argument(
                "Domain equilibrium model energy convention mismatch");
        const auto& species = thermoModel_->speciesNames();
        if (!species.empty())
            stateContract_.setSpeciesLayout(
                species,static_cast<std::size_t>(nSpecies_));
    }

    bool recoverConservedField(
        const std::vector<ConservedVars>& source,
        const std::vector<PrimitiveVars>& seeds,
        std::vector<CompleteCellState>& target,
        int* failedCell = nullptr,
        thermo::EquilibriumStatus* failedStatus = nullptr,
        thermo::EquilibriumDiagnostics* failedDiagnostics = nullptr,
        const std::vector<CompleteCellState>* previousStates = nullptr,
        const std::vector<thermo::EquilibriumSeed>* equilibriumSeeds = nullptr) const
    {
        const int total = Nx_ * Ny_ * Nz_;
        if (source.size() != static_cast<std::size_t>(total)
            || seeds.size() != static_cast<std::size_t>(total))
            return false;
        target.assign(total,CompleteCellState(nSpecies_));
        std::vector<thermo::EquilibriumStatus> statuses(
            total,thermo::EquilibriumStatus::INVALID_CONFIGURATION);
        std::vector<thermo::EquilibriumDiagnostics> diagnostics(total);
        int firstInvalid = total;
        #pragma omp parallel for reduction(min:firstInvalid) schedule(static)
        for (int k = 0; k < total; ++k)
        {
#ifdef RFAM_THERMO_DIAGNOSTICS
            if (k == 1448 || k == 1449)
                std::cerr << "[RCF-DIAG] k=" << k
                          << " rho=" << source[k].rho
                          << " rhoE=" << source[k].rhoE << "\n";
#endif
            thermo::EquilibriumSeed phaseSeed;
            const thermo::EquilibriumSeed* phaseSeedPtr = nullptr;
            if (equilibriumSeeds
                && equilibriumSeeds->size() == source.size())
            {
                phaseSeed = (*equilibriumSeeds)[k];
                phaseSeedPtr = &phaseSeed;
            }
            else if (previousStates
                && previousStates->size() == source.size())
            {
                const thermo::EquilibriumState& previousQ =
                    (*previousStates)[k].Q;
                if (previousQ.region == thermo::EquilibriumRegion::TWO_PHASE
                    && std::isfinite(previousQ.T)
                    && std::isfinite(previousQ.vaporMassFraction))
                {
                    phaseSeed.previousRegion = previousQ.region;
                    phaseSeed.previousTemperature = previousQ.T;
                    phaseSeed.previousVaporMassFraction =
                        previousQ.vaporMassFraction;
                    phaseSeedPtr = &phaseSeed;
                }
            }
            CompleteCellConversionResult converted =
                CompleteStateConversion::fromConserved(
                    source[k],seeds[k],*thermoModel_,phaseSeedPtr);
            statuses[k] = converted.status;
            diagnostics[k] = std::move(converted.diagnostics);
            if (!converted.success()) firstInvalid = std::min(firstInvalid,k);
            else target[k] = std::move(*converted.state);
        }
        if (failedCell) *failedCell = firstInvalid == total ? -1 : firstInvalid;
        if (firstInvalid != total)
        {
            if (failedStatus) *failedStatus = statuses[firstInvalid];
            if (failedDiagnostics)
                *failedDiagnostics = std::move(diagnostics[firstInvalid]);
        }
        return firstInvalid == total;
    }

    // Used only while constructing a detached restart Domain. Live restart
    // publication remains owned by the generation coordinator.
    void restoreDetachedRestartStepLedger(
        std::uint64_t acceptedSteps,
        std::uint64_t rejectedSteps) noexcept
    {
        acceptedSteps_ = acceptedSteps;
        rejectedSteps_ = rejectedSteps;
    }

    bool recoverPrimitiveField(
        const std::vector<PrimitiveVars>& source,
        std::vector<CompleteCellState>& target,
        int* failedCell = nullptr,
        thermo::EquilibriumStatus* failedStatus = nullptr,
        thermo::EquilibriumDiagnostics* failedDiagnostics = nullptr) const
    {
        const int total = Nx_ * Ny_ * Nz_;
        if (source.size() != static_cast<std::size_t>(total)) return false;
        target.assign(total,CompleteCellState(nSpecies_));
        std::vector<thermo::EquilibriumStatus> statuses(
            total,thermo::EquilibriumStatus::INVALID_CONFIGURATION);
        std::vector<thermo::EquilibriumDiagnostics> diagnostics(total);
        int firstInvalid = total;
        #pragma omp parallel for reduction(min:firstInvalid) schedule(static)
        for (int k = 0; k < total; ++k)
        {
            CompleteCellConversionResult converted =
                CompleteStateConversion::fromPrimitive(
                    source[k],*thermoModel_);
            statuses[k] = converted.status;
            diagnostics[k] = std::move(converted.diagnostics);
            if (!converted.success()) firstInvalid = std::min(firstInvalid,k);
            else target[k] = std::move(*converted.state);
        }
        if (failedCell) *failedCell = firstInvalid == total ? -1 : firstInvalid;
        if (firstInvalid != total)
        {
            if (failedStatus) *failedStatus = statuses[firstInvalid];
            if (failedDiagnostics)
                *failedDiagnostics = std::move(diagnostics[firstInvalid]);
        }
        return firstInvalid == total;
    }

    void publishCompleteField(std::vector<CompleteCellState>& candidate)
    {
        cells_.swap(candidate);
        Unew_.resize(cells_.size());
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < N(); ++k) Unew_[k] = cells_[k].U;
        stateReady_ = true;
        primitiveStageActive_ = false;
    }

    bool replaceConservedFieldImpl(
        const std::vector<ConservedVars>& source,
        const std::vector<double>& temperatureSeeds,
        bool requireActiveTransaction)
    {
        if (tripped_) return false;
        const char* component = requireActiveTransaction
            ? "Domain::replaceConservedFieldInStepTransaction"
            : "Domain::replaceConservedField";
        if (stepTransactionActive_ != requireActiveTransaction
            || (requireActiveTransaction && !stateReady_))
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,component,
                requireActiveTransaction
                    ? "inner conserved replacement requires an active transaction"
                    : "conserved replacement requires no active transaction");
            return false;
        }
        if (source.size() != cells_.size()
            || (!temperatureSeeds.empty()
                && temperatureSeeds.size() != source.size()))
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,component,
                "replacement field dimensions do not match Domain");
            return false;
        }
        std::vector<PrimitiveVars> seeds(
            source.size(),PrimitiveVars(nSpecies_));
        for (std::size_t k = 0; k < source.size(); ++k)
        {
            if (stateReady_) seeds[k] = cells_[k].W;
            if (!temperatureSeeds.empty()
                && std::isfinite(temperatureSeeds[k]))
                seeds[k].T = temperatureSeeds[k];
        }
        std::vector<CompleteCellState> candidate;
        int failedCell = -1;
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics diagnostics;
        if (!recoverConservedField(
                source,seeds,candidate,&failedCell,
                &failedStatus,&diagnostics,&cells_))
        {
            setFailure(FailureCode::STATE_RECOVERY_FAILED,failedCell,
                component,"replacement conserved field cannot be recovered",
                failedStatus,diagnostics);
            return false;
        }
        publishCompleteField(candidate);
        return true;
    }

    void preparePrimitiveStage()
    {
        if (stepTransactionActive_)
            throw std::logic_error(
                "cannot stage primitive state during a step transaction");
        if (primitiveStageActive_) return;
        primitiveStage_.resize(cells_.size(),PrimitiveVars(nSpecies_));
        for (std::size_t k = 0; k < cells_.size(); ++k)
            primitiveStage_[k] = cells_[k].W;
        primitiveStageActive_ = true;
    }

    bool setCompleteCell(int id,const PrimitiveVars& primitive)
    {
        if (!stateReady_ || stepTransactionActive_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,id,
                "Domain::setCell",
                "cell update requires a ready state outside a transaction");
            return false;
        }
        CompleteCellConversionResult converted =
            CompleteStateConversion::fromPrimitive(primitive,*thermoModel_);
        if (!converted.success())
        {
            setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,id,
                "Domain::setCell",
                "primitive state is not a complete admissible state",
                converted.status,converted.diagnostics);
            return false;
        }
        cells_[id] = std::move(*converted.state);
        Unew_[id] = cells_[id].U;
        return true;
    }

    bool setCompleteCell(
        int id,const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium)
    {
        if (!stateReady_ || stepTransactionActive_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,id,
                "Domain::setCell",
                "cell update requires a ready state outside a transaction");
            return false;
        }
        CompleteCellConversionResult converted =
            CompleteStateConversion::fromPrimitiveAndEquilibrium(
                primitive,equilibrium,*thermoModel_);
        if (!converted.success())
        {
            setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,id,
                "Domain::setCell",
                "primitive/equilibrium state pair is inconsistent",
                converted.status,converted.diagnostics);
            return false;
        }
        cells_[id] = std::move(*converted.state);
        Unew_[id] = cells_[id].U;
        return true;
    }

public:
    // ---- 1D constructor (legacy, unchanged API) ----
    Domain(
        int N,double L,int nSpecies,
        std::shared_ptr<const thermo::EquilibriumThermoModel> model = {})
        : Nx_(N), Ny_(1), Nz_(1), nSpecies_(nSpecies),
          Lx_(L), Ly_(0.0), Lz_(0.0), dx_(L/N), dy_(0.0), dz_(0.0),
          cells_(N,CompleteCellState(nSpecies)),
          Unew_(N,ConservedVars(nSpecies)),tripped_(false),
          primitiveStage_(N,PrimitiveVars(nSpecies))
    { configureThermoModel(std::move(model)); }

    // ---- 2D constructor ----
    Domain(
        int Nx,int Ny,double Lx,double Ly,int nSpecies,
        std::shared_ptr<const thermo::EquilibriumThermoModel> model = {})
        : Nx_(Nx), Ny_(Ny), Nz_(1), nSpecies_(nSpecies),
          Lx_(Lx), Ly_(Ly), Lz_(0.0), dx_(Lx/Nx), dy_(Ly/Ny), dz_(0.0),
          cells_(Nx*Ny,CompleteCellState(nSpecies)),
          Unew_(Nx*Ny,ConservedVars(nSpecies)),tripped_(false),
          primitiveStage_(Nx*Ny,PrimitiveVars(nSpecies))
    { configureThermoModel(std::move(model)); }

    // ---- 3D constructor ----
    Domain(
        int Nx,int Ny,int Nz,double Lx,double Ly,double Lz,int nSpecies,
        std::shared_ptr<const thermo::EquilibriumThermoModel> model = {})
        : Nx_(Nx), Ny_(Ny), Nz_(Nz), nSpecies_(nSpecies),
          Lx_(Lx), Ly_(Ly), Lz_(Lz), dx_(Lx/Nx), dy_(Ly/Ny), dz_(Lz/Nz),
          cells_(Nx*Ny*Nz,CompleteCellState(nSpecies)),
          Unew_(Nx*Ny*Nz,ConservedVars(nSpecies)),tripped_(false),
          primitiveStage_(Nx*Ny*Nz,PrimitiveVars(nSpecies))
    { configureThermoModel(std::move(model)); }

    // ---- Grid metadata ----
    // 1D legacy accessors (unchanged):
    int    N()  const { return Nx_ * Ny_ * Nz_; }   // total cell count
    double L()  const { return Lx_; }
    double dx() const { return dx_; }
    double x(int i) const { return (i + 0.5) * dx_; }

    // 2D accessors:
    int    Nx() const { return Nx_; }
    int    Ny() const { return Ny_; }
    double Lx() const { return Lx_; }
    double Ly() const { return Ly_; }

    // Public linear-index mapping. Residual/rate fields are stored with the
    // same layout as the state (k*Nx*Ny + j*Nx + i), so an operator can index
    // its own per-cell rate without owning a copy of the grid metadata.
    int    index(int i) const { return i; }
    int    index(int i,int j) const { return idx(i,j); }
    int    index(int i,int j,int k) const { return idx(i,j,k); }

    // 3D accessors:
    int    Nz() const { return Nz_; }
    double Lz() const { return Lz_; }
    double dz() const { return dz_; }
    double z(int /*i*/, int /*j*/, int k) const { return (k + 0.5) * dz_; }

    // 3D state access:
    const ConservedVars& U(int i,int j,int k) const
    { return cells_[idx(i,j,k)].U; }
    const PrimitiveVars& W(int i,int j,int k) const
    { return cells_[idx(i,j,k)].W; }
    const thermo::EquilibriumState& Q(int i,int j,int k) const
    { return cells_[idx(i,j,k)].Q; }
    PrimitiveVars& stagedW(int i,int j,int k)
    { preparePrimitiveStage(); return primitiveStage_[idx(i,j,k)]; }
    ConservedVars&       Unew(int i, int j, int k)       { return Unew_[idx(i,j,k)]; }
    const ConservedVars& Unew(int i, int j, int k) const { return Unew_[idx(i,j,k)]; }
    double dy() const { return dy_; }
    double x(int i, int /*j*/) const { return (i + 0.5) * dx_; }
    double y(int /*i*/, int j)  const { return (j + 0.5) * dy_; }

    bool   tripped()  const { return tripped_; }
    bool   stateReady() const { return stateReady_; }
    int    nSpecies() const { return nSpecies_; }
    const FailureInfo& failure() const { return failure_; }
    std::uint64_t acceptedSteps() const { return acceptedSteps_; }
    std::uint64_t rejectedSteps() const { return rejectedSteps_; }

    // Clear the trip flag and failure record. Intended ONLY for integrators
    // that roll back to a consistent step-start snapshot (transaction) and
    // retry with a reduced step; the caller is responsible for restoring the
    // conserved field before resuming operations.
    void resetTrip()
    {
        tripped_ = false;
        failure_ = FailureInfo{};
    }
    EnergyConvention energyConvention() const
    { return stateContract_.energyConvention(); }
    const thermo::EquilibriumThermoModel& thermoModel() const
    { return *thermoModel_; }
    std::shared_ptr<const thermo::EquilibriumThermoModel>
    thermoModelOwner() const
    { return thermoModel_; }
    const std::string& thermoModelId() const { return thermoModel_->modelId(); }

    void recordFailure(
        FailureCode code,int cell,
        const std::string& component,const std::string& detail)
    { setFailure(code,cell,component,detail); }

    // ---- Species names (EOS reback opt-in) ----
    void setSpecies(const std::vector<std::string>& names)
    {
        if (stateReady_ || stepTransactionActive_)
        {
            if (stateContract_.speciesNames() == names) return;
            throw std::logic_error(
                "cannot change Domain species after state publication");
        }
        if (externalThermoModel_
            && thermoModel_->speciesNames() != names)
            throw std::invalid_argument(
                "Domain species do not match the equilibrium model");
        stateContract_.setSpeciesLayout(
            names,static_cast<std::size_t>(nSpecies_));
        if (!externalThermoModel_)
            thermoModel_ = thermo::makeLegacySinglePhaseThermoModel(names);
    }
    const std::vector<std::string>& speciesNames() const
    { return stateContract_.speciesNames(); }

    // ---- State access — 1D legacy (read-only) ----
    const ConservedVars& U(int i) const { return cells_[i].U; }
    const PrimitiveVars& W(int i) const { return cells_[i].W; }
    const thermo::EquilibriumState& Q(int i) const { return cells_[i].Q; }
    PrimitiveVars& stagedW(int i)
    { preparePrimitiveStage(); return primitiveStage_[i]; }
    ConservedVars&       Unew(int i)       { return Unew_[i]; }
    const ConservedVars& Unew(int i) const { return Unew_[i]; }

    // ---- State access — 2D ----
    const ConservedVars& U(int i,int j) const
    { return cells_[idx(i,j)].U; }
    const PrimitiveVars& W(int i,int j) const
    { return cells_[idx(i,j)].W; }
    const thermo::EquilibriumState& Q(int i,int j) const
    { return cells_[idx(i,j)].Q; }
    PrimitiveVars& stagedW(int i,int j)
    { preparePrimitiveStage(); return primitiveStage_[idx(i,j)]; }
    ConservedVars&       Unew(int i, int j)       { return Unew_[idx(i,j)]; }
    const ConservedVars& Unew(int i, int j) const { return Unew_[idx(i,j)]; }

    // Set a boundary cell: writes W and updates U (and Unew) from P2C.
    // Call before beginStep() so the updated state is the t^n base.
    bool setCell(int i,int j,const PrimitiveVars& w)
    { return setCompleteCell(idx(i,j),w); }
    bool setCell(
        int i,int j,const PrimitiveVars& w,
        const thermo::EquilibriumState& q)
    { return setCompleteCell(idx(i,j),w,q); }

    // 3D version
    bool setCell(int i,int j,int k,const PrimitiveVars& w)
    { return setCompleteCell(idx(i,j,k),w); }
    bool setCell(
        int i,int j,int k,const PrimitiveVars& w,
        const thermo::EquilibriumState& q)
    { return setCompleteCell(idx(i,j,k),w,q); }

    // Validate every sparse replacement before publishing any of them. This
    // is the boundary-condition counterpart of replaceCompleteField(): a bad
    // nonfirst cell may trip the Domain, but it cannot leave a partially
    // updated authoritative or candidate field behind.
    bool replaceCompleteCellsImpl(
        const std::vector<CompleteCellReplacement>& replacements,
        bool requireActiveTransaction)
    {
        if (tripped_) return false;
        if (!stateReady_
            || stepTransactionActive_ != requireActiveTransaction)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::replaceCompleteCells",
                requireActiveTransaction
                    ? "sparse inner replacement requires an active transaction"
                    : "sparse replacement requires a ready state outside a transaction");
            return false;
        }

        std::vector<int> cellIds;
        cellIds.reserve(replacements.size());
        std::vector<unsigned char> seen(
            static_cast<std::size_t>(N()),0);
        for (const CompleteCellReplacement& replacement : replacements)
        {
            if (replacement.i < 0 || replacement.i >= Nx_
                || replacement.j < 0 || replacement.j >= Ny_
                || replacement.k < 0 || replacement.k >= Nz_)
            {
                setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                    "Domain::replaceCompleteCells",
                    "sparse replacement cell is outside the Domain");
                return false;
            }
            const int cell = idx(
                replacement.i,replacement.j,replacement.k);
            if (seen[static_cast<std::size_t>(cell)] != 0)
            {
                setFailure(FailureCode::TRANSACTION_PROTOCOL,cell,
                    "Domain::replaceCompleteCells",
                    "sparse replacement contains a duplicate cell");
                return false;
            }
            seen[static_cast<std::size_t>(cell)] = 1;
            cellIds.push_back(cell);
        }

        std::vector<CompleteCellState> candidateCells;
        std::vector<ConservedVars> candidateUnew;
        candidateCells.reserve(replacements.size());
        candidateUnew.reserve(replacements.size());
        for (std::size_t replacementIndex = 0;
             replacementIndex < replacements.size();
             ++replacementIndex)
        {
            const CompleteCellReplacement& replacement =
                replacements[replacementIndex];
            CompleteCellConversionResult converted =
                CompleteStateConversion::fromPrimitiveAndEquilibrium(
                    replacement.primitive,
                    replacement.equilibrium,
                    *thermoModel_);
            if (!converted.success())
            {
                setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,
                    cellIds[replacementIndex],
                    "Domain::replaceCompleteCells",
                    "sparse primitive/equilibrium replacement is inconsistent",
                    converted.status,converted.diagnostics);
                return false;
            }
            candidateUnew.push_back(converted.state->U);
            candidateCells.push_back(std::move(*converted.state));
        }

        static_assert(
            std::is_nothrow_swappable<CompleteCellState>::value,
            "atomic sparse publication requires no-throw complete-state swap");
        static_assert(
            std::is_nothrow_swappable<ConservedVars>::value,
            "atomic sparse publication requires no-throw conserved-state swap");
        using std::swap;
        for (std::size_t replacementIndex = 0;
             replacementIndex < replacements.size();
             ++replacementIndex)
        {
            const int cell = cellIds[replacementIndex];
            swap(cells_[cell],candidateCells[replacementIndex]);
            swap(Unew_[cell],candidateUnew[replacementIndex]);
        }
        primitiveStageActive_ = false;
        return true;
    }

    bool replaceCompleteCells(
        const std::vector<CompleteCellReplacement>& replacements)
    {
        return replaceCompleteCellsImpl(replacements,false);
    }

    // Explicit inner-iteration boundary publication. The transaction-start
    // snapshot still owns rollback, so a rejected dual-time physical step
    // restores both interior and boundary states atomically.
    bool replaceCompleteCellsInStepTransaction(
        const std::vector<CompleteCellReplacement>& replacements)
    {
        return replaceCompleteCellsImpl(replacements,true);
    }

    // ---- Time-step lifecycle ----
    void beginStep()
    {
        if (!stateReady_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::beginStep","Domain state is not initialized");
            return;
        }
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < N(); ++k) Unew_[k] = cells_[k].U;
    }

    bool beginStepTransaction()
    {
        if (tripped_) return false;
        if (!stateReady_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::beginStepTransaction","Domain state is not initialized");
            return false;
        }
        if (stepTransactionActive_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::beginStepTransaction","nested step transaction");
            return false;
        }
        cellsStepStart_ = cells_;
        UnewStepStart_ = Unew_;
        stepTransactionActive_ = true;
        return true;
    }

    void acceptStepTransaction()
    {
        if (!stepTransactionActive_)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::acceptStepTransaction","no active step transaction");
            return;
        }
        stepTransactionActive_ = false;
        ++acceptedSteps_;
    }

    // Restore the exact transaction-start U/Unew/W/Q snapshot while keeping
    // the transaction open. Steady nonlinear line searches use this between
    // trial corrections; it is neither an accepted nor a rejected solver step.
    bool restoreStepTransactionStartState()
    {
        if(!stepTransactionActive_
            ||cellsStepStart_.size()!=cells_.size()
            ||UnewStepStart_.size()!=Unew_.size())
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::restoreStepTransactionStartState",
                "active transaction snapshot is unavailable");
            return false;
        }
        cells_=cellsStepStart_;
        Unew_=UnewStepStart_;
        primitiveStageActive_=false;
        return true;
    }

    void rollbackStepTransaction(bool countRejected = true) noexcept
    {
        if (!stepTransactionActive_) return;
        cells_.swap(cellsStepStart_);
        Unew_.swap(UnewStepStart_);
        stepTransactionActive_ = false;
        if(countRejected) ++rejectedSteps_;
    }

    bool stepTransactionActive() const { return stepTransactionActive_; }

    bool commit()
    {
        if (tripped_ || !stateReady_) return false;
        const int total = N();
        int firstInvalid = total;
        #pragma omp parallel for reduction(min:firstInvalid) schedule(static)
        for (int k = 0; k < total; ++k)
        {
            const ConservedVars& c = Unew_[k];
            bool invalid = !std::isfinite(c.rho) || c.rho <= 0.0
                || !std::isfinite(c.rhoE) || c.rhoU.size() != 3
                || c.rhoY.size() != static_cast<std::size_t>(nSpecies_)
                || !std::isfinite(c.rhoZ) || c.rhoZ < 0.0 || c.rhoZ > c.rho
                || !std::isfinite(c.rhoC) || c.rhoC < 0.0 || c.rhoC > c.rho;
            for (double momentum : c.rhoU)
                if (!std::isfinite(momentum)) invalid = true;
            if(std::isfinite(c.rho)&&c.rho>0.0)
            {
                double speciesMass=0.0;
                for(double rhoY:c.rhoY)
                {
                    if(!std::isfinite(rhoY)||rhoY < 0.0) invalid=true;
                    speciesMass+=rhoY;
                }
                if(!std::isfinite(speciesMass)
                   ||std::abs(speciesMass-c.rho)>1e-6*c.rho) invalid=true;
            }
            if (invalid) firstInvalid = std::min(firstInvalid,k);
        }
        if (firstInvalid != total)
        {
            setFailure(FailureCode::INVALID_CONSERVED_STATE,firstInvalid,
                "Domain::commit","candidate conserved state violates invariants");
            return false;
        }

        std::vector<PrimitiveVars> seeds(total,PrimitiveVars(nSpecies_));
        for (int k = 0; k < total; ++k) seeds[k] = cells_[k].W;
        std::vector<CompleteCellState> candidate;
        int failedCell = -1;
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics diagnostics;
        if (!recoverConservedField(
                Unew_,seeds,candidate,&failedCell,
                &failedStatus,&diagnostics,&cells_))
        {
            setFailure(FailureCode::STATE_RECOVERY_FAILED,failedCell,
                "Domain::commit","rho-e-Y state recovery rejected candidate",
                failedStatus,diagnostics);
            return false;
        }

        cells_.swap(candidate);
        primitiveStageActive_ = false;
        return true;
    }

    // ---- RK2 helpers ----
    void saveStateForRK2() const
    {
        if (Usaved_.size() != cells_.size())
            Usaved_.resize(cells_.size());
        for (std::size_t k = 0; k < cells_.size(); ++k)
            Usaved_[k] = cells_[k].U;
    }

    void averageWithStepStartState()
    {
        if (!stepTransactionActive_ || cellsStepStart_.size() != Unew_.size())
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::averageWithStepStartState",
                "RK2 average requires an active step transaction");
            return;
        }
        const int total = N();
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < total; ++k)
        {
            const ConservedVars& start = cellsStepStart_[k].U;
            Unew_[k].rho  = 0.5 * (start.rho  + Unew_[k].rho);
            Unew_[k].rhoE = 0.5 * (start.rhoE + Unew_[k].rhoE);
            for (std::size_t d = 0; d < Unew_[k].rhoU.size(); ++d)
                Unew_[k].rhoU[d] =
                    0.5 * (start.rhoU[d] + Unew_[k].rhoU[d]);
            for (std::size_t s = 0; s < Unew_[k].rhoY.size(); ++s)
                Unew_[k].rhoY[s] =
                    0.5 * (start.rhoY[s] + Unew_[k].rhoY[s]);
            Unew_[k].rhoZ = 0.5 * (start.rhoZ + Unew_[k].rhoZ);
            Unew_[k].rhoC = 0.5 * (start.rhoC + Unew_[k].rhoC);
        }
    }

    // ---- 3D Cyclic (periodic) BC ----
    // Apply rotational periodicity in the z (azimuthal) direction.
    // angle_deg: sector angle (e.g. 120 for 3-pintle).
    // Copies k=Nz-2 → k=0 halo and k=1 → k=Nz-1 halo,
    // rotating the tangential velocity components (U[1],U[2]) by ±angle_deg.
    void applyCyclicZImpl(double angle_deg,bool requireActiveTransaction)
    {
        if (!stateReady_
            || stepTransactionActive_ != requireActiveTransaction
            || Nz_ < 3)
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                requireActiveTransaction
                    ?"Domain::applyCyclicZInStepTransaction"
                    :"Domain::applyCyclicZ",
                requireActiveTransaction
                    ?"cyclic inner update requires an active transaction"
                    :"cyclic update requires a ready 3D state outside a transaction");
            return;
        }
        const double theta = angle_deg * M_PI / 180.0;
        const double c = std::cos(theta), s = std::sin(theta);
        std::vector<CompleteCellState> candidate = cells_;
        std::vector<ConservedVars> candidateUnew = Unew_;
        int firstInvalid = N();
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics failedDiagnostics;
        for (int i = 0; i < Nx_; ++i)
        for (int j = 0; j < Ny_; ++j)
        {
            // k=0 halo <- upper interior, rotated from +theta back to 0.
            {
                const int source = idx(i,j,Nz_-2);
                const int target = idx(i,j,0);
                PrimitiveVars w = cells_[source].W;
                double uy = w.U[1], uz = w.U[2];
                w.U[1] =  c*uy + s*uz;
                w.U[2] = -s*uy + c*uz;
                CompleteCellConversionResult converted =
                    CompleteStateConversion::fromPrimitiveAndEquilibrium(
                        w,cells_[source].Q,*thermoModel_);
                if (!converted.success())
                {
                    if (target < firstInvalid)
                    {
                        firstInvalid = target;
                        failedStatus = converted.status;
                        failedDiagnostics = std::move(converted.diagnostics);
                    }
                }
                else
                {
                    candidate[target] = std::move(*converted.state);
                    candidateUnew[target] = candidate[target].U;
                }
            }
            // k=Nz-1 halo <- lower interior, rotated from 0 to +theta.
            {
                const int source = idx(i,j,1);
                const int target = idx(i,j,Nz_-1);
                PrimitiveVars w = cells_[source].W;
                double uy = w.U[1], uz = w.U[2];
                w.U[1] =  c*uy - s*uz;
                w.U[2] =  s*uy + c*uz;
                CompleteCellConversionResult converted =
                    CompleteStateConversion::fromPrimitiveAndEquilibrium(
                        w,cells_[source].Q,*thermoModel_);
                if (!converted.success())
                {
                    if (target < firstInvalid)
                    {
                        firstInvalid = target;
                        failedStatus = converted.status;
                        failedDiagnostics = std::move(converted.diagnostics);
                    }
                }
                else
                {
                    candidate[target] = std::move(*converted.state);
                    candidateUnew[target] = candidate[target].U;
                }
            }
        }
        if (firstInvalid != N())
        {
            setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,firstInvalid,
                "Domain::applyCyclicZ","cyclic halo conversion failed",
                failedStatus,failedDiagnostics);
            return;
        }
        cells_.swap(candidate);
        Unew_.swap(candidateUnew);
        primitiveStageActive_ = false;
    }

    void applyCyclicZ(double angle_deg)
    { applyCyclicZImpl(angle_deg,false); }

    void applyCyclicZInStepTransaction(double angle_deg)
    { applyCyclicZImpl(angle_deg,true); }

    void averageWithSavedState()
    {
        const int total = N();
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < total; ++k)
        {
            Unew_[k].rho  = 0.5 * (Usaved_[k].rho  + Unew_[k].rho);
            Unew_[k].rhoE = 0.5 * (Usaved_[k].rhoE + Unew_[k].rhoE);
            for (size_t d = 0; d < Unew_[k].rhoU.size(); ++d)
                Unew_[k].rhoU[d] = 0.5 * (Usaved_[k].rhoU[d] + Unew_[k].rhoU[d]);
            for (size_t s = 0; s < Unew_[k].rhoY.size(); ++s)
                Unew_[k].rhoY[s] = 0.5 * (Usaved_[k].rhoY[s] + Unew_[k].rhoY[s]);
            Unew_[k].rhoZ = 0.5 * (Usaved_[k].rhoZ + Unew_[k].rhoZ);
            Unew_[k].rhoC = 0.5 * (Usaved_[k].rhoC + Unew_[k].rhoC);
        }
    }

    // ---- C2P / P2C (loop over all cells) ----
    bool conservedToPrimitive()
    {
        if (tripped_ || !stateReady_ || stepTransactionActive_) return false;
        std::vector<ConservedVars> source(cells_.size(),ConservedVars(nSpecies_));
        std::vector<PrimitiveVars> seeds(cells_.size(),PrimitiveVars(nSpecies_));
        for (std::size_t k = 0; k < cells_.size(); ++k)
        {
            source[k] = cells_[k].U;
            seeds[k] = cells_[k].W;
        }
        std::vector<CompleteCellState> candidate;
        int failedCell = -1;
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics diagnostics;
        if (!recoverConservedField(
                source,seeds,candidate,&failedCell,
                &failedStatus,&diagnostics,&cells_))
        {
            setFailure(FailureCode::STATE_RECOVERY_FAILED,failedCell,
                "Domain::conservedToPrimitive",
                "authoritative conserved state cannot be recovered",
                failedStatus,diagnostics);
            return false;
        }
        cells_.swap(candidate);
        primitiveStageActive_ = false;
        return true;
    }

    bool primitiveToConserved()
    {
        if (tripped_ || stepTransactionActive_ || !primitiveStageActive_)
            return false;
        std::vector<CompleteCellState> candidate;
        int failedCell = -1;
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics diagnostics;
        if (!recoverPrimitiveField(
                primitiveStage_,candidate,&failedCell,
                &failedStatus,&diagnostics))
        {
            setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,failedCell,
                "Domain::primitiveToConserved",
                "primitive initialization violates the state contract",
                failedStatus,diagnostics);
            return false;
        }
        publishCompleteField(candidate);
        return true;
    }

    bool replaceConservedField(
        const std::vector<ConservedVars>& source,
        const std::vector<double>& temperatureSeeds = {})
    {
        return replaceConservedFieldImpl(source,temperatureSeeds,false);
    }

    // Matrix-free Newton/Krylov residual evaluation needs to publish valid
    // U+epsilon*v trial fields without escaping the enclosing nonlinear-step
    // transaction. A later rollback restores the exact U/Unew/W/Q snapshot.
    bool replaceConservedFieldInStepTransaction(
        const std::vector<ConservedVars>& source,
        const std::vector<double>& temperatureSeeds = {})
    {
        return replaceConservedFieldImpl(source,temperatureSeeds,true);
    }

    bool replaceConservedField(
        const std::vector<ConservedVars>& source,
        const std::vector<thermo::EquilibriumSeed>& equilibriumSeeds)
    {
        if (tripped_ || stepTransactionActive_) return false;
        if (source.size() != cells_.size()
            || equilibriumSeeds.size() != source.size())
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::replaceConservedField",
                "replacement equilibrium-seed dimensions do not match Domain");
            return false;
        }
        std::vector<PrimitiveVars> seeds(
            source.size(),PrimitiveVars(nSpecies_));
        for (std::size_t cell = 0; cell < source.size(); ++cell)
        {
            const thermo::EquilibriumSeed& seed = equilibriumSeeds[cell];
            if (!std::isfinite(seed.pressure) || seed.pressure <= 0.0
                || !std::isfinite(seed.temperature)
                || seed.temperature <= 0.0)
            {
                setFailure(FailureCode::TRANSACTION_PROTOCOL,
                    static_cast<int>(cell),"Domain::replaceConservedField",
                    "restart equilibrium seed must contain finite positive p/T");
                return false;
            }
            seeds[cell].p = seed.pressure;
            seeds[cell].T = seed.temperature;
        }
        std::vector<CompleteCellState> candidate;
        int failedCell = -1;
        thermo::EquilibriumStatus failedStatus =
            thermo::EquilibriumStatus::INVALID_CONFIGURATION;
        thermo::EquilibriumDiagnostics diagnostics;
        if (!recoverConservedField(
                source,seeds,candidate,&failedCell,
                &failedStatus,&diagnostics,nullptr,&equilibriumSeeds))
        {
            setFailure(FailureCode::STATE_RECOVERY_FAILED,failedCell,
                "Domain::replaceConservedField",
                "replacement conserved field cannot be recovered from "
                "the authoritative p/T seeds",
                failedStatus,diagnostics);
            return false;
        }
        publishCompleteField(candidate);
        return true;
    }

    bool replaceCompleteField(
        const std::vector<PrimitiveVars>& primitive,
        const std::vector<thermo::EquilibriumState>& equilibrium)
    {
        if (tripped_ || stepTransactionActive_) return false;
        if (primitive.size() != cells_.size()
            || equilibrium.size() != cells_.size())
        {
            setFailure(FailureCode::TRANSACTION_PROTOCOL,-1,
                "Domain::replaceCompleteField",
                "replacement field dimensions do not match Domain");
            return false;
        }
        const int total = N();
        std::vector<CompleteCellState> candidate(
            total,CompleteCellState(nSpecies_));
        std::vector<thermo::EquilibriumStatus> statuses(
            total,thermo::EquilibriumStatus::INVALID_CONFIGURATION);
        std::vector<thermo::EquilibriumDiagnostics> diagnostics(total);
        int firstInvalid = total;
        #pragma omp parallel for reduction(min:firstInvalid) schedule(static)
        for (int k = 0; k < total; ++k)
        {
            CompleteCellConversionResult converted =
                CompleteStateConversion::fromPrimitiveAndEquilibrium(
                    primitive[k],equilibrium[k],*thermoModel_);
            statuses[k] = converted.status;
            diagnostics[k] = std::move(converted.diagnostics);
            if (!converted.success()) firstInvalid = std::min(firstInvalid,k);
            else candidate[k] = std::move(*converted.state);
        }
        if (firstInvalid != total)
        {
            setFailure(FailureCode::PRIMITIVE_CONVERSION_FAILED,firstInvalid,
                "Domain::replaceCompleteField",
                "replacement primitive/equilibrium field is inconsistent",
                statuses[firstInvalid],diagnostics[firstInvalid]);
            return false;
        }
        publishCompleteField(candidate);
        return true;
    }
};

} // End namespace Foam

#endif
