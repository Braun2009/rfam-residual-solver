#ifndef RealFluidState_H
#define RealFluidState_H

#ifdef RFAM_HAVE_TEQP
#include "TeqpMultifluidThermo.H"
#endif

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "HybridHelmholtzThermo.H"
#include "PengRobinsonEOS.H"
#include "SpeciesRegistry.H"

namespace Foam
{

class RealFluidState
{
public:
    enum class Status
    {
        SUCCESS,
        INVERSION_FAILED,
        NEGATIVE_CV,
        NEGATIVE_CP,
        NEGATIVE_SOUNDSPEED,
        NON_PHYSICAL
    };

private:
    // GERG-2008 residual table backend (opt-in).  Declared first so the
    // constructor can initialize gergThermo_ before any scalar members.
    const GERGResidualTable* gergTable_ = nullptr;
    const GERGSaturationTable* saturationTable_ = nullptr;
    SpeciesRegistry gergRegistry_;
    HybridHelmholtzThermo gergThermo_;
    long long gergFallbackCount_ = 0;

    double p_, T_, rho_, e_, h_, Cp_, Cv_, Z_, a_, mu_, kappa_;
    std::vector<double> Y_;
    // Opt-in species-aware NASA-7 thermo (Stage 3 step 1). speciesNames_ is the
    // caller's exact rhoY / global-registry order; empty + useNasa7_==false =>
    // legacy constant-Cp placeholder, so the MVP path stays byte-for-byte.
    // Set only via setSpecies(). [[NASA7Thermo]]
    std::vector<std::string> speciesNames_;
    bool useNasa7_;
    double dpdrho_T_, dpdT_rho_, dedrho_T_, dedT_rho_;
    PengRobinsonEOS eos_;
    bool cached_;
    Status status_;

    static std::string activePureSpecies(
        const std::vector<std::string>& names,
        const std::vector<double>& Y, double threshold = 0.999);
    Status updatePTYGERG(double p, double T, const std::vector<double>& Y);
    Status updateRhoTYGERG(double rho, double T, const std::vector<double>& Y);
    Status fillFromHybrid(const HybridHelmholtzThermo::State& st);
    // Minimal HEM: two-phase recovery for a (rho, e) state inside the dome.
    // Returns SUCCESS with twoPhase_=true, or NON_PHYSICAL.
    Status updateRhoEYTwoPhase(double rho, double e,
                               const std::vector<double>& Y,
                               bool allowFullScan = true);
    // Condensed-phase (solid CO2/H2O) recovery for cold diluted states
    // whose single-phase root is mechanically unstable or absent (e.g. a
    // 185 K O2/C3H8/H2O/CO2 mixture: H2O freezes, CO2 sublimes).  The gas
    // EOS operates on the O2/C3H8/... sub-composition; the solid supplies
    // volume, specific heat and formation-free mass/energy.  Returns
    // SUCCESS with condensed_=true, or NON_PHYSICAL/INVERSION_FAILED.
    Status updateRhoEYCondensed(double rho, double e,
                                const std::vector<double>& Y);
    mutable bool condensed_ = false;
    mutable double condensedMassFraction_ =
        std::numeric_limits<double>::quiet_NaN();
    mutable std::vector<double> condensedMassFractions_;  // per-species
    mutable std::vector<double> gasMassFractions_;        // gas sub-comp.
    mutable double gasRho_ = std::numeric_limits<double>::quiet_NaN();
    mutable double gasE_ = std::numeric_limits<double>::quiet_NaN();
    bool saturationRoots(const std::vector<std::string>& names,
                         const std::vector<double>& Ypure,
                         double T, double pSat, double& rhoL, double& rhoV,
                         double& eL, double& eV) const;
    bool inSaturationDome(double rho, double e, double xO2) const;
    mutable bool twoPhase_ = false;
    mutable bool eProjected_ = false;
    mutable double eProjectionDelta_ = 0.0;
    mutable bool twoPhaseHintValid_ = false;
    mutable double twoPhaseHintT_ = std::numeric_limits<double>::quiet_NaN();
    mutable double twoPhaseHintQ_ = std::numeric_limits<double>::quiet_NaN();
    // Composition roundoff repair audit (2026-08-07): set when the
    // conserved-field mass noise was clipped/renormalised before the EOS
    // closure.  Exposed so the model layer can mark the backend status.
    mutable bool repairedComposition_ = false;
    // Energy-not-bracketed diagnostics (2026-08-07, cold-region class 2)
    mutable bool energyNotBracketedDiagnostic_ = false;
    mutable double energyNotBracketedColdGap_ = 0.0;
    mutable double energyNotBracketedDTEq_ = 0.0;
    // Mechanical-stability diagnostic (2026-08-07, cold-region class 3):
    // set when the committed single-phase root has (dp/drho)_T <= 0
    // (metastable/unstable branch); recorded, not fatal.
    mutable bool phaseUnstable_ = false;
    mutable double quality_ = std::numeric_limits<double>::quiet_NaN();
    mutable double pSat_ = std::numeric_limits<double>::quiet_NaN();

    static constexpr double SMALL_DENSITY = 1e-8;
    static constexpr double SMALL_PRESSURE = 1e3;
    static constexpr double SMALL_TEMPERATURE = 1.0;

    bool computeDerivatives();
    Status enforceBounds(double soundSpeedSquared);

public:
    RealFluidState();

    Status updatePTY(double p, double T, const std::vector<double>& Y);
    Status updateRhoTY(double rho, double T, const std::vector<double>& Y);
    Status updateRhoEY(double rho, double e, const std::vector<double>& Y, double T0 = 300.0);

    static Status validateThermodynamicScalars(
        double Cv, double Cp, double soundSpeedSquared);

    // Phase 7 Step 2: make this state (its EOS/mixing rule) species-aware over
    // the caller's exact species order — for reacting flow this is the global
    // SpeciesRegistry order, i.e. the SAME order as ConservedVars::rhoY / Y.
    // OPT-IN: the default constructor stays 2-species {O2,C3H8}, so callers that
    // never call setSpecies (e.g. the MVP) are byte-for-byte unchanged. Throws
    // on an unknown/non-physical species (no silent fallback).
    void setSpecies(const std::vector<std::string>& names)
    {
        eos_.setSpeciesNames(names);
        speciesNames_ = names;   // remember exact order for NASA-7 mixture thermo
        useNasa7_ = true;        // opt in to real Cp(T)/h(T); legacy path stays 1000
        if (Y_.size() != names.size()) Y_.assign(names.size(), 0.0);
    }

    // Register the GERG-2008 residual table (hot-start production backend).
    // OPT-IN: without it the PR path (including its two-phase-dome fail-close)
    // is unchanged.  With it, single-phase pure/near-pure states use the GERG
    // residual (defined for every rho, including the dense liquid branch the
    // PR isochoric path cannot represent).  PR fallback after a GERG failure
    // is never silent: it increments gergFallbackCount() for explicit audit.
    void setGERGTable(const GERGResidualTable* table);
    long long gergFallbackCount() const { return gergFallbackCount_; }

    // Register the GERG saturation-line table (bubble pressures).  With both
    // tables registered, a single-phase recovery failure that lands inside
    // the two-phase dome is answered by the minimal HEM path instead of
    // failing: T -> Tsat(p or e), quality from the energy lever rule.
    // The HEM result carries a two-phase flag (isTwoPhase()) and quality.
    void setSaturationTable(const GERGSaturationTable* table);
    // Two-phase continuation hint: the previous (T, quality) of this cell
    // seeds the next HEM recovery (fast local scan).  A hint is NEVER
    // authoritative: if the local scan fails, the full stability scan
    // falls back, and the final state still passes the residual/density
    // checks.  Cleared automatically when the hint misses.
    void setTwoPhaseHint(double T, double quality)
    {
        twoPhaseHintT_ = T;
        twoPhaseHintQ_ = quality;
        twoPhaseHintValid_ = std::isfinite(T) && std::isfinite(quality);
    }
    bool hasTwoPhaseHint() const { return twoPhaseHintValid_; }
    bool isTwoPhase() const { return twoPhase_; }
    bool isCondensed() const { return condensed_; }
    double condensedMassFraction() const { return condensedMassFraction_; }
    const std::vector<double>& condensedMassFractions() const
    { return condensedMassFractions_; }
    const std::vector<double>& gasMassFractions() const
    { return gasMassFractions_; }
    double gasRho() const { return gasRho_; }
    double gasE() const { return gasE_; }
    bool energyProjected() const { return eProjected_; }
    double energyProjectionDelta() const { return eProjectionDelta_; }
    // Recovery profiling (2026-08-04): counters for the rho-e-Y recovery
    // path — Newton iterations vs bracket fallback, and HEM full scans.
    struct RecoveryProfile
    {
        long long recoveryCalls = 0;        // updateRhoEY invocations
        long long newtonIterations = 0;     // total Newton iterations
        long long newtonSuccesses = 0;      // solved by Newton alone
        long long bracketFallbacks = 0;     // bracket path used
        long long bracketIterations = 0;    // total bracket (bisection) iters
        long long hemFullScans = 0;         // HEM full stability scans
        long long hemHintHits = 0;          // HEM continuation-hint hits
        double newtonAvg = 0.0;             // avg Newton iters per Newton solve
        double bracketRate = 0.0;           // bracket fallback rate
    };
    mutable RecoveryProfile profile_;
    const RecoveryProfile& recoveryProfile() const { return profile_; }
    void resetRecoveryProfile() const { profile_ = RecoveryProfile{}; }


    double quality() const { return quality_; }
    double pSat() const { return pSat_; }

    // Mixture molecular weight [kg/mol] for the currently stored Y (diagnostic;
    // used by the species-mapping smoke test to verify setSpecies routing).
    double MW() const { return eos_.MW(Y_); }
#ifdef RFAM_HAVE_TEQP
    // GERG-2008 multifluid mixture recovery (teqp).  Lazily constructed
    // singleton; data root from RFAM_TEQP_DATA env or the default.
    static const TeqpMultifluidThermo& teqpInstance();
    static bool teqpAvailable();
#endif
    // True iff the last recovery clipped/renormalised roundoff noise
    bool repairedComposition() const { return repairedComposition_; }
    bool energyNotBracketed() const { return energyNotBracketedDiagnostic_; }
    double energyNotBracketedColdGap() const
    { return energyNotBracketedColdGap_; }
    double energyNotBracketedDTEq() const
    { return energyNotBracketedDTEq_; }
    bool phaseUnstable() const { return phaseUnstable_; }
    double mwFor(const std::vector<double>& Y) const
    { return eos_.MW(Y); }

    // Thermodynamic properties
    double p() const { return p_; }
    double T() const { return T_; }
    double rho() const { return rho_; }
    double e() const { return e_; }
    double h() const { return h_; }
    double Cp() const { return Cp_; }
    double Cv() const { return Cv_; }
    double Z() const { return Z_; }
    double a() const { return a_; }
    double gamma() const { return Cp_ / Cv_; }

    // Thermodynamic derivatives (for preconditioning)
    double dpdrho_T() const { return dpdrho_T_; }
    double dpdT_rho() const { return dpdT_rho_; }
    double dedrho_T() const { return dedrho_T_; }
    double dedT_rho() const { return dedT_rho_; }

    // Flow properties
    double getMach(const std::vector<double>& U) const {
        double u_mag = std::sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
        return u_mag / (a_ + 1e-10);
    }

    // Diagnostics
    Status status() const { return status_; }
    bool isValid() const { return status_ == Status::SUCCESS; }
    // Reduced coordinates based on the EOS pseudo-critical mixing rule.
    // These are diagnostics, not a multicomponent phase-envelope definition.
    double Tr() const;
    double Pr() const;
    bool isNearCritical(double Tr_tol = 0.05, double Pr_tol = 0.2) const;

    void print(std::ostream& os) const;
};

} // End namespace Foam

#endif
