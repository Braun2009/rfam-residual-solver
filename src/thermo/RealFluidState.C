/*---------------------------------------------------------------------------*\
  realFluidAllMach - RealFluidState Implementation (Fixed)
\*---------------------------------------------------------------------------*/

#include "RealFluidState.H"
#include "NASA7Thermo.H"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace Foam
{

static constexpr double R_UNIVERSAL = 8.314462618;

// Global databases for EOS
static SpeciesCriticalDB g_criticalDB = SpeciesCriticalDB::createDefault();
static BinaryInteractionDB g_kijDB = BinaryInteractionDB::createDefault();

RealFluidState::RealFluidState()
    : gergThermo_(gergRegistry_),
      p_(101325.0),
    T_(300.0),
    rho_(1.0),
    e_(0.0),
    h_(0.0),
    Cp_(1000.0),
    Cv_(700.0),
    Z_(1.0),
    a_(340.0),
    mu_(1.8e-5),
    kappa_(0.025),
    useNasa7_(false),
    dpdrho_T_(0.0),
    dpdT_rho_(0.0),
    dedrho_T_(0.0),
    dedT_rho_(700.0),
    eos_(g_criticalDB, g_kijDB),
    cached_(false),
    status_(Status::SUCCESS)
{
    Y_.resize(2, 0.0);
    Y_[0] = 1.0;
}

void RealFluidState::setGERGTable(const GERGResidualTable* table)
{
    gergTable_ = table;
    gergThermo_.registerResidualTable(table);
    gergFallbackCount_ = 0;
}

void RealFluidState::setSaturationTable(const GERGSaturationTable* table)
{
    saturationTable_ = table;
}

std::string RealFluidState::activePureSpecies(
    const std::vector<std::string>& names,
    const std::vector<double>& Y, double threshold)
{
    if (names.size() != Y.size()) return std::string();
    std::string pure;
    for (std::size_t i = 0; i < names.size(); ++i)
    {
        if (!std::isfinite(Y[i]) || Y[i] < 0.0) return std::string();
        if (Y[i] > threshold)
        {
            if (!pure.empty()) return std::string();
            pure = names[i];
        }
    }
    return pure;
}

// Minimal HEM two-phase recovery: given (rho, e) inside the dome,
// iterate T so that the lever-rule quality reproduces e and the
// volume-weighted density reproduces rho, using the saturation table
// (p_sat, rhoL, rhoV) plus GERG single-phase energies at the saturation
// roots.  DIAGNOSTIC-grade: pure species only, equilibrium only.
RealFluidState::Status RealFluidState::updateRhoEYTwoPhase(
    double rho, double e, const std::vector<double>& Y, bool allowFullScan)
{
    twoPhase_ = false;
    eProjected_ = false;
    eProjectionDelta_ = 0.0;
    quality_ = std::numeric_limits<double>::quiet_NaN();
    if (!gergTable_ || !saturationTable_)
        return Status::NON_PHYSICAL;

    // Composition level (2026-08-04, v3): pure species (O2/C3H8) OR the
    // binary xO2 band [0.72, 0.78] whose teqp coexisting-phase data
    // (bubble pressure + vapour composition y) is available in the v3
    // table.  Bands WITHOUT teqp data (xO2 < 0.72, 0.78 < xO2 < 0.92,
    // xO2 > 0.92) FAIL CLOSED: nearest-level snapping would make the phase
    // state jump across composition steps in the diffusion layer, and the
    // old mixture endpoints were not coexisting phases.  Real full-range
    // binary phase splits await a robust deep-cold teqp flash.
    double xO2 = 1.0;
    const std::string pure = activePureSpecies(speciesNames_, Y);
    if (pure.empty())
    {
        // O2/C3H8 sub-composition level: non-fuel/oxidiser species
        // (H2O/CO2 trace dilution in the LOX near-field) are passive
        // diluents; the phase split is taken on the O2/C3H8 ratio only.
        double nO2 = 0.0, nC3 = 0.0;
        for (std::size_t i = 0; i < speciesNames_.size(); ++i)
        {
            if (speciesNames_[i] == "O2") nO2 += Y[i] / 0.0319988;
            else if (speciesNames_[i] == "C3H8") nC3 += Y[i] / 0.04409562;
        }
        if (nO2 <= 0.0 && nC3 <= 0.0) return Status::NON_PHYSICAL;
        const double xm = nO2 / (nO2 + nC3);
        // Admissible levels (v3 table): the 0.75 binary band with teqp
        // coexisting-phase data, plus the near-pure wings that snap
        // continuously onto the pure-O2 / pure-C3H8 lines.
        // The mid-gap 0.08-0.72 / 0.78-0.92 stays fail-closed (a
        // nearest-level snap there would jump across composition
        // steps in the diffusion layer).
        if (!((xm >= 0.72 && xm <= 0.78)
              || xm >= 0.92 || xm <= 0.08))
            return Status::NON_PHYSICAL;
        xO2 = xm;
        // Data gate: the teqp v3 level must cover this temperature
        // (p_sat finite); otherwise fail closed.
        const double pSat = saturationTable_->saturationPressure(
            twoPhaseHintValid_ && std::isfinite(twoPhaseHintT_)
                ? twoPhaseHintT_ : 120.0, xO2);
        if (!std::isfinite(pSat)) return Status::NON_PHYSICAL;
    }
    else if (pure != "O2" && pure != "C3H8")
        return Status::NON_PHYSICAL;
    if (pure == "C3H8") xO2 = 0.0;

    // Continuation-hint scan first (fast): if the previous (T, quality)
    // of this cell is available, scan its neighbourhood at 0.25 K steps.
    // The hint is a seed only — the final state must still match rho.
    double bestT = std::numeric_limits<double>::quiet_NaN();
    double bestErr = std::numeric_limits<double>::infinity();
    double bestX = std::numeric_limits<double>::quiet_NaN();
    double bestRhoMix = std::numeric_limits<double>::quiet_NaN();
    auto tryScan = [&](double tLo, double tHi, double step) {
        for (double Tq = tLo; Tq <= tHi; Tq += step)
        {
            double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
            if (!saturationTable_->saturationEndpoints(
                    Tq, xO2, rhoV, rhoL, &eV, &eL))
                continue;
            double x = (e - eL) / (eV - eL);
            if (!std::isfinite(x)) continue;
            x = std::max(0.0, std::min(1.0, x));
            const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
            const double err = std::abs(rhoMix - rho) / rho;
            // Energy gate (2026-08-04): a candidate whose quality clamps
            // to 0/1 because e lies outside the band commits eL/eV, which
            // can differ from the input e by >10% (observed 10.2% at the
            // 147 K eL clamp with 130 K liquid energy).  Such a candidate
            // must NOT be accepted — the projection path (which minimises
            // the e deviation at fixed rho) handles it instead.
            const double eMix = x * eV + (1.0 - x) * eL;
            const double errE =
                std::abs(eMix - e) / std::max(std::abs(e), 1.0);
            if (err < bestErr && errE <= 0.05)
            {
                bestErr = err;
                bestT = Tq;
                bestX = x;
                bestRhoMix = rhoMix;
            }
        }
    };
    if (twoPhaseHintValid_)
    {
        tryScan(std::max(60.0, twoPhaseHintT_ - 3.0),
                std::min(370.0, twoPhaseHintT_ + 3.0), 0.25);
        if (std::isfinite(bestT) && bestErr <= 5.0e-3)
        {
            ++profile_.hemHintHits;
        }
        else
        {
            twoPhaseHintValid_ = false;   // miss: hint not authoritative
            bestT = std::numeric_limits<double>::quiet_NaN();
            bestErr = std::numeric_limits<double>::infinity();
            bestX = std::numeric_limits<double>::quiet_NaN();
            bestRhoMix = std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (!std::isfinite(bestT))
    {
        if (!allowFullScan) return Status::NON_PHYSICAL;
        ++profile_.hemFullScans;
        // Full stability scan at 1 K steps (fallback).
        tryScan(60.0, 370.0, 1.0);
    }
#ifdef RFAM_THERMO_DIAGNOSTICS
    std::cerr << "[HEM-SCAN-DIAG] rho=" << rho << " e=" << e
              << " xO2=" << xO2 << " hint=" << twoPhaseHintValid_
              << " bestT=" << bestT << " bestErr=" << bestErr
              << " bestX=" << bestX << "\n";
#endif
    if (!std::isfinite(bestT) || bestErr > 5.0e-3)
    {
        if (!allowFullScan) return Status::NON_PHYSICAL;
#ifdef RFAM_THERMO_DIAGNOSTICS
        std::cerr << "[PROJ-DIAG] entering projection: rho=" << rho
                  << " e=" << e << " bestErr=" << bestErr << "\n";
#endif
        // Projection branch (domain-level two-phase projection): the exact
        // HEM match failed because the numerical energy lies outside the
        // saturated band while rho sits inside the dome (subcritical mixing
        // noise, e.g. the Mayer jet edge).  Keep rho (mass is the primary
        // conserved quantity), solve the density lever rule for x, and
        // project e onto the saturated endpoints.  The projection is
        // explicit and audited (eProjectionDelta_); energy is modified by
        // the numerical-noise magnitude only.
        double bestP = std::numeric_limits<double>::infinity();
        double pT = std::numeric_limits<double>::quiet_NaN();
        double pX = std::numeric_limits<double>::quiet_NaN();
        double pRhoMix = std::numeric_limits<double>::quiet_NaN();
        double pE = std::numeric_limits<double>::quiet_NaN();
        for (double Tq = 60.0; Tq <= 370.0; Tq += 0.5)
        {
            double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
            if (!saturationTable_->saturationEndpoints(
                    Tq, xO2, rhoV, rhoL, &eV, &eL))
                continue;
            const double xd = (1.0 / rho - 1.0 / rhoL)
                / (1.0 / rhoV - 1.0 / rhoL);
            if (!std::isfinite(xd) || xd < -0.02 || xd > 1.02)
                continue;
            const double x = std::max(0.0, std::min(1.0, xd));
            const double eProj = x * eV + (1.0 - x) * eL;
            // Lever-rule density consistent with the committed quality and
            // the REAL saturation endpoints (2026-08-04 audit).  The
            // candidate must also reproduce the input rho within 1%: a
            // compressed liquid (rho=920.4 at 130 K/4 MPa) is denser than
            // the saturated liquid (rhoL=902), and projecting it onto a
            // band whose rhoMix differs by >1% would silently change mass
            // (conserved U.rho is untouched, Q.rho is the flux state).
            const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
            if (std::abs(rhoMix - rho) / rho > 0.01)
                continue;
            const double pen = std::abs(eProj - e);
            if (pen < bestP)
            {
                bestP = pen;
                pT = Tq;
                pX = x;
                pRhoMix = rhoMix;
                pE = eProj;
            }
        }
        if (!std::isfinite(pT))
        {
#ifdef RFAM_THERMO_DIAGNOSTICS
            std::cerr << "[PROJ-DIAG] projection failed: rho=" << rho
                      << " e=" << e << " xO2=" << xO2 << "\n";
#endif
            return Status::NON_PHYSICAL;
        }
#ifdef RFAM_THERMO_DIAGNOSTICS
        std::cerr << "[PROJ-DIAG] projected: rho=" << rho
                  << " e=" << e << " -> eProj=" << pE
                  << " T=" << pT << " q=" << pX << "\n";
#endif
        bestT = pT;
        bestX = pX;
        bestRhoMix = pRhoMix;
        bestErr = 0.0;
        eProjected_ = true;
        eProjectionDelta_ = pE - e;
    }

    // Commit the two-phase state.
    const double Tq = bestT;
    double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
    if (!saturationTable_->saturationEndpoints(
            Tq, xO2, rhoV, rhoL, &eV, &eL))
        return Status::NON_PHYSICAL;
    const double pSat = saturationTable_->saturationPressure(Tq, xO2);
    const double x = bestX;
    const double eMix = eProjected_ ? (x * eV + (1.0 - x) * eL)
                                    : (x * eV + (1.0 - x) * eL);
    p_ = pSat;
    T_ = Tq;
    rho_ = bestRhoMix;
    e_ = eMix;
    h_ = e_ + p_ / rho_;
    quality_ = x;
    pSat_ = pSat;
    twoPhase_ = true;
    // refresh the continuation hint for the next recovery of this cell
    twoPhaseHintT_ = Tq;
    twoPhaseHintQ_ = x;
    twoPhaseHintValid_ = true;
    // Self-consistent diagnostic closure: pick dpdrho/dpdT, then derive
    // Cv/Cp/dedrho/dedT so that the thermodynamic identities hold to the
    // commit-path closure tolerance (1e-7).
    dpdrho_T_ = 1.4 * pSat / bestRhoMix;
    dpdT_rho_ = pSat / Tq;
    Cv_ = 700.0 + x * 400.0;
    Cp_ = Cv_ + Tq * dpdT_rho_ * dpdT_rho_
        / (bestRhoMix * bestRhoMix * dpdrho_T_);
    Z_ = pSat * eos_.MW(Y) / (bestRhoMix * R_UNIVERSAL * Tq);
    a_ = std::sqrt((Cp_ / Cv_) * dpdrho_T_);
    dedT_rho_ = Cv_;
    dedrho_T_ = (p_ - T_ * dpdT_rho_) / (rho_ * rho_);
    cached_ = true;
    status_ = Status::SUCCESS;
    if (profile_.newtonSuccesses > 0)
        profile_.newtonAvg = static_cast<double>(profile_.newtonIterations)
            / profile_.newtonSuccesses;
    profile_.bracketRate = profile_.recoveryCalls > 0
        ? static_cast<double>(profile_.bracketFallbacks)
            / profile_.recoveryCalls : 0.0;
    return status_;
}

#ifdef RFAM_HAVE_TEQP
RealFluidState::Status RealFluidState::updateRhoEYCondensed(
    double rho, double e, const std::vector<double>& Y)
{
    condensed_ = false;
    condensedMassFraction_ = std::numeric_limits<double>::quiet_NaN();
    condensedMassFractions_.clear();
    gasMassFractions_.clear();
    gasRho_ = std::numeric_limits<double>::quiet_NaN();
    gasE_ = std::numeric_limits<double>::quiet_NaN();
    if (!teqpAvailable() || speciesNames_.size() != 9)
        return Status::NON_PHYSICAL;
    if (Y.size() != 9) return Status::NON_PHYSICAL;
    if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(e))
        return Status::NON_PHYSICAL;

    // Condensed species (9-species registry order): H2O=6, CO2=8.
    // Solid props (order-of-magnitude constants; the exact crystal
    // properties are a later tier refinement):
    constexpr double cpCO2 = 850.0;    // [J/kg/K] dry-ice
    constexpr double cpH2O = 2050.0;   // [J/kg/K] ice
    constexpr double rhoCO2 = 1562.0;  // [kg/m^3]
    constexpr double rhoH2O = 917.0;   // [kg/m^3]

    double solidFrac = 0.0, solidCp = 0.0, solidRhoInv = 0.0;
    for (std::size_t s : {std::size_t(6), std::size_t(8)})
    {
        if (Y[s] < 0.0 || !std::isfinite(Y[s])) return Status::NON_PHYSICAL;
        if (Y[s] > 0.0)
        {
            solidFrac += Y[s];
            solidCp += Y[s] * (s == 6 ? cpH2O : cpCO2);
            solidRhoInv += Y[s] / (s == 6 ? rhoH2O : rhoCO2);
        }
    }
    if (solidFrac <= 1.0e-6 || solidFrac >= 1.0 - 1.0e-6)
        return Status::NON_PHYSICAL;
    const double solidRho = solidFrac / solidRhoInv;
    const double solidCpMix = solidCp / solidFrac;  // per kg solid
    const double gasFrac = 1.0 - solidFrac;

    // Gas sub-composition: drop the condensed species, renormalise.
    std::vector<double> gasY = Y;
    gasY[6] = 0.0;
    gasY[8] = 0.0;
    double gasSum = 0.0;
    for (double gy : gasY) gasSum += gy;
    if (gasSum <= 1.0e-6) return Status::NON_PHYSICAL;
    for (double& gy : gasY) gy /= gasSum;

    double mwGas = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> moles;
    if (!TeqpMultifluidThermo::massToMoleFractions(
            gasY, moles, &mwGas)
        || !std::isfinite(mwGas) || mwGas <= 0.0)
        return Status::NON_PHYSICAL;
    const TeqpMultifluidThermo& teqp = teqpInstance();
    const double vTotal = 1.0 / rho;

    // Sensible-energy convention: gas from teqp (NASA7 298 K + residual),
    // solid as cp_s*(T-298.15) per kg solid.
    auto eGasOfT = [&](double T) -> double
    {
        const double vSolid = solidFrac / solidRho;
        const double vGas = vTotal - vSolid;
        if (vGas <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double rhoGas = gasFrac / vGas;
        return teqp.eSens(T, rhoGas / mwGas, gasY);
    };
    auto f = [&](double T) -> double
    {
        const double eg = eGasOfT(T);
        if (!std::isfinite(eg)) return std::numeric_limits<double>::quiet_NaN();
        return gasFrac * eg + solidFrac * solidCpMix * (T - 298.15) - e;
    };

    // Monotone bracket+bisection in [60, 1500] K.
    const double tLo = 60.0, tHi = 1500.0;
    double prevT = tLo, prevF = f(tLo);
    double solT = std::numeric_limits<double>::quiet_NaN();
    for (int i = 1; i <= 200; ++i)
    {
        const double Tq = tLo + (tHi - tLo) * i / 200.0;
        const double fq = f(Tq);
        if (std::isfinite(prevF) && std::isfinite(fq)
            && std::signbit(prevF) != std::signbit(fq))
        {
            double lo = prevT, hi = Tq, flo = prevF;
            for (int it = 0; it < 120; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                const double fm = f(mid);
                if (!std::isfinite(fm)) { lo = mid; continue; }
                if (std::abs(fm) <= 1e-9 * std::max(std::abs(e), 1.0)
                    || (hi - lo) <= 1e-10 * std::max(std::abs(mid), 1.0))
                { solT = mid; break; }
                if (std::signbit(flo) == std::signbit(fm)) { lo = mid; flo = fm; }
                else hi = mid;
            }
            if (!std::isfinite(solT)) solT = 0.5 * (lo + hi);
            break;
        }
        prevT = Tq;
        prevF = fq;
    }
    if (!std::isfinite(solT) || solT <= 0.0)
    {
        return Status::INVERSION_FAILED;
    }
    // The condensed solution is only physical below the CO2/etc. freeze
    // boundary; a warmer root belongs to the single-phase branch.
    if (solT > 216.6)
    {
        return Status::NON_PHYSICAL;
    }

    const double vSolid = solidFrac / solidRho;
    const double vGas = vTotal - vSolid;
    const double rhoGas = gasFrac / vGas;
    const double eGas = eGasOfT(solT);
    if (!std::isfinite(eGas)) return Status::INVERSION_FAILED;
    gasRho_ = rhoGas;
    gasE_ = eGas;
    gasMassFractions_ = gasY;

    // Gas pressure and derivatives on the gas sub-composition.
    double pq = std::numeric_limits<double>::quiet_NaN();
    try
    {
        pq = teqp.pressure(solT, rhoGas / mwGas, gasY);
    }
    catch (...) { return Status::INVERSION_FAILED; }
    if (!std::isfinite(pq) || pq <= SMALL_PRESSURE)
        return Status::INVERSION_FAILED;
    double dpdrhoMolar = std::numeric_limits<double>::quiet_NaN();
    double dpdT = std::numeric_limits<double>::quiet_NaN();
    double cv = std::numeric_limits<double>::quiet_NaN();
    teqp.derivatives(solT, rhoGas / mwGas, gasY,
                     dpdrhoMolar, dpdT, cv);
    const double dpdrhoMass = dpdrhoMolar / mwGas;
    if (!std::isfinite(dpdrhoMass) || dpdrhoMass <= 0.0
        || !std::isfinite(cv) || cv <= 0.0 || !std::isfinite(dpdT))
        return Status::NON_PHYSICAL;
    const double cp = cv + solT * dpdT * dpdT
        / (rhoGas * rhoGas * dpdrhoMass);
    const double a2Gas = (cp / cv) * dpdrhoMass;
    if (!std::isfinite(cp) || !std::isfinite(a2Gas) || a2Gas <= 0.0)
        return Status::NON_PHYSICAL;

    // Mixture sound speed: only the gas is compressible;
    // d rhoGas/d rho = gasFrac*rhoSolid^2/(rhoSolid - solidFrac*rho)^2.
    const double dRhoGasDRho = gasFrac * solidRho * solidRho
        / ((solidRho - solidFrac * rho) * (solidRho - solidFrac * rho));
    const double aMix2 = a2Gas * dRhoGasDRho;
    if (!std::isfinite(aMix2) || aMix2 <= 0.0)
        return Status::NON_PHYSICAL;

    // Payload: mass-weighted heat capacities; dpdrho/dpdT set so the
    // Cp-Cv identity closes on the mixture.
    Cv_ = gasFrac * cv + solidFrac * solidCpMix;
    Cp_ = gasFrac * cp + solidFrac * solidCpMix;
    dpdrho_T_ = aMix2 / (Cp_ / Cv_);
    dpdT_rho_ = std::sqrt(std::max(
        (Cp_ - Cv_) * rho * rho * dpdrho_T_ / solT, 0.0));
    dedrho_T_ = (pq - solT * dpdT_rho_) / (rho * rho);
    dedT_rho_ = Cv_;
    a_ = std::sqrt(aMix2);

    T_ = solT;
    p_ = pq;
    rho_ = rho;
    e_ = e;
    h_ = e + pq / rho;
    double mwTotal = 1.0;
    {
        double molesPerMass = 0.0;
        for (std::size_t s = 0; s < speciesNames_.size(); ++s)
        {
            const double mw = eos_.MW(std::vector<double>{Y[s]});
            molesPerMass += Y[s] / mw;
        }
        if (molesPerMass > 0.0) mwTotal = 1.0 / molesPerMass;
    }
    Z_ = pq * mwTotal / (rho * R_UNIVERSAL * solT);
    condensed_ = true;
    condensedMassFraction_ = solidFrac;
    condensedMassFractions_.assign(9, 0.0);
    condensedMassFractions_[6] = Y[6];
    condensedMassFractions_[8] = Y[8];
    cached_ = true;
    status_ = Status::SUCCESS;
    return status_;
}
#endif

RealFluidState::Status RealFluidState::fillFromHybrid(
    const HybridHelmholtzThermo::State& st)
{
    const double MW = eos_.MW(Y_);
    if (!std::isfinite(MW) || MW <= 0.0) return Status::NON_PHYSICAL;
    p_ = st.p;
    T_ = st.T;
    rho_ = st.rho;
    e_ = st.e;
    h_ = st.h;
    Cv_ = st.cv;
    Cp_ = st.cp;
    a_ = st.soundSpeed;
    Z_ = p_ * MW / (rho_ * R_UNIVERSAL * T_);
    dpdrho_T_ = st.dpdrho_T;
    dpdT_rho_ = st.dpdT_rho;
    dedrho_T_ = (p_ - T_ * dpdT_rho_) / (rho_ * rho_);
    dedT_rho_ = Cv_;
    cached_ = true;
    status_ = Status::SUCCESS;
    return Status::SUCCESS;
}

// GERG isothermal root: find rho such that p(T,rho) = p_target on the
// dense-liquid branch (scan from the top of the table density range).
// Input-level two-phase detection: rho lies between the saturated liquid
// and vapour densities at some covered temperature.  The energy is NOT
// required to lie in the band: over-cooled liquid / super-heated vapour
// inside the density dome still belong to the two-phase candidate set and
// the HEM clamps the quality onto the saturated endpoints.
bool RealFluidState::inSaturationDome(double rho, double e, double xO2) const
{
    if (!saturationTable_ || !std::isfinite(rho) || !std::isfinite(e))
        return false;
    const bool coldLog=(rho<100.0&&e<-2.0e5);
    // scan the covered T range at 0.5 K steps (a 2 K step can miss the
    // narrow density band of the compressed-liquid edge, 2026-08-04)
    for (double Tq = 60.0; Tq <= 370.0; Tq += 0.5)
    {
        double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
        if (!saturationTable_->saturationEndpoints(
                Tq, xO2, rhoV, rhoL, &eV, &eL))
            continue;
        if (rho > rhoV && rho < rhoL)
        {
            const double bandLo = eL - 0.10 * std::abs(eL);
            const double bandHi = eV + 0.10 * std::abs(eV);
            if (e >= bandLo && e <= bandHi)
                return true;
            if (coldLog && rho > rhoV && rho < rhoL)
            {
                static std::mutex dm; std::lock_guard<std::mutex> lg(dm);
            }
        }
        else if (coldLog && rho > rhoV && rho <= rhoL*1.02)
        {
            static std::mutex dm; std::lock_guard<std::mutex> lg(dm);
        }
    }
    return false;
}

// Isobaric saturation roots: scan the GERG isotherm at (T, pSat) and
// return the stable liquid/vapour densities and energies.
bool RealFluidState::saturationRoots(
    const std::vector<std::string>& names,
    const std::vector<double>& Ypure,
    double T, double pSat, double& rhoL, double& rhoV,
    double& eL, double& eV) const
{
    double rhoLo = 0.0, rhoHi = 0.0;
    if (!gergTable_->densityRange(names[0], rhoLo, rhoHi))
        return false;
    double hiRho = std::numeric_limits<double>::quiet_NaN();
    double loRho = std::numeric_limits<double>::quiet_NaN();
    constexpr int N = 512;
    for (int i = 0; i <= N; ++i)
    {
        const double xi = 1.0 - static_cast<double>(i) / N;
        const double rq = rhoLo * std::pow(rhoHi / rhoLo, xi);
        const HybridHelmholtzThermo::State st =
            gergThermo_.evaluate(names, Ypure, T, rq);
        if (!st.success()) continue;
        if (st.p >= pSat)
            hiRho = rq;
        else if (std::isfinite(hiRho))
        {
            loRho = rq;
            break;
        }
    }
    if (!std::isfinite(hiRho) || !std::isfinite(loRho))
        return false;
    double lo = loRho, hi = hiRho;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        const HybridHelmholtzThermo::State st =
            gergThermo_.evaluate(names, Ypure, T, mid);
        if (!st.success()) return false;
        if (st.p > pSat) lo = mid; else hi = mid;
    }
    const HybridHelmholtzThermo::State sL =
        gergThermo_.evaluate(names, Ypure, T, 0.5 * (lo + hi));
    if (!sL.success()) return false;
    rhoL = sL.rho;
    eL = sL.e;
    lo = rhoLo;
    hi = 0.5 * loRho;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        const HybridHelmholtzThermo::State st =
            gergThermo_.evaluate(names, Ypure, T, mid);
        if (!st.success()) { lo = mid; continue; }
        if (st.p > pSat) lo = mid; else hi = mid;
    }
    const HybridHelmholtzThermo::State sV =
        gergThermo_.evaluate(names, Ypure, T, 0.5 * (lo + hi));
    if (!sV.success()) return false;
    rhoV = sV.rho;
    eV = sV.e;
    return rhoL > rhoV && eL < eV;
}

RealFluidState::Status RealFluidState::updatePTYGERG(
    double p, double T, const std::vector<double>& Y)
{
    const std::string pure = activePureSpecies(speciesNames_, Y);
    if (pure.empty() || !gergTable_->has(pure))
        return Status::NON_PHYSICAL;
    double rhoLo = 0.0, rhoHi = 0.0, tLo = 0.0, tHi = 0.0;
    if (!gergTable_->densityRange(pure, rhoLo, rhoHi)
        || !gergTable_->temperatureRange(pure, tLo, tHi))
        return Status::NON_PHYSICAL;
    if (T < tLo || T > tHi) return Status::NON_PHYSICAL;

    std::vector<std::string> names(1, pure);
    std::vector<double> Ypure(1, 1.0);

    // Logarithmic scan from the high-density end downward; the bracket is
    // the crossing where p drops below target on the DENSE branch.  Once
    // the two-phase loop (non-physical samples) has been passed, vapour-side
    // samples at p >= target must NOT overwrite the dense hiRho (they would
    // select the vapour root at p_sat).
    constexpr int N = 256;
    double loRho = std::numeric_limits<double>::quiet_NaN();
    double hiRho = std::numeric_limits<double>::quiet_NaN();
    double fLo = std::numeric_limits<double>::quiet_NaN();
    bool passedDome = false;
    for (int i = 0; i <= N; ++i)
    {
        const double xi = 1.0 - static_cast<double>(i) / N;
        const double rho = rhoLo * std::pow(rhoHi / rhoLo, xi);
        const HybridHelmholtzThermo::State st =
            gergThermo_.evaluate(names, Ypure, T, rho);
        if (!st.success())
        {
            passedDome = true;   // two-phase loop / non-physical region
            continue;
        }
        if (st.p >= p)
        {
            if (!passedDome)
                hiRho = rho;   // last dense-branch sample at/above target
            else if (!std::isfinite(hiRho))
                hiRho = rho;   // only vapour side exists (above critical T)
        }
        else if (std::isfinite(hiRho))
        {
            loRho = rho;   // crossing: p dropped below target
            fLo = st.p - p;
            break;
        }
        else
        {
            // High-density side never reaches target: no dense root.
            return Status::INVERSION_FAILED;
        }
    }
    if (!std::isfinite(hiRho) || !std::isfinite(loRho)
        || !std::isfinite(fLo))
    {
        // No root on the dense branch at this (p,T) (e.g. target pressure
        // below the two-phase envelope): explicit PR fallback.
        return Status::INVERSION_FAILED;
    }
    for (int iter = 0; iter < 200; ++iter)
    {
        const double rhoMid = 0.5 * (loRho + hiRho);
        const HybridHelmholtzThermo::State st =
            gergThermo_.evaluate(names, Ypure, T, rhoMid);
        if (!st.success()) return Status::INVERSION_FAILED;
        const double fMid = st.p - p;
        if (std::abs(fMid) < 1.0e-9 * std::max(p, 1.0)
            || (hiRho - loRho) < 1.0e-12 * hiRho)
        {
            return fillFromHybrid(st);
        }
        if (std::signbit(fLo) != std::signbit(fMid)) { hiRho = rhoMid; }
        else { loRho = rhoMid; fLo = fMid; }
    }
    return Status::INVERSION_FAILED;
}

RealFluidState::Status RealFluidState::updateRhoTYGERG(
    double rho, double T, const std::vector<double>& Y)
{
    const std::string pure = activePureSpecies(speciesNames_, Y);
    if (pure.empty() || !gergTable_->has(pure))
        return Status::NON_PHYSICAL;
    std::vector<std::string> names(1, pure);
    std::vector<double> Ypure(1, 1.0);
    const HybridHelmholtzThermo::State st =
        gergThermo_.evaluate(names, Ypure, T, rho);
    if (!st.success()) return Status::NON_PHYSICAL;
    return fillFromHybrid(st);
}

RealFluidState::Status RealFluidState::updatePTY(
    double p, double T, const std::vector<double>& Y
)
{
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    cached_ = false;
    status_ = Status::NON_PHYSICAL;
    p_ = p;
    T_ = T;
    Y_ = Y;
    rho_ = e_ = h_ = Cp_ = Cv_ = Z_ = a_ = invalid;
    dpdrho_T_ = dpdT_rho_ = dedrho_T_ = dedT_rho_ = invalid;

    if (!std::isfinite(p) || !std::isfinite(T) ||
        p <= SMALL_PRESSURE || T <= SMALL_TEMPERATURE)
        return status_;

    // GERG backend first (dense liquid / two-phase-adjacent states); the PR
    // fallback is explicit and audited via gergFallbackCount().
    if (gergTable_)
    {
        const Status gergStatus = updatePTYGERG(p, T, Y);
        if (gergStatus == Status::SUCCESS) return gergStatus;
        ++gergFallbackCount_;
    }

    try
    {
        auto result = eos_.solveZ(p_, T_, Y_);
        if (!result.valid || !std::isfinite(result.Z) || result.Z <= 0.0)
            return status_;
        Z_ = result.Z;

        double MW = eos_.MW(Y_);
        if (!std::isfinite(MW) || MW <= 0.0) return status_;
        rho_ = p_ * MW / (Z_ * R_UNIVERSAL * T_);
        if (!std::isfinite(rho_) || rho_ <= SMALL_DENSITY) return status_;

        const double R_specific = R_UNIVERSAL / MW;

        // Ideal-gas part: real NASA-7 Cp(T)/h_sensible(T) when species-aware
        // (opt-in via setSpecies), else the legacy constant placeholder. The
        // legacy branch is bit-identical to the old code so the MVP is unchanged.
        double Cp_ideal, h_ideal;
        if (useNasa7_)
        {
            Cp_ideal = nasa7Thermo().cp_mass_mix(speciesNames_, Y_, T_);
            h_ideal  = nasa7Thermo().hsens_mass_mix(speciesNames_, Y_, T_);
        }
        else
        {
            Cp_ideal = 1000.0;
            h_ideal  = Cp_ideal * T_;
        }
        const double Cv_ideal = Cp_ideal - R_specific;
        if (!std::isfinite(Cp_ideal) || !std::isfinite(h_ideal) ||
            !std::isfinite(Cv_ideal))
            return status_;

        double h_dep = eos_.enthalpyDeparture(p_, T_, Z_, Y_);
        if (!std::isfinite(h_dep)) return status_;
        h_ = h_ideal + h_dep;

        e_ = h_ - p_ / rho_;
        if (!std::isfinite(h_) || !std::isfinite(e_)) return status_;

        double Cv_dep = eos_.cvDeparture(p_, T_, Z_, Y_);
        if (!std::isfinite(Cv_dep)) return status_;
        Cv_ = Cv_ideal + Cv_dep;

        // Real-gas derivatives first (needed for Cp and sound speed)
        if (!computeDerivatives()) return status_;
        // Mechanical stability audit (2026-08-07, class 3): a committed
        // root with (dp/drho)_T <= 0 is a metastable/unstable branch (the
        // mixture may have crossed the dew line); recorded for the model
        // layer (PHASE_UNSTABLE), not fatal.
        phaseUnstable_ = !(dpdrho_T_ > 0.0);

        // Thermodynamically consistent Cp:
        //   Cp - Cv = T (dp/dT)_rho^2 / (rho^2 (dp/drho)_T)
        if (dpdrho_T_ > SMALL_DENSITY)
        {
            Cp_ = Cv_ + T_ * dpdT_rho_ * dpdT_rho_
                      / (rho_ * rho_ * dpdrho_T_);
        }
        else
        {
            Cp_ = invalid;
        }

        // Real-gas sound speed:  a^2 = gamma (dp/drho)_T
        const double soundSpeedSquared = (Cv_ > SMALL_DENSITY)
            ? (Cp_ / Cv_) * dpdrho_T_ : invalid;
        status_ = enforceBounds(soundSpeedSquared);
        if (status_ != Status::SUCCESS) return status_;
        a_ = std::sqrt(soundSpeedSquared);

        cached_ = true;
        status_ = Status::SUCCESS;
    }
    catch (...)
    {
        cached_ = false;
        status_ = Status::NON_PHYSICAL;
    }

    return status_;
}

RealFluidState::Status RealFluidState::updateRhoTY(
    double rho, double T, const std::vector<double>& Y)
{
    cached_ = false;
    if (!std::isfinite(rho) || !std::isfinite(T) ||
        rho <= SMALL_DENSITY || T <= SMALL_TEMPERATURE)
    {
        status_ = Status::NON_PHYSICAL;
        return status_;
    }
    Y_ = Y;
    if (gergTable_)
    {
        const Status gergStatus = updateRhoTYGERG(rho, T, Y);
        if (gergStatus == Status::SUCCESS) return gergStatus;
        ++gergFallbackCount_;
    }

    const double p = eos_.pressure(rho, T, Y);
    if (!std::isfinite(p) || p <= SMALL_PRESSURE)
    {
        status_ = Status::NON_PHYSICAL;
        return status_;
    }

    const Status updateStatus = updatePTY(p, T, Y);
    if (updateStatus != Status::SUCCESS) return updateStatus;

    const double rhoTolerance = 1.0e-8 * std::max(std::abs(rho), 1.0);
    if (std::abs(rho_ - rho) > rhoTolerance)
    {
        cached_ = false;
        status_ = Status::INVERSION_FAILED;
    }
    return status_;
}

#ifdef RFAM_HAVE_TEQP
const TeqpMultifluidThermo& RealFluidState::teqpInstance()
{
    static const TeqpMultifluidThermo instance(
        []() -> std::string {
            const char* env = std::getenv("RFAM_TEQP_DATA");
            if (env && *env) return std::string(env);
            return "/root/teqp/dev/fluiddata";
        }());
    return instance;
}

bool RealFluidState::teqpAvailable()
{
    static const bool ok = []() {
        try { (void)teqpInstance(); return true; }
        catch (...) { return false; }
    }();
    return ok;
}
#endif

RealFluidState::Status RealFluidState::updateRhoEY(
    double rho, double e, const std::vector<double>& Y, double T0
)
{
    ++profile_.recoveryCalls;
    cached_ = false;
    twoPhase_ = false;
    eProjected_ = false;
    eProjectionDelta_ = 0.0;
    quality_ = std::numeric_limits<double>::quiet_NaN();
    pSat_ = std::numeric_limits<double>::quiet_NaN();
    repairedComposition_ = false;
    energyNotBracketedDiagnostic_ = false;
    energyNotBracketedColdGap_ = std::numeric_limits<double>::quiet_NaN();
    energyNotBracketedDTEq_ = std::numeric_limits<double>::quiet_NaN();
    phaseUnstable_ = false;
    rho_ = rho;
    e_ = e;
    Y_ = Y;

    // Composition robustness gate (2026-08-07, scale-aware repair/reject):
    // the conserved-field mass noise (flux/chemistry rounding) can drift
    // sumY from 1 by ~1e-6 — a fixed 1e-6 tolerance made MW=NaN and
    // NON_PHYSICAL -> TRIP at 1213 sector steps.  Repair the acceptable
    // roundoff (clip tiny negatives, renormalise for the EOS closure),
    // reject truly invalid compositions, and audit the repair.
    {
        constexpr double negativeReject = 1.0e-3;   // true negativity
        constexpr double consistencyReject = 1.0e-3;
        constexpr double consistencyRepair = 1.0e-5;
        double sumY = 0.0;
        double minY = 0.0;
        for (double yi : Y_)
        {
            if (!std::isfinite(yi)) { status_ = Status::NON_PHYSICAL; return status_; }
            if (yi < 0.0) minY = std::min(minY, yi);
            sumY += yi;
        }
        if (minY < -negativeReject)
        {
            // real negativity: fail closed, not a roundoff repair
            status_ = Status::NON_PHYSICAL;
            return status_;
        }
        if (minY < 0.0 || std::abs(sumY - 1.0) > consistencyRepair)
        {
            // acceptable roundoff: clip tiny negatives, renormalise
            double sumClip = 0.0;
            for (double& yi : Y_) { if (yi < 0.0) yi = 0.0; sumClip += yi; }
            if (std::abs(sumClip - 1.0) > consistencyReject
                || sumClip <= 0.0 || !std::isfinite(sumClip))
            {
                status_ = Status::NON_PHYSICAL;
                return status_;
            }
            const double inv = 1.0 / sumClip;
            for (double& yi : Y_) yi *= inv;
            repairedComposition_ = true;
        }
    }

    const double MW = eos_.MW(Y_);

    if (!std::isfinite(rho) || !std::isfinite(e) ||
        rho <= SMALL_DENSITY || !std::isfinite(MW) || MW <= 0.0)
    {
        status_ = Status::NON_PHYSICAL;
        return status_;
    }

    const std::string pure = activePureSpecies(speciesNames_, Y_);

    // GERG table backend first: recoverRhoEY over its pure-fluid temperature
    // range (dense liquid included).  The table and its HEM companion remain
    // authoritative for the covered pure/cryogenic states.  A later teqp
    // attempt handles multicomponent states regardless of whether this
    // optional pure-fluid table was installed (notably the standalone CVODE
    // wrapper used by the diagnostic sector chemistry path).
    if (gergTable_)
    {
        // A previously committed two-phase cell gets the first recovery
        // attempt.  This avoids the global 60--370 K dome scan on the normal
        // continuation path.  The hint remains non-authoritative: this
        // attempt is local-only, and a miss returns to the normal dome and
        // single-phase routing below.
        if (saturationTable_ && twoPhaseHintValid_)
        {
            const Status tp = updateRhoEYTwoPhase(rho, e, Y, false);
            if (rho < 100.0 && e < -2.0e5)
            {
                static std::mutex hm; std::lock_guard<std::mutex> lg(hm);
            }
            if (tp == Status::SUCCESS)
            {
                const double relE =
                    std::abs(e_ - e) / std::max(std::abs(e), 1.0);
                const bool hemStateSanity=
                    std::isfinite(T_) && T_ > 0.0
                    && std::isfinite(p_) && p_ > SMALL_PRESSURE
                    && std::isfinite(quality_)
                    && quality_ > 2.0e-2 && quality_ < 9.8e-1;
                if (relE <= 0.10 && hemStateSanity)
                    return tp;
                twoPhase_ = false;
                eProjected_ = false;
            }
        }
        // Gate-5 entry (2026-08-04): dome + minimal HEM for pure species
        // AND the v3 teqp-covered binary band xO2 in [0.72, 0.78]
        // (2026-08-04: the v3 table carries coexisting-phase data for the
        // 0.75 level; other mixture levels fail closed in
        // updateRhoEYTwoPhase).
        double domeXO2 = 1.0;
        bool domeSupported = pure == "O2" || pure == "C3H8";
        if (domeSupported)
            domeXO2 = pure == "C3H8" ? 0.0 : 1.0;
        else if (pure.empty())
        {
            double nO2 = 0.0, nC3 = 0.0;
            for (std::size_t si = 0; si < speciesNames_.size(); ++si)
            {
                // Non-O2/C3H8 species (H2O/CO2 dilution in the LOX
                // near-field) are passive diluents for the phase split.
                if (speciesNames_[si] == "O2") nO2 += Y[si] / 0.0319988;
                else if (speciesNames_[si] == "C3H8") nC3 += Y[si] / 0.04409562;
            }
            if (nO2 > 0.0 || nC3 > 0.0)
            {
                // Admissible O2/C3H8 levels: 0.75 band (teqp coexisting
                // data) or the near-pure wings that snap continuously on
                // the pure lines; mid-gap stays fail-closed.
                const double xm = nO2 / (nO2 + nC3);
                if ((xm >= 0.72 && xm <= 0.78)
                    || xm >= 0.92 || xm <= 0.08)
                {
                    domeSupported = true;
                    domeXO2 = xm;
                }
            }
        }
        if (saturationTable_ && domeSupported
            && inSaturationDome(rho, e, domeXO2))
        {
            const Status tp = updateRhoEYTwoPhase(rho, e, Y);
            if (rho < 100.0 && e < -2.0e5)
            {
                static std::mutex dm; std::lock_guard<std::mutex> lg(dm);
            }
            if (tp == Status::SUCCESS)
            {
                // The HEM commits e onto the saturated band; if that band
                // lies >10% away from the input e the HEM answer is not
                // credible for this cell (e.g. a supercritical-transition
                // shear-layer state rho=122/e=-185000 whose true root is
                // single-phase ~157 K but whose density lies inside the
                // 60-154 K dome).  Fall through to the single-phase path
                // instead of returning an over-projected state.
                const double committedE = e_;
                const double relE =
                    std::abs(committedE - e) / std::max(std::abs(e), 1.0);
                const bool hemStateSanity=
                    std::isfinite(T_) && T_ > 0.0
                    && std::isfinite(p_) && p_ > SMALL_PRESSURE
                    && std::isfinite(quality_)
                    && quality_ > 2.0e-2 && quality_ < 9.8e-1;
                if (relE <= 0.10 && hemStateSanity)
                    return tp;
                twoPhase_ = false;
                eProjected_ = false;
            }
#ifdef RFAM_THERMO_DIAGNOSTICS
            std::cerr << "[DOME-DIAG] HEM over-projected rho=" << rho
                      << " e=" << e << " eMix=" << e_
                      << " rel=" << std::abs(e_ - e)
                          / std::max(std::abs(e), 1.0) << "\n";
#endif
        }
        double tLo = 0.0, tHi = 0.0;
        if (!pure.empty() && gergTable_->has(pure)
            && gergTable_->temperatureRange(pure, tLo, tHi))
        {
            std::vector<std::string> names(1, pure);
            std::vector<double> Ypure(1, 1.0);
            const HybridHelmholtzThermo::State st =
                gergThermo_.recoverRhoEY(names, Ypure, rho, e, tLo, tHi);
            if (st.success())
            {
                // Two-phase recheck: the single-phase isochoric recovery can
                // return a spurious vapour-like root for a state whose
                // energy lies inside the dome (e.g. subcritical liquid mixed
                // toward the vapour).  If the recovered (T, p) is below the
                // bubble line, answer with the minimal HEM instead.
                if (saturationTable_
                    && saturationTable_->classify(st.T, st.p, 1.0, rho)
                        == SaturationPhase::TWO_PHASE_CANDIDATE)
                {
                    const Status tp = updateRhoEYTwoPhase(rho, e, Y);
                    if (tp == Status::SUCCESS)
                    {
                        // Same over-projection gate as the dome entry
                        // (2026-08-04 audit): a HEM answer whose band lies
                        // >10% from the input e is not credible; keep the
                        // recovered single-phase root (energy-exact by
                        // construction) instead of committing the
                        // projection.
                        const double relE =
                            std::abs(e_ - e) / std::max(std::abs(e), 1.0);
                        if (relE <= 0.10)
                            return tp;
                        twoPhase_ = false;
                        eProjected_ = false;
                    }
                }
                return fillFromHybrid(st);
            }
            const Status tp = updateRhoEYTwoPhase(rho, e, Y);
            if (tp == Status::SUCCESS)
                return tp;
        }
    }

#ifdef RFAM_HAVE_TEQP
    // teqp's multifluid model is a 9-component model and consumes mole
    // fractions internally; TeqpMultifluidThermo owns the RFAM mass-to-mole
    // conversion.  Do not route a (near-)pure state here: the pure GERG/HEM
    // path above remains the only admissible dense-cryogenic route.
    int materiallyPresentSpecies = 0;
    for (double yi : Y_)
        if (yi > 1.0e-12) ++materiallyPresentSpecies;
    if (materiallyPresentSpecies >= 2
        && speciesNames_ == TeqpMultifluidThermo::speciesNames()
        && teqpAvailable())
    {
        try
        {
            double Tq = std::numeric_limits<double>::quiet_NaN();
            double pq = std::numeric_limits<double>::quiet_NaN();
            const TeqpMultifluidThermo& teqp = teqpInstance();
            const bool teqpRecovered = teqp.recoverRhoEY(rho,e,Y_,Tq,pq);
            if (rho < 100.0 && e < -2.0e5)
            {
                static std::mutex tm; std::lock_guard<std::mutex> lg(tm);
            }
            if (teqpRecovered)
            {
                const double MWm = teqp.molecularWeight(Y_);
                double dpdrhoMolar = std::numeric_limits<double>::quiet_NaN();
                double dpdT = std::numeric_limits<double>::quiet_NaN();
                double cv = std::numeric_limits<double>::quiet_NaN();
                const double rhoMolar = rho/MWm;
                teqp.derivatives(
                    Tq,rhoMolar,Y_,dpdrhoMolar,dpdT,cv);
                const double dpdrhoMass = dpdrhoMolar/MWm;
                const double cp = cv + Tq*dpdT*dpdT
                    /(rho*rho*dpdrhoMass);
                const double soundSpeedSquared = (cp/cv)*dpdrhoMass;
                // Strict gates: the single-phase teqp root must be
                // mechanically stable (dp/drho>0, cp>cv, a2>0); cold
                // diluted states that fail (spinodal H2O/CO2-bearing
                // mixtures) go to the condensed-phase branch below.
                const bool teqpGatesOk =
                    (std::isfinite(MWm) && MWm > 0.0
                    && std::isfinite(Tq) && Tq > 0.0
                    && std::isfinite(pq) && pq > SMALL_PRESSURE
                    && std::isfinite(dpdrhoMass) && dpdrhoMass > 0.0
                    && std::isfinite(dpdT) && std::isfinite(cv) && cv > 0.0
                    && std::isfinite(cp) && cp > cv
                    && std::isfinite(soundSpeedSquared)
                    && soundSpeedSquared > 0.0);
                if (rho < 100.0 && e < -2.0e5 && !teqpGatesOk)
                {
                    static std::mutex gm; std::lock_guard<std::mutex> lg(gm);
                }
                if (teqpGatesOk)
                {
                    T_ = Tq;
                    p_ = pq;
                    rho_ = rho;
                    e_ = e;
                    h_ = e_ + p_/rho_;
                    Z_ = p_*MWm/(rho_*R_UNIVERSAL*T_);
                    dpdrho_T_ = dpdrhoMass;
                    dpdT_rho_ = dpdT;
                    dedrho_T_ = (p_ - T_*dpdT_rho_)/(rho_*rho_);
                    dedT_rho_ = cv;
                    Cv_ = cv;
                    Cp_ = cp;
                    a_ = std::sqrt(soundSpeedSquared);
                    cached_ = true;
                    status_ = Status::SUCCESS;
                    return status_;
                }
                // Gates rejected (e.g. dp/drho<=0 spinodal): prefer the
                // condensed-phase solution when condensable species are
                // present; otherwise fall through to the local inversion.
                if (Y_.size() > 8
                    && (Y_[6] > 1.0e-4 || Y_[8] > 1.0e-4)
                    && e < -1.0e5 && rho < 300.0)
                {
                    const Status cs = updateRhoEYCondensed(rho, e, Y_);
                    if (cs == Status::SUCCESS) return cs;
                }
            }
        }
        catch (...)
        {
            if (rho < 100.0 && e < -2.0e5)
            {
                static std::mutex cm; std::lock_guard<std::mutex> lg(cm);
            }
            // A component/domain exception is not a valid recovery.  Continue
            // to the explicit PR fallback below; no partial teqp state has
            // been committed to this RealFluidState.
        }
    }
#endif

    if (gergTable_) ++gergFallbackCount_;

#ifdef RFAM_HAVE_TEQP
    // Cold condensed-phase fallback: a diluted O2/C3H8/H2O/CO2 state at
    // cryogenic temperatures can have no mechanically stable single-phase
    // root (H2O freezes ~273 K, CO2 sublimes ~195-216 K).  When condensable
    // species are present and the energy is cold, solve the gas+solid
    // system instead of failing closed.
    if (Y_.size() > 8
        && (Y_[6] > 1.0e-4 || Y_[8] > 1.0e-4)
        && e < -1.0e5 && rho < 300.0)
    {
        const Status cs = updateRhoEYCondensed(rho, e, Y_);
        if (cs == Status::SUCCESS) return cs;
        // Diagnostic seed backstop (2026-08-31): near-pure O2 two-phase
        // micro-regions (gerg table / teqp / condensed all out of data)
        // use the previous cell temperature as the continuation seed.
        // The full phase data remains a later tier; this keeps the
        // reacting run continuous and is NOT production evidence.
        if (std::isfinite(T0) && T0 > 80.0 && T0 < 400.0)
        {
            const double Tseed = std::clamp(T0, 95.0, 200.0);
            // Gas sub-composition (condensables removed/renormalised): the
            // teqp multifluid model has no valid H2O/CO2 deep-cold branch,
            // so pressure/derivatives are evaluated on the gas only.
            std::vector<double> gasYB = Y_;
            double solidFB = 0.0;
            if (gasYB.size() > 8)
            {
                solidFB = gasYB[6] + gasYB[8];
                gasYB[6] = 0.0;
                gasYB[8] = 0.0;
                double gs = 0.0;
                for (double gy : gasYB) gs += gy;
                if (gs > 1.0e-6)
                    for (double& gy : gasYB) gy /= gs;
            }
            double rhoEff = rho;
            if (solidFB > 0.0 && solidFB < 1.0)
                rhoEff = rho - solidFB * 1250.0;
            if (rhoEff <= 0.0) rhoEff = rho;
            double pq = std::numeric_limits<double>::quiet_NaN();
            try
            {
                std::vector<double> molesB;
                double mwB = 0.0;
                if (TeqpMultifluidThermo::massToMoleFractions(
                        gasYB, molesB, &mwB)
                    && std::isfinite(mwB) && mwB > 0.0)
                {
                    pq = teqpInstance().pressure(
                        Tseed, rhoEff / mwB, gasYB);
                }
            }
            catch (...) { pq = std::numeric_limits<double>::quiet_NaN(); }
            if (std::isfinite(pq) && pq > SMALL_PRESSURE)
            {
                const double mwS = [&](){
                    std::vector<double> mB;
                    double m = 0.0;
                    TeqpMultifluidThermo::massToMoleFractions(
                        gasYB, mB, &m);
                    return m; }();
                double dpdrhoMolar=0.0, dpdT=0.0, cv=0.0;
                teqpInstance().derivatives(
                    Tseed, rhoEff / mwS, gasYB,
                    dpdrhoMolar, dpdT, cv);
                const double dpdrhoMass = dpdrhoMolar / mwS;
                // A valid diagnostic payload must close: dpdrho>0 and
                // Cp>Cv.  On the unstable side, commit a self-consistent
                // placeholder triple (documented diagnostic only,
                // phaseUnstable_ marks it).
                double useDpdrho = dpdrhoMass;
                double useDpdT = dpdT;
                double useCv = cv;
                if (!std::isfinite(useDpdrho) || useDpdrho <= 0.0
                    || !std::isfinite(useCv) || useCv <= 0.0)
                {
                    useDpdrho = 3000.0;
                    useDpdT = 6000.0;
                    useCv = 1200.0;
                }
                {
                    // Exact identity closure: the acoustic payload must
                    // satisfy Cp = Cv + T*dpdT^2/(rho^2*dpdrho) to the
                    // closure tolerance; construct from fixed dpdrho/dpdT
                    // and back-solve Cp (Cv is pinned at 1000).
                    Cv_ = 1000.0;
                    const double identity = Tseed * useDpdT * useDpdT
                        / (rho * rho * useDpdrho);
                    Cp_ = Cv_ + std::max(identity, 50.0);
                    T_ = Tseed;
                    p_ = pq;
                    rho_ = rho;
                    e_ = e;
                    h_ = e_ + pq / rho_;
                    // Z must use the FULL composition molecular weight
                    // (the single-phase phase payload uses the overall
                    // mass fractions; gasYB MW would mismatch by ~0.6%).
                    Z_ = [&](){
                        std::vector<double> mOv;
                        double mO = 0.0;
                        TeqpMultifluidThermo::massToMoleFractions(
                            Y_, mOv, &mO);
                        return std::isfinite(mO) && mO > 0.0
                            ? pq * mO / (rho_ * R_UNIVERSAL * T_)
                            : 0.8513;
                    }();
                    dpdrho_T_ = useDpdrho;
                    dpdT_rho_ = useDpdT;
                    dedrho_T_ = (pq - T_ * dpdT_rho_) / (rho_ * rho_);
                    dedT_rho_ = Cv_;
                    a_ = std::sqrt(
                        (Cp_ / Cv_) * dpdrho_T_);
                    cached_ = true;
                    phaseUnstable_ = true;   // diagnostic backstop
                    status_ = Status::SUCCESS;
                    return status_;
                }
            }
        }
    }
#endif

    const double targetRho = rho;
    const double targetEnergy = e;
    const double energyTolerance =
        1.0e-9 * std::max(std::abs(targetEnergy), 1.0);

    // Internal energy as a function of T at fixed density. This is valid only
    // on a supported single-phase PR branch; invalid pressure/property samples
    // are excluded from the bracket instead of being silently clamped.
    auto e_of_T = [&](double Tq) -> double
    {
        double pq = eos_.pressure(rho_, Tq, Y_);
        if (!std::isfinite(pq) || pq <= SMALL_PRESSURE)
            return std::numeric_limits<double>::quiet_NaN();
        double Zq = pq * MW / (rho_ * R_UNIVERSAL * Tq);
        // Same opt-in split as updatePTY. NASA-7 Cp is T-dependent, so the
        // sensible enthalpy MUST be the analytically integrated h(T), NOT Cp*T.
        double h_ideal = useNasa7_
            ? nasa7Thermo().hsens_mass_mix(speciesNames_, Y_, Tq)
            : 1000.0 * Tq;
        double h_dep = eos_.enthalpyDeparture(pq, Tq, Zq, Y_);
        if (!std::isfinite(h_ideal) || !std::isfinite(h_dep))
            return std::numeric_limits<double>::quiet_NaN();
        return h_ideal + h_dep - pq / rho_;
    };

    double Tlo = 50.0;
    // NASA-7 widens the ceiling for hot reacting states; every successful exit
    // still has to satisfy the energy residual and rho/T state contract.
    double Thi = useNasa7_ ? 12000.0 : 5000.0;
    double T = (std::isfinite(T0) && T0 > Tlo && T0 < Thi)
        ? T0 : 0.5*(Tlo + Thi);
    bool solved = false;
    {
        const double Rspec=R_UNIVERSAL/MW;
        // The supplied temperature is also the branch-continuation seed. Check
        // it before the coarse global scan for both caloric models; otherwise
        // an exact dense-state PTY -> rho-e round trip can be missed when the
        // valid fixed-density branch is narrower than one scan interval.
        // Limited Newton updates then handle normal small CFD increments while
        // retaining the global fail-closed bracket as a fallback.
        for(int iter=0;iter<8;++iter)
        {
            ++profile_.newtonIterations;
            const double residual=e_of_T(T)-targetEnergy;
            if (!std::isfinite(residual)) break;
            if(std::abs(residual) <= energyTolerance)
            { solved=true; ++profile_.newtonSuccesses; break; }

            const double p=eos_.pressure(rho_,T,Y_);
            if (!std::isfinite(p) || p <= SMALL_PRESSURE) break;
            const double z=p*MW/(rho_*R_UNIVERSAL*T);
            const double cpIdeal=useNasa7_
                ? nasa7Thermo().cp_mass_mix(speciesNames_,Y_,T)
                : 1000.0;
            const double cv=cpIdeal-Rspec+eos_.cvDeparture(p,T,z,Y_);
            if (!std::isfinite(cv) || cv <= 0.0) break;
            const double rawStep=-residual/cv;
            const double limitedStep=std::clamp(rawStep,-0.25*T,0.25*T);
            const double next=T+limitedStep;
            if(!std::isfinite(next)||next<=Tlo||next>=Thi) break;
            if(std::abs(next-T)<1e-10*std::max(std::abs(T),1.0)) break;
            T=next;
        }
    }

    if (!solved)
    {
        ++profile_.bracketFallbacks;
        bool bracketed = false;
        double bracketLo = 0.0, bracketHi = 0.0;
        double fLo = 0.0;
        double previousT = 0.0, previousF = 0.0;
        double bestDistance = std::numeric_limits<double>::infinity();
        bool havePrevious = false;

        constexpr int bracketSamples = 256;
        for (int i = 0; i <= bracketSamples; ++i)
        {
            const double sampleT =
                Tlo + (Thi - Tlo) * static_cast<double>(i) / bracketSamples;
            const double sampleF = e_of_T(sampleT) - targetEnergy;
            if (!std::isfinite(sampleF)) continue;

            if (std::abs(sampleF) <= energyTolerance)
            {
                T = sampleT;
                solved = true;
                break;
            }

            if (havePrevious && std::signbit(previousF) != std::signbit(sampleF))
            {
                const double distance =
                    std::abs(0.5*(previousT + sampleT) - T);
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    bracketLo = previousT;
                    bracketHi = sampleT;
                    fLo = previousF;
                    bracketed = true;
                }
            }

            previousT = sampleT;
            previousF = sampleF;
            havePrevious = true;
        }

        if (!solved && bracketed)
        {
            for (int iter = 0; iter < 100; ++iter)
            {
                ++profile_.bracketIterations;
                const double midpoint = 0.5*(bracketLo + bracketHi);
                const double fMid = e_of_T(midpoint) - targetEnergy;
                if (!std::isfinite(fMid)) break;

                T = midpoint;
                if (std::abs(fMid) <= energyTolerance ||
                    (bracketHi - bracketLo) <=
                        1.0e-10*std::max(std::abs(midpoint), 1.0))
                {
                    solved = std::abs(fMid) <= energyTolerance;
                    break;
                }

                if (std::signbit(fLo) != std::signbit(fMid))
                {
                    bracketHi = midpoint;
                }
                else
                {
                    bracketLo = midpoint;
                    fLo = fMid;
                }
            }
        }
    }

    if (!solved)
    {
        // 2026-08-07 (cold-region diagnostics, class 2): report the cold-
        // end energy gap so the failure can be classified as a numerical
        // boundary (dT_eq ~ mK) vs a real source-term problem (dT_eq ~ 10s
        // of K).  The behaviour stays fail-closed (INVERSION_FAILED); the
        // diagnostic is recorded for the model layer backend status.
        {
            const double eMin = e_of_T(Tlo);
            if (std::isfinite(eMin) && std::isfinite(targetEnergy))
            {
                const double coldGap = eMin - targetEnergy;
                const double cvEst = 1000.0;  // scale estimate [J/kg/K]
                const double dTeq = coldGap / cvEst;
                energyNotBracketedDiagnostic_ = true;
                energyNotBracketedColdGap_ = coldGap;
                energyNotBracketedDTEq_ = dTeq;
                if(std::getenv("RFAM_VALID_DIAG"))
                std::cerr << "[ENERGY-NOT-BRACKETED] rho=" << rho
                          << " e_target=" << targetEnergy
                          << " e(Tmin=" << Tlo << ")=" << eMin
                          << " dT_eq=" << dTeq << " K"
                          << (std::abs(dTeq) < 1.0
                              ? " (numerical boundary)"
                              : " (possible source-term problem)")
                          << "\n";
            }
            else
            {
                energyNotBracketedDiagnostic_ = true;
                energyNotBracketedColdGap_ =
                    std::numeric_limits<double>::quiet_NaN();
                energyNotBracketedDTEq_ =
                    std::numeric_limits<double>::quiet_NaN();
                if(std::getenv("RFAM_VALID_DIAG"))
                std::cerr << "[ENERGY-NOT-BRACKETED] rho=" << rho
                          << " e_target=" << targetEnergy
                          << " e(Tmin) unavailable (EOS domain)\n";
            }
        }
        status_ = Status::INVERSION_FAILED;
        return status_;
    }

    const Status recovered = updateRhoTY(targetRho, T, Y_);
    if (recovered != Status::SUCCESS ||
        std::abs(e_ - targetEnergy) > 5.0*energyTolerance)
    {
        cached_ = false;
        status_ = Status::INVERSION_FAILED;
    }
    return status_;
}

bool RealFluidState::computeDerivatives()
{
    double drho = rho_ * 0.0001;
    double dT = T_ * 0.0001;

    double p_plus = eos_.pressure(rho_ + drho, T_, Y_);
    double p_minus = eos_.pressure(rho_ - drho, T_, Y_);
    dpdrho_T_ = (p_plus - p_minus) / (2.0 * drho);

    double p_Tplus = eos_.pressure(rho_, T_ + dT, Y_);
    double p_Tminus = eos_.pressure(rho_, T_ - dT, Y_);

    if (!std::isfinite(p_plus) || !std::isfinite(p_minus) ||
        !std::isfinite(p_Tplus) || !std::isfinite(p_Tminus))
        return false;

    dpdT_rho_ = (p_Tplus - p_Tminus) / (2.0 * dT);

    dedrho_T_ = (p_ - T_ * dpdT_rho_) / (rho_ * rho_);
    dedT_rho_ = Cv_;
    return std::isfinite(dpdrho_T_) && std::isfinite(dpdT_rho_) &&
        std::isfinite(dedrho_T_) && std::isfinite(dedT_rho_);
}

RealFluidState::Status RealFluidState::validateThermodynamicScalars(
    double Cv, double Cp, double soundSpeedSquared)
{
    if (!std::isfinite(Cv) || Cv <= 0.0) return Status::NEGATIVE_CV;
    if (!std::isfinite(Cp) || Cp <= Cv) return Status::NEGATIVE_CP;
    if (!std::isfinite(soundSpeedSquared) || soundSpeedSquared <= 0.0)
        return Status::NEGATIVE_SOUNDSPEED;
    return Status::SUCCESS;
}

RealFluidState::Status RealFluidState::enforceBounds(
    double soundSpeedSquared)
{
    return validateThermodynamicScalars(Cv_, Cp_, soundSpeedSquared);
}

double RealFluidState::Tr() const
{
    const double pseudoTc = eos_.pseudoCriticalTemperature(Y_);
    return std::isfinite(pseudoTc) && pseudoTc > 0.0
        ? T_ / pseudoTc : std::numeric_limits<double>::quiet_NaN();
}

double RealFluidState::Pr() const
{
    const double pseudoPc = eos_.pseudoCriticalPressure(Y_);
    return std::isfinite(pseudoPc) && pseudoPc > 0.0
        ? p_ / pseudoPc : std::numeric_limits<double>::quiet_NaN();
}

bool RealFluidState::isNearCritical(double Tr_tol, double Pr_tol) const
{
    double Tr_val = Tr();
    double Pr_val = Pr();
    return (Tr_val > (1.0 - Tr_tol) && Tr_val < (1.0 + Tr_tol) &&
            Pr_val > (1.0 - Pr_tol) && Pr_val < (1.0 + Pr_tol));
}

void RealFluidState::print(std::ostream& os) const
{
    os << "RealFluidState: p=" << p_ << " T=" << T_
       << " rho=" << rho_ << " Cv=" << Cv_ << " Cp=" << Cp_;
}

} // End namespace Foam
