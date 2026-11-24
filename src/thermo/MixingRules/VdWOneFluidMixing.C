#include "VdWOneFluidMixing.H"
#include "PengRobinsonAlpha.H"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Foam
{

static constexpr double R_UNIVERSAL = 8.314462618;

static double invalidValue()
{
    return std::numeric_limits<double>::quiet_NaN();
}

VdWOneFluidMixing::VdWOneFluidMixing(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij
)
:
    criticalData_(criticalData),
    kij_(kij),
    multiSpeciesEnabled_(false)   // legacy default: pseudo-critical path, byte-for-byte MVP
{
    // Default 2-species set {O2, C3H8} preserves the legacy EOS path used by
    // the bit-for-bit MVP/precond reproductions. Route through setSpeciesNames
    // so the cache is the SINGLE source of truth (no separate default path).
    setSpeciesNames({"O2", "C3H8"});
    // NOTE: multiSpeciesEnabled_ stays false here. It is flipped to true only
    // by enableMultiSpeciesMixing() called from PengRobinsonEOS::setSpeciesNames
    // (the external opt-in path). This ensures the internal ctor call above
    // does NOT change the legacy EOS behaviour.
}

void VdWOneFluidMixing::setSpeciesNames(const std::vector<std::string>& names)
{
    // Explicit failure, no silent fallback: every name MUST resolve to a
    // physically-valid critical-property entry. A typo (e.g. "CO2" vs a
    // mis-encoded variant) throws here instead of silently poisoning MW/a/b
    // for the whole simulation. [[rfam-species-energy-foundation]]
    if (names.empty())
        throw std::runtime_error("VdWOneFluidMixing::setSpeciesNames: empty species list");

    std::vector<CriticalProperties> cache;
    cache.reserve(names.size());
    for (const std::string& nm : names)
    {
        if (!criticalData_.has(nm))
            throw std::runtime_error(
                "VdWOneFluidMixing::setSpeciesNames: unknown species '" + nm +
                "' (not in SpeciesCriticalDB)");
        const CriticalProperties& p = criticalData_.get(nm);
        if (!(p.MW > 0.0) || !(p.Tc > 0.0) || !(p.Pc > 0.0))
            throw std::runtime_error(
                "VdWOneFluidMixing::setSpeciesNames: non-physical critical "
                "properties for '" + nm + "' (MW/Tc/Pc must be > 0)");
        cache.push_back(p);   // omega may be negative (e.g. H2), do not gate on it
    }

    speciesNames_ = names;
    cachedProps_  = cache;
}

bool VdWOneFluidMixing::isValidComposition(
    const std::vector<double>& Y) const
{
    if (Y.size() != cachedProps_.size() || Y.empty()) return false;

    double sumY = 0.0;
    for (double y : Y)
    {
        if (!std::isfinite(y) || y < 0.0) return false;
        sumY += y;
    }

    // 2026-08-07: scale-aware tolerance (absTol + relTol*scale).  The old
    // fixed 1e-6 was hit by the conserved-field mass noise at 1213 sector
    // steps (sumY drift ~1.1e-6) -> MW=NaN -> NON_PHYSICAL -> TRIP.
    // Scale = max(|sumY|, 1); the 1e-4 relative component matches the
    // numerical-noise level (a 1e-4 composition error shifts MW by
    // ~3e-6 kg/mol = 0.01%, far below the state-commit tolerances).
    // The primary repair path is RealFluidState::updateRhoEY (clip +
    // renormalise with audit); this is the defensive backstop for any
    // direct EOS caller.
    const double scale = std::max(std::abs(sumY), 1.0);
    return std::isfinite(sumY)
        && std::abs(sumY - 1.0) <= 1.0e-10 + 1.0e-4 * scale;
}

double VdWOneFluidMixing::Tc(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    double Tc_mix = 0.0;
    double sum_x = 0.0;

    for (size_t i = 0; i < cachedProps_.size() && i < Y.size(); i++)
    {
        const CriticalProperties& props = cachedProps_[i];
        double x_i = Y[i] / props.MW;
        Tc_mix += x_i * props.Tc;
        sum_x += x_i;
    }

    return Tc_mix / (sum_x + 1e-20);
}

double VdWOneFluidMixing::Pc(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    double Pc_mix = 0.0;
    double sum_x = 0.0;

    for (size_t i = 0; i < cachedProps_.size() && i < Y.size(); i++)
    {
        const CriticalProperties& props = cachedProps_[i];
        double x_i = Y[i] / props.MW;
        Pc_mix += x_i * props.Pc;
        sum_x += x_i;
    }

    return Pc_mix / (sum_x + 1e-20);
}

double VdWOneFluidMixing::omega(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    double omega_mix = 0.0;
    double sum_x = 0.0;

    for (size_t i = 0; i < cachedProps_.size() && i < Y.size(); i++)
    {
        const CriticalProperties& props = cachedProps_[i];
        double x_i = Y[i] / props.MW;
        omega_mix += x_i * props.omega;
        sum_x += x_i;
    }

    return omega_mix / (sum_x + 1e-20);
}

double VdWOneFluidMixing::MW(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    double MW_inv = 0.0;

    for (size_t i = 0; i < cachedProps_.size() && i < Y.size(); i++)
    {
        const CriticalProperties& props = cachedProps_[i];
        MW_inv += Y[i] / props.MW;
    }

    if (MW_inv < 1e-20) return invalidValue();

    return 1.0 / MW_inv;
}

double VdWOneFluidMixing::a_mix(const std::vector<double>& Y, double T) const
{
    if (!isValidComposition(Y) || !std::isfinite(T) || T <= 0.0)
        return invalidValue();
    // PLACEHOLDER (flagged, see [[rfam-species-energy-foundation]]): pure-species
    // a from index 1 (C3H8 in the default {O2,C3H8} set -> legacy value byte-for-
    // byte). For a >2-species global set index 1 is arbitrary (e.g. N2), so
    // absolute temperature from the EOS reback is NOT quantitatively trustworthy;
    // structure/conservation ARE. Upgrade path (Stage 3): true multi-species
    //   a_mix = SUM_i SUM_j x_i x_j sqrt(a_i a_j)(1-k_ij).
    (void)Y;
    const CriticalProperties& props = cachedProps_.size() > 1 ? cachedProps_[1] : cachedProps_[0];
    double Tr = T / props.Tc;
    double kappa = 0.37464 + 1.54226*props.omega - 0.26992*props.omega*props.omega;
    double alpha = std::pow(1.0 + kappa*(1.0 - std::sqrt(Tr)), 2);
    double a = 0.45724 * R_UNIVERSAL*R_UNIVERSAL*props.Tc*props.Tc / props.Pc * alpha;

    return a;
}

double VdWOneFluidMixing::b_mix(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    // PLACEHOLDER (see a_mix): pure-species b from index 1. Default {O2,C3H8}
    // -> C3H8, legacy value byte-for-byte. Upgrade path: b_mix = SUM_i x_i b_i.
    (void)Y;
    const CriticalProperties& props = cachedProps_.size() > 1 ? cachedProps_[1] : cachedProps_[0];
    double b = 0.07780 * R_UNIVERSAL*props.Tc / props.Pc;

    return b;
}

// ---------------------------------------------------------------------------
// Stage 3 step 3: true quadratic van der Waals one-fluid mixing rules.
//
// Per-species helpers use the exact derivatives of the same alpha function
// used to evaluate a_i(T). The shared implementation is also used by the
// pseudo-critical EOS path so the caloric and pressure surfaces cannot drift.
// ---------------------------------------------------------------------------

namespace {
// a_i(T) = 0.45724*R^2*Tc_i^2/Pc_i * (1 + kappa_i*(1-sqrt(T/Tc_i)))^2
static double speciesA(const CriticalProperties& p, double T)
{
    return evaluatePRAlpha(p.Tc, p.Pc, p.omega, T).a;
}

// b_i = 0.07780*R*Tc_i/Pc_i  (temperature-independent)
static double speciesB(const CriticalProperties& p)
{
    return 0.07780 * R_UNIVERSAL * p.Tc / p.Pc;
}

static double speciesDaDT(const CriticalProperties& p, double T)
{
    return evaluatePRAlpha(p.Tc, p.Pc, p.omega, T).dadT;
}

static double speciesD2aDT2(const CriticalProperties& p, double T)
{
    return evaluatePRAlpha(p.Tc, p.Pc, p.omega, T).d2adT2;
}

// Convert mass fractions Y to mole fractions x. MW must be pre-computed.
static std::vector<double> massToMole(
    const std::vector<CriticalProperties>& props,
    const std::vector<double>& Y)
{
    size_t n = props.size();
    std::vector<double> x(n);
    double sumInv = 0.0;
    for (size_t i = 0; i < n && i < Y.size(); ++i)
        sumInv += Y[i] / props[i].MW;
    if (sumInv < 1e-30) { x.assign(n, 1.0/n); return x; }
    for (size_t i = 0; i < n && i < Y.size(); ++i)
        x[i] = (Y[i] / props[i].MW) / sumInv;
    return x;
}
} // anonymous namespace

double VdWOneFluidMixing::a_quad(
    const std::vector<double>& Y, double T) const
{
    if (!isValidComposition(Y) || !std::isfinite(T) || T <= 0.0)
        return invalidValue();
    // a_mix = sum_i sum_j x_i*x_j*a_ij,  a_ij = sqrt(a_i*a_j)*(1-k_ij)
    std::vector<double> x = massToMole(cachedProps_, Y);
    size_t n = cachedProps_.size();

    double amix = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        double ai = speciesA(cachedProps_[i], T);
        for (size_t j = 0; j < n; ++j)
        {
            double aj  = speciesA(cachedProps_[j], T);
            double kij = kij_.get(speciesNames_[i], speciesNames_[j]);
            amix += x[i] * x[j] * std::sqrt(ai * aj) * (1.0 - kij);
        }
    }
    return amix;
}

double VdWOneFluidMixing::b_lin(const std::vector<double>& Y) const
{
    if (!isValidComposition(Y)) return invalidValue();
    // b_mix = sum_i x_i*b_i  (linear rule, standard for PR)
    std::vector<double> x = massToMole(cachedProps_, Y);
    size_t n = cachedProps_.size();

    double bmix = 0.0;
    for (size_t i = 0; i < n; ++i)
        bmix += x[i] * speciesB(cachedProps_[i]);
    return bmix;
}

double VdWOneFluidMixing::dadT_quad(
    const std::vector<double>& Y, double T) const
{
    if (!isValidComposition(Y) || !std::isfinite(T) || T <= 0.0)
        return invalidValue();
    // da_mix/dT = sum_i sum_j x_i*x_j * da_ij/dT
    // da_ij/dT = (1-kij)/(2*sqrt(ai*aj)) * (da_i/dT * aj + ai * da_j/dT)
    std::vector<double> x = massToMole(cachedProps_, Y);
    size_t n = cachedProps_.size();

    double dadT = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        double ai    = speciesA(cachedProps_[i], T);
        double daidT = speciesDaDT(cachedProps_[i], T);
        for (size_t j = 0; j < n; ++j)
        {
            double aj    = speciesA(cachedProps_[j], T);
            double dajdT = speciesDaDT(cachedProps_[j], T);
            double sqaij = std::sqrt(ai * aj);
            if (sqaij < 1e-30) continue;
            double kij   = kij_.get(speciesNames_[i], speciesNames_[j]);
            dadT += x[i] * x[j] * (1.0 - kij)
                  * (daidT * aj + ai * dajdT) / (2.0 * sqaij);
        }
    }
    return dadT;
}

double VdWOneFluidMixing::d2adT2_quad(
    const std::vector<double>& Y, double T) const
{
    if (!isValidComposition(Y) || !std::isfinite(T) || T <= 0.0)
        return invalidValue();
    // d2a_mix/dT2 = sum_i sum_j x_i*x_j * d2a_ij/dT2
    // Using product/chain rule on da_ij/dT = (1-kij)*(ai*daj+aj*dai)/(2*sqaij):
    //   d2a_ij/dT2 = (1-kij) * [
    //       (d2ai*aj + 2*dai*daj + ai*d2aj) / (2*sqaij)
    //     - (ai*daj+aj*dai)^2 / (4*(ai*aj)^(3/2))
    //   ]
    std::vector<double> x = massToMole(cachedProps_, Y);
    size_t n = cachedProps_.size();

    double d2adT2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        double ai     = speciesA(cachedProps_[i], T);
        double daidT  = speciesDaDT(cachedProps_[i], T);
        double d2aidT = speciesD2aDT2(cachedProps_[i], T);
        for (size_t j = 0; j < n; ++j)
        {
            double aj     = speciesA(cachedProps_[j], T);
            double dajdT  = speciesDaDT(cachedProps_[j], T);
            double d2ajdT = speciesD2aDT2(cachedProps_[j], T);
            double aij    = ai * aj;
            double sqaij  = std::sqrt(aij);
            if (sqaij < 1e-30) continue;
            double kij    = kij_.get(speciesNames_[i], speciesNames_[j]);
            double cross  = ai*dajdT + aj*daidT;
            d2adT2 += x[i] * x[j] * (1.0 - kij) * (
                (d2aidT*aj + 2.0*daidT*dajdT + ai*d2ajdT) / (2.0*sqaij)
                - cross*cross / (4.0*sqaij*aij)
            );
        }
    }
    return d2adT2;
}

} // End namespace Foam
