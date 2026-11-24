/*---------------------------------------------------------------------------*  realFluidAllMach - van der Waals One-Fluid Mixing Rules
\*---------------------------------------------------------------------------*/

#ifndef VdWOneFluidMixing_H
#define VdWOneFluidMixing_H

#include "SpeciesCriticalDB.H"
#include "BinaryInteractionDB.H"
#include <vector>
#include <string>

namespace Foam
{

class VdWOneFluidMixing
{
private:
    const SpeciesCriticalDB& criticalData_;
    const BinaryInteractionDB& kij_;

    std::vector<std::string> speciesNames_;
    // Validated critical props, 1:1 with speciesNames_ (populated in
    // setSpeciesNames). Compute methods read this cache, NOT per-call DB
    // lookups, so an unknown species is caught ONCE, up front, not silently.
    std::vector<CriticalProperties> cachedProps_;

    // Stage 3 step 3: true quadratic VdW mixing. Set false in ctor (legacy
    // path: pseudo-critical a/b from Tc_mix/Pc_mix -> MVP iter=293 unchanged).
    // Set true ONLY via enableMultiSpeciesMixing(), which is called from
    // PengRobinsonEOS::setSpeciesNames (the external opt-in path). The
    // internal ctor call to setSpeciesNames does NOT flip this flag so that
    // the default {O2,C3H8} EOS path stays byte-for-byte identical.
    bool multiSpeciesEnabled_;

public:
    VdWOneFluidMixing(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij
    );

    void setSpeciesNames(const std::vector<std::string>& names);

    // Stage 3 step 3: flip to true quadratic mixing. Called by
    // PengRobinsonEOS::setSpeciesNames (external opt-in). NOT called from the
    // internal ctor setSpeciesNames so the default {O2,C3H8} path stays legacy.
    void enableMultiSpeciesMixing() { multiSpeciesEnabled_ = true; }
    bool multiSpeciesEnabled() const { return multiSpeciesEnabled_; }
    bool isValidComposition(const std::vector<double>& Y) const;

    double Tc(const std::vector<double>& Y) const;
    double Pc(const std::vector<double>& Y) const;
    double omega(const std::vector<double>& Y) const;
    double MW(const std::vector<double>& Y) const;

    // Legacy scalar a/b (placeholder — still used by old a_mix()/b_mix() callers).
    double a_mix(const std::vector<double>& Y, double T) const;
    double b_mix(const std::vector<double>& Y) const;

    // Stage 3 step 3: true quadratic VdW mixing (only when multiSpeciesEnabled_).
    // Per-species alpha: alpha_i = (1 + kappa_i*(1-sqrt(T/Tc_i)))^2
    // a_ij = sqrt(a_i*a_j)*(1 - kij)
    // a_mix = sum_i sum_j x_i*x_j*a_ij         (mole-fraction quadratic)
    // b_mix = sum_i x_i*b_i                    (linear)
    // da_mix/dT  = sum_i sum_j x_i*x_j * d(a_ij)/dT   (analytic)
    // d2a_mix/dT2 = analytic second derivative
    double a_quad(const std::vector<double>& Y, double T) const;
    double b_lin(const std::vector<double>& Y) const;
    double dadT_quad(const std::vector<double>& Y, double T) const;
    double d2adT2_quad(const std::vector<double>& Y, double T) const;
};

} // End namespace Foam

#endif
