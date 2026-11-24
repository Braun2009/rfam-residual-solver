#ifndef NASA7Thermo_H
#define NASA7Thermo_H

/*---------------------------------------------------------------------------*\
    NASA7Thermo - ideal-gas heat capacity / enthalpy from NASA-7 polynomials.

    PURPOSE (Stage 3, step 1 - thermodynamic foundation):
    Replaces the hard-coded placeholder `Cp_ideal = 1000 J/(kg K)` in
    RealFluidState with real, TEMPERATURE-dependent, COMPOSITION-dependent
    ideal-gas thermodynamics. This is what presses the reacting-box ignition
    temperature back from the 15218 K placeholder to a physical value.

    DATA SOURCE: GRI-Mech 3.0 thermo (all 9 RFAM species are in it -> ONE table,
    no second set of coefficients drifting; cf. the two-MW-table pitfall).
    Each species: two temperature ranges [Tlow, Tcommon] and [Tcommon, Thigh]
    (Tcommon = 1000 K), 7 coefficients each:
        Cp/R    = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        H/(R T) = a1 + a2/2 T + a3/3 T^2 + a4/4 T^3 + a5/5 T^4 + a6/T   (absolute, incl. hf)
    Outside [Tlow, Thigh] the nearest range's polynomial is extrapolated
    (documented; NASA-7 stays well-behaved a little past its range).

    ENERGY CONVENTION (must match RealFluidState / SpeciesRegistry):
    RealFluidState::e is SENSIBLE energy and does NOT include formation
    enthalpy. hsens(T) = H(T) - H(Tref) with Tref = 298.15 K subtracts the
    formation enthalpy automatically (H(298.15) == hf per species), so the
    formation enthalpy stays OUT of the flow energy and is released explicitly
    by chemistry as qdot. [[rfam-species-energy-foundation]]

    MW authority: molecular weights come from SpeciesRegistry (the single
    authority), NOT a private table here.

    Indexed by NAME (aligned with SpeciesRegistry) - no implicit species order.
    Pure std C++ (header-only). double == OpenFOAM scalar, std::string == word.
\*---------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include "SandiegoNasa7Data.H"
#include "SpeciesRegistry.H"

namespace Foam
{

class NASA7Thermo
{
public:
    struct Poly
    {
        double Tlow, Tcommon, Thigh;
        double lo[7];   // valid on [Tlow, Tcommon]
        double hi[7];   // valid on [Tcommon, Thigh]
    };

private:
    // MUST match RealFluidState::R_UNIVERSAL so Cp/enthalpy are self-consistent.
    static constexpr double Runiv = 8.314462618;   // [J/(mol K)]
    static constexpr double Tref  = 298.15;         // sensible-enthalpy reference [K]

    std::vector<std::string> names_;
    std::vector<Poly>        poly_;   // 1:1 with names_
    SpeciesRegistry          reg_;    // single MW authority (no private MW table)

    int idx(const std::string& name) const
    {
        for (int i = 0; i < (int)names_.size(); ++i)
            if (names_[i] == name) return i;
        throw std::runtime_error("NASA7Thermo: unknown species '" + name + "'");
    }

    // Clamp T into the polynomial's valid range. ABOVE Thigh (below Tlow) the
    // quartic Cp is NOT extrapolated (it can turn negative and break the
    // monotonicity of e(T) that the energy inversion relies on); instead Cp is
    // frozen at the boundary value and h is extended LINEARLY. This is the
    // standard robust treatment and keeps Cp > 0 at any temperature.
    static double clampT(const Poly& p, double T)
    {
        return (T < p.Tlow) ? p.Tlow : (T > p.Thigh ? p.Thigh : T);
    }
    static const double* rangeCoeffs(const Poly& p, double Tc)   // Tc must be in-range
    {
        return (Tc < p.Tcommon) ? p.lo : p.hi;
    }
    static double cpR_at(const double* a, double T)              // Cp/R at in-range T
    {
        return a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])));
    }
    static double hOverR_at(const double* a, double T)           // H/R [K] at in-range T
    {
        return T*(a[0] + T*(a[1]/2.0 + T*(a[2]/3.0 + T*(a[3]/4.0 + T*(a[4]/5.0))))) + a[5];
    }

public:
    static constexpr double sensibleReferenceTemperature = 298.15;

    NASA7Thermo() { build(); }

    bool has(const std::string& name) const
    {
        for (const auto& n : names_) if (n == name) return true;
        return false;
    }

    bool temperatureRange(
        const std::string& name,
        double& minimum,
        double& maximum) const
    {
        for (std::size_t i = 0; i < names_.size(); ++i)
        {
            if (names_[i] != name) continue;
            minimum = poly_[i].Tlow;
            maximum = poly_[i].Thigh;
            return true;
        }
        return false;
    }

    const std::vector<std::string>& speciesNames() const { return names_; }
    const std::vector<Poly>& polynomials() const { return poly_; }
    const SpeciesRegistry& speciesRegistry() const { return reg_; }

    static std::string provenanceManifest(
        const std::vector<std::string>& names,
        const std::vector<Poly>& polynomials,
        const std::vector<SpeciesRegistry::Entry>& registryEntries)
    {
        std::string manifest = "schema=rfam.nasa7-gri30/v1;";
        thermoProvenance::appendDouble(manifest, "Runiv", Runiv);
        thermoProvenance::appendDouble(manifest, "Tref", Tref);
        thermoProvenance::appendSize(manifest, "species.count", names.size());
        thermoProvenance::appendSize(
            manifest, "polynomial.count", polynomials.size());
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            const std::string prefix = "species[" + std::to_string(i) + "].";
            thermoProvenance::appendString(
                manifest, prefix + "name", names[i]);
            if (i >= polynomials.size()) continue;
            const Poly& polynomial = polynomials[i];
            thermoProvenance::appendDouble(
                manifest, prefix + "Tlow", polynomial.Tlow);
            thermoProvenance::appendDouble(
                manifest, prefix + "Tcommon", polynomial.Tcommon);
            thermoProvenance::appendDouble(
                manifest, prefix + "Thigh", polynomial.Thigh);
            for (std::size_t coefficient = 0; coefficient < 7; ++coefficient)
            {
                thermoProvenance::appendDouble(
                    manifest,
                    prefix + "lo[" + std::to_string(coefficient) + "]",
                    polynomial.lo[coefficient]);
                thermoProvenance::appendDouble(
                    manifest,
                    prefix + "hi[" + std::to_string(coefficient) + "]",
                    polynomial.hi[coefficient]);
            }
        }
        thermoProvenance::appendString(
            manifest,
            "speciesRegistry.sha256",
            SpeciesRegistry::provenanceSha256(registryEntries));
        return manifest;
    }

    static std::string provenanceSha256(
        const std::vector<std::string>& names,
        const std::vector<Poly>& polynomials,
        const std::vector<SpeciesRegistry::Entry>& registryEntries)
    {
        return thermoProvenance::sha256Hex(
            provenanceManifest(names, polynomials, registryEntries));
    }

    std::string provenanceManifest() const
    {
        return provenanceManifest(names_, poly_, reg_.entries());
    }

    std::string provenanceSha256() const
    {
        return provenanceSha256(names_, poly_, reg_.entries());
    }

    // molecular weight [kg/mol] from the registry (single authority)
    double MW(const std::string& name) const
    {
        return reg_.molecularWeight(reg_.index(name));
    }

    // ---- per-species (molar basis) ----

    // Cp [J/(mol K)] — frozen at the boundary value outside [Tlow,Thigh]
    double cp_mole(const std::string& name, double T) const
    {
        const Poly& p = poly_[idx(name)];
        const double Tc = clampT(p, T);
        return Runiv * cpR_at(rangeCoeffs(p, Tc), Tc);
    }

    // absolute molar enthalpy (INCLUDES formation enthalpy) [J/mol].
    // Outside the valid range: h is continued LINEARLY with the boundary Cp so
    // that dh/dT == cp holds everywhere and e(T) stays strictly monotone.
    double h_mole(const std::string& name, double T) const
    {
        const Poly& p = poly_[idx(name)];
        const double Tc = clampT(p, T);
        const double* a = rangeCoeffs(p, Tc);
        double h = Runiv * hOverR_at(a, Tc);
        if (T != Tc) h += Runiv * cpR_at(a, Tc) * (T - Tc);   // linear tail
        return h;
    }

    // sensible molar enthalpy referenced to Tref=298.15 K [J/mol]
    // (subtracts H(Tref) == hf -> formation enthalpy is removed)
    double hsens_mole(const std::string& name, double T) const
    {
        return h_mole(name, T) - h_mole(name, Tref);
    }

    // standard molar entropy [J/(mol K)] — NASA-7 polynomial: S/R = a1 ln(T) + a2 T + ...
    // Note: this is the standard entropy S°(T); for equilibrium constant calculation
    // we need ΔS° = S°_products - S°_reactants at temperature T.
    double s_mole(const std::string& name, double T) const
    {
        const Poly& p = poly_[idx(name)];
        const double Tc = clampT(p, T);
        const double* a = rangeCoeffs(p, Tc);
        double sR = a[0]*std::log(Tc) + Tc*(a[1] + Tc*(a[2]/2.0 + Tc*(a[3]/3.0 + Tc*(a[4]/4.0)))) + a[6];
        double s = Runiv * sR;
        // Linear extrapolation outside range: ds/dT = Cp(Tclamp)/T
        if (T != Tc)
        {
            double cp_at_bound = Runiv * cpR_at(a, Tc);
            s += cp_at_bound / Tc * (T - Tc);
        }
        return s;
    }

    // ---- mixture (mass basis), for the EOS ideal-gas part ----

    // mass-specific mixture Cp [J/(kg K)] = sum_i Y_i * cp_mole_i / MW_i
    double cp_mass_mix
    (
        const std::vector<std::string>& names,
        const std::vector<double>&      Y,
        double                          T
    ) const
    {
        double cp = 0.0;
        for (size_t i = 0; i < names.size(); ++i)
        {
            if (Y[i] == 0.0) continue;
            cp += Y[i] * cp_mole(names[i], T) / MW(names[i]);
        }
        return cp;
    }

    // mass-specific SENSIBLE mixture enthalpy [J/kg] (ref 298.15 K)
    double hsens_mass_mix
    (
        const std::vector<std::string>& names,
        const std::vector<double>&      Y,
        double                          T
    ) const
    {
        double h = 0.0;
        for (size_t i = 0; i < names.size(); ++i)
        {
            if (Y[i] == 0.0) continue;
            h += Y[i] * hsens_mole(names[i], T) / MW(names[i]);
        }
        return h;
    }

private:
    void add
    (
        const std::string& n,
        double Tlow, double Tcom, double Thigh,
        std::initializer_list<double> lo,
        std::initializer_list<double> hi
    )
    {
        Poly p; p.Tlow = Tlow; p.Tcommon = Tcom; p.Thigh = Thigh;
        int k = 0; for (double v : lo) p.lo[k++] = v;
        k = 0;     for (double v : hi) p.hi[k++] = v;
        names_.push_back(n);
        poly_.push_back(p);
    }

    // GRI-Mech 3.0 NASA-7 coefficients. lo=[Tlow,1000], hi=[1000,Thigh].
    // Verified by hand vs textbook Cp(298 K): O2 918, CO2 844, H2O 1864 J/(kg K).
    void build()
    {
        add("O2", 200.0, 1000.0, 3500.0,
            { 3.78245636e+00, -2.99673416e-03,  9.84730201e-06, -9.68129509e-09,  3.24372837e-12, -1.06394356e+03,  3.65767573e+00},
            { 3.28253784e+00,  1.48308754e-03, -7.57966669e-07,  2.09470555e-10, -2.16717794e-14, -1.08845772e+03,  5.45323129e+00});

        add("N2", 300.0, 1000.0, 5000.0,
            { 3.29867700e+00,  1.40824040e-03, -3.96322200e-06,  5.64151500e-09, -2.44485400e-12, -1.02089990e+03,  3.95037200e+00},
            { 2.92664000e+00,  1.48797680e-03, -5.68476000e-07,  1.00970380e-10, -6.75335100e-15, -9.22797700e+02,  5.98052800e+00});

        add("H2", 200.0, 1000.0, 3500.0,
            { 2.34433112e+00,  7.98052075e-03, -1.94781510e-05,  2.01572094e-08, -7.37611761e-12, -9.17935173e+02,  6.83010238e-01},
            { 3.33727920e+00, -4.94024731e-05,  4.99456778e-07, -1.79566394e-10,  2.00255376e-14, -9.50158922e+02, -3.20502331e+00});

        add("CH4", 200.0, 1000.0, 3500.0,
            { 5.14987613e+00, -1.36709788e-02,  4.91800599e-05, -4.84743026e-08,  1.66693956e-11, -1.02466476e+04, -4.64130376e+00},
            { 7.48514950e-02,  1.33909467e-02, -5.73285809e-06,  1.22292535e-09, -1.01815230e-13, -9.46834459e+03,  1.84373180e+01});

        add("C2H6", 200.0, 1000.0, 3500.0,
            { 4.29142492e+00, -5.50154270e-03,  5.99438288e-05, -7.08466285e-08,  2.68685771e-11, -1.15222055e+04,  2.66682316e+00},
            { 1.07188150e+00,  2.16852677e-02, -1.00256067e-05,  2.21412001e-09, -1.90002890e-13, -1.14263932e+04,  1.51156107e+01});

        add("C3H8", 300.0, 1000.0, 5000.0,
            { 9.33553810e-01,  2.64245790e-02,  6.10597270e-06, -2.19774990e-08,  9.51492530e-12, -1.39585200e+04,  1.92016910e+01},
            { 7.53413680e+00,  1.88722390e-02, -6.27184910e-06,  9.14756490e-10, -4.78380690e-14, -1.64675160e+04, -1.78923490e+01});

        add("H2O", 200.0, 1000.0, 3500.0,
            { 4.19864056e+00, -2.03643410e-03,  6.52040211e-06, -5.48797062e-09,  1.77197817e-12, -3.02937267e+04, -8.49032208e-01},
            { 3.03399249e+00,  2.17691804e-03, -1.64072518e-07, -9.70419870e-11,  1.68200992e-14, -3.00042971e+04,  4.96677010e+00});

        add("CO", 200.0, 1000.0, 3500.0,
            { 3.57953347e+00, -6.10353680e-04,  1.01681433e-06,  9.07005884e-10, -9.04424499e-13, -1.43440860e+04,  3.50840928e+00},
            { 2.71518561e+00,  2.06252743e-03, -9.98825771e-07,  2.30053008e-10, -2.03647716e-14, -1.41518724e+04,  7.81868772e+00});

        add("CO2", 200.0, 1000.0, 3500.0,
            { 2.35677352e+00,  8.98459677e-03, -7.12356269e-06,  2.45919022e-09, -1.43699548e-13, -4.83719697e+04,  9.90105222e+00},
            { 3.85746029e+00,  4.41437026e-03, -2.21481404e-06,  5.23490188e-10, -4.72084164e-14, -4.87591660e+04,  2.27163806e+00});

        // ---- mechanism species beyond the GRI nine -------------------------
        // The multi-species EOS evaluates every species of the loaded mechanism
        // through this table, and the San Diego mechanism carries 50 more than
        // the GRI set. Without them the EOS threw 'unknown species', which the
        // legacy model's catch-all converted into INVALID_CONFIGURATION, so a
        // 59-species run failed at supply initialization with a message that
        // named neither the species nor the table. The polynomials are the same
        // ones the chemistry operator integrates (SandiegoNasa7Data.H), so the
        // energy the EOS inverts and the heat release the kinetics reports come
        // from one source.
        for(int i=0;i<SandiegoNasa7::nSpecies;++i)
        {
            const auto& s=SandiegoNasa7::species[i];
            if(has(s.name)) continue;
            add(s.name,s.Tlow,s.Tcommon,s.Thigh,
                { s.lo[0],s.lo[1],s.lo[2],s.lo[3],s.lo[4],s.lo[5],s.lo[6] },
                { s.hi[0],s.hi[1],s.hi[2],s.hi[3],s.hi[4],s.hi[5],s.hi[6] });
        }
    }
};

// Process-wide singleton (coefficients are read-only, so this is safe & cheap).
inline const NASA7Thermo& nasa7Thermo()
{
    static const NASA7Thermo instance;
    return instance;
}

} // End namespace Foam

#endif
