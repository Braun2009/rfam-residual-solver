#ifndef SpeciesRegistry_H
#define SpeciesRegistry_H

/*---------------------------------------------------------------------------*\
    SpeciesRegistry - the SINGLE authority for species identity in RFAM.

    PURPOSE (Stage 1.5 - Species/Energy Consistency):
    Every module that touches species (RealFluidState/EOS, ChemistryModel,
    ConservedVars rhoY ordering, SpeciesMapper) resolves name <-> global index,
    molecular weight, formation enthalpy and element composition THROUGH THIS
    REGISTRY. No module may invent its own implicit ordering (the old
    "species[0] is the fuel" convention is forbidden).

    The global species set is hard-coded here and kept ALIGNED with
    SpeciesCriticalDB::createDefault() so that a global index maps to the same
    species the EOS critical-property lookup uses.

    Energy convention (see [[Stage1.5 StepD]]): formationEnthalpy is provided
    for chemistry/heat-release bookkeeping. RealFluidState::e is SENSIBLE energy
    and does NOT include formation enthalpy; chemistry returns qdot explicitly.

    Pure std C++ : int  == OpenFOAM label, std::string == OpenFOAM word.
\*---------------------------------------------------------------------------*/

#include <string>
#include "SandiegoNasa7Data.H"
#include <vector>

#include "ThermoDataProvenance.H"

namespace Foam
{

class SpeciesRegistry
{
public:
    struct Entry
    {
        std::string name;
        double MW;        // molecular weight [kg/mol]
        double hf;        // standard enthalpy of formation [J/mol] @ 298.15 K
        int    nC, nH, nO, nN;   // element composition (atoms per molecule)
    };

private:
    std::vector<Entry> entries_;

    int findOrNeg(const std::string& name) const
    {
        for (int i = 0; i < (int)entries_.size(); ++i)
            if (entries_[i].name == name) return i;
        return -1;
    }

public:
    SpeciesRegistry() { buildDefault(); }

    // ---- authoritative lookups ----
    int size() const { return (int)entries_.size(); }

    // name -> global index ; throws via -1 sentinel check at call sites
    int index(const std::string& name) const { return findOrNeg(name); }
    bool has(const std::string& name) const { return findOrNeg(name) >= 0; }

    const std::string& name(int g) const { return entries_.at(g).name; }
    double molecularWeight(int g) const   { return entries_.at(g).MW; }   // [kg/mol]
    double formationEnthalpy(int g) const { return entries_.at(g).hf; }   // [J/mol]
    int nC(int g) const { return entries_.at(g).nC; }
    int nH(int g) const { return entries_.at(g).nH; }
    int nO(int g) const { return entries_.at(g).nO; }
    int nN(int g) const { return entries_.at(g).nN; }

    const std::vector<Entry>& entries() const { return entries_; }

    static std::string provenanceManifest(const std::vector<Entry>& entries)
    {
        std::string manifest = "schema=rfam.species-registry/v1;";
        thermoProvenance::appendSize(manifest, "entry.count", entries.size());
        for (std::size_t i = 0; i < entries.size(); ++i)
        {
            const Entry& entry = entries[i];
            const std::string prefix = "entry[" + std::to_string(i) + "].";
            thermoProvenance::appendString(
                manifest, prefix + "name", entry.name);
            thermoProvenance::appendDouble(
                manifest, prefix + "MW", entry.MW);
            thermoProvenance::appendDouble(
                manifest, prefix + "hf", entry.hf);
            thermoProvenance::appendSigned(
                manifest, prefix + "nC", entry.nC);
            thermoProvenance::appendSigned(
                manifest, prefix + "nH", entry.nH);
            thermoProvenance::appendSigned(
                manifest, prefix + "nO", entry.nO);
            thermoProvenance::appendSigned(
                manifest, prefix + "nN", entry.nN);
        }
        return manifest;
    }

    static std::string provenanceSha256(const std::vector<Entry>& entries)
    {
        return thermoProvenance::sha256Hex(provenanceManifest(entries));
    }

    std::string provenanceManifest() const
    {
        return provenanceManifest(entries_);
    }

    std::string provenanceSha256() const
    {
        return provenanceSha256(entries_);
    }

    // Global default registry (ALIGNED with SpeciesCriticalDB::createDefault()).
    // hf [J/mol] standard values; element composition exact.
    void buildDefault()
    {
        entries_.clear();
        // name  MW[kg/mol]  hf[J/mol]  C H O N
        // MW values match the chemistry models' constants to full precision so
        // there is ONE molecular-weight authority (TwoStepCO etc. must migrate
        // to query these; until then the values are kept bit-identical to avoid
        // a second MW table drifting -> element-balance residue). [[Stage1.5]]
        entries_.push_back({"O2",   0.0319988,        0.0,  0, 0, 2, 0});
        entries_.push_back({"N2",   0.0280134,        0.0,  0, 0, 0, 2});
        entries_.push_back({"H2",   0.00201588,       0.0,  0, 2, 0, 0});
        entries_.push_back({"CH4",  0.01604246, -74600.0,  1, 4, 0, 0});
        entries_.push_back({"C2H6", 0.03006904, -84000.0,  2, 6, 0, 0});
        entries_.push_back({"C3H8", 0.04409562,-103850.0,  3, 8, 0, 0});
        entries_.push_back({"H2O",  0.01801528,-241820.0,  0, 2, 1, 0});
        entries_.push_back({"CO",   0.0280101, -110530.0,  1, 0, 1, 0});
        entries_.push_back({"CO2",  0.0440095, -393510.0,  1, 0, 2, 0});
        // Mechanism species beyond the historical nine. Appended, never
        // inserted, so the indices above stay fixed for every caller that
        // assumes the alignment with SpeciesCriticalDB. Needed because the
        // multi-species EOS requires a molecular weight for every species in
        // the mechanism: without these the 59-species San Diego mechanism
        // failed configuration outright (INVALID_CONFIGURATION).
        for(int i=0;i<SandiegoNasa7::nSpecies;++i)
        {
            const auto& s=SandiegoNasa7::species[i];
            if(findOrNeg(s.name)>=0) continue;
            entries_.push_back({s.name,s.mw,nasa7FormationEnthalpy(s),
                                s.nC,s.nH,s.nO,s.nN});
        }
    }

    // Standard enthalpy of formation at 298.15 K from the NASA-7 polynomial
    // h(T) = R T (a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T);
    // the polynomial's h already carries the formation term.
    static double nasa7FormationEnthalpy(const SandiegoNasa7::SandiegoSpeciesData& s)
    {
        constexpr double R=8.314462618;
        const double T=298.15;
        const double* a=(T<=s.Tcommon)?s.lo:s.hi;
        return R*T*(a[0]+a[1]*T/2.0+a[2]*T*T/3.0+a[3]*T*T*T/4.0
                    +a[4]*T*T*T*T/5.0+a[5]/T);
    }

    // Convenience: the canonical global species-name list (insertion order).
    std::vector<std::string> names() const
    {
        std::vector<std::string> out;
        out.reserve(entries_.size());
        for (const auto& e : entries_) out.push_back(e.name);
        return out;
    }
};

} // End namespace Foam

#endif
