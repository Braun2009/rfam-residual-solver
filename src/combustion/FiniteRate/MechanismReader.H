#ifndef MechanismReader_H
#define MechanismReader_H

/*---------------------------------------------------------------------------*\
    MechanismReader - minimal YAML-subset chemistry mechanism file parser.

    Reads a simple YAML-like mechanism file (no third-party libraries).
    Sufficient to reproduce TwoStepCO (C3H8/O2, Jones & Lindstedt 1988)
    and any global multi-step mechanism with Arrhenius rates.

    File format example (see tests/twostep_co.yaml):
    ---
    mechanism: TwoStepCO
    species:
      - name: C3H8
        MW: 0.0440956    # [kg/mol]
        hf: -103850.0    # standard formation enthalpy [J/mol]
      - name: O2
        MW: 0.0319988
        hf: 0.0
      ...

    reactions:
      - name: R1_fuel
        A: 1.0e10
        n: 0.0
        Ea: 125520.0     # [J/mol]
        DH: -1.19502e6   # heat of reaction [J/mol_of_limiting_reactant]
        conc_exp:        # concentration exponents (Arrhenius rate)
          C3H8: 0.1
          O2:   1.65
        stoich:          # net stoichiometry (positive=product, negative=reactant)
          C3H8: -1.0
          O2:   -3.5
          CO:   +3.0
          H2O:  +4.0
      ...
    ---

    Pure std C++, no dependencies. Returns a MechanismData struct.
\*---------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace Foam
{

// ---- Data structures returned by the reader ----

struct SpeciesData
{
    std::string name;
    double MW = 0.0;    // [kg/mol]
    double hf = 0.0;    // standard formation enthalpy [J/mol]
};

struct ReactionData
{
    std::string name;
    double A = 0.0;      // pre-exponential [SI]
    double n = 0.0;      // temperature exponent
    double Ea = 0.0;     // activation energy [J/mol]
    double DH = 0.0;     // heat of reaction [J/mol_limiting]
    bool reversible = false;

    std::map<std::string, double> concExp;    // concentration exponents
    std::map<std::string, double> stoich;     // net stoichiometry
};

struct MechanismData
{
    std::string name;
    std::vector<SpeciesData>  species;
    std::vector<ReactionData> reactions;

    int speciesIndex(const std::string& n) const
    {
        for (int i = 0; i < (int)species.size(); ++i)
            if (species[i].name == n) return i;
        return -1;
    }
};

// ---- Parser ----

class MechanismReader
{
    // ---- string helpers ----
    static std::string trim(const std::string& s)
    {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }
    static bool startsWith(const std::string& s, const std::string& pre)
    {
        return s.size() >= pre.size() && s.substr(0, pre.size()) == pre;
    }
    static double parseDouble(const std::string& s)
    {
        const std::string value = trim(s);
        std::size_t parsed = 0;
        double result = 0.0;
        try
        {
            result = std::stod(value, &parsed);
        }
        catch (const std::exception&)
        {
            throw std::runtime_error(
                "MechanismReader: invalid numeric value '" + value + "'");
        }
        if (parsed != value.size())
            throw std::runtime_error(
                "MechanismReader: invalid numeric value '" + value + "'");
        return result;
    }
    static bool parseBool(const std::string& s)
    {
        const std::string value = trim(s);
        if (value == "true") return true;
        if (value == "false") return false;
        throw std::runtime_error(
            "MechanismReader: boolean value must be 'true' or 'false'");
    }
    // "key: value" -> {key, value}  or  "  - name: X" -> {"name","X"}
    static bool parseKeyVal(const std::string& line,
                            std::string& key, std::string& val)
    {
        size_t pos = line.find(':');
        if (pos == std::string::npos) return false;
        key = trim(line.substr(0, pos));
        val = trim(line.substr(pos + 1));
        if (key.size() > 0 && key[0] == '-') key = trim(key.substr(1));
        return true;
    }

    enum Section { NONE, SPECIES, REACTIONS, REACTION_BODY,
                   CONC_EXP, STOICH };

public:
    static MechanismData read(const std::string& filename)
    {
        std::ifstream f(filename);
        if (!f.good())
            throw std::runtime_error("MechanismReader: cannot open '" + filename + "'");

        MechanismData mech;
        Section section = NONE;
        SpeciesData  curSpec;
        ReactionData curRxn;
        struct SpeciesPresence
        {
            bool MW = false;
            bool hf = false;
        } specPresence;
        struct ReactionPresence
        {
            bool A = false;
            bool n = false;
            bool Ea = false;
            bool DH = false;
            bool reversible = false;
            bool concExp = false;
            bool stoich = false;
        } rxnPresence;
        bool inSpec = false;
        bool inRxn = false;
        bool hasMechanism = false;
        bool hasSpeciesSection = false;
        bool hasReactionSection = false;

        auto flushSpecies = [&]()
        {
            if (!inSpec) return;
            if (curSpec.name.empty())
                throw std::runtime_error(
                    "MechanismReader: species name must not be empty");
            if (!specPresence.MW)
                throw std::runtime_error(
                    "MechanismReader: species '" + curSpec.name
                    + "' is missing required field 'MW'");
            if (!specPresence.hf)
                throw std::runtime_error(
                    "MechanismReader: species '" + curSpec.name
                    + "' is missing required field 'hf'");
            mech.species.push_back(curSpec);
            curSpec = SpeciesData{};
            specPresence = SpeciesPresence{};
            inSpec = false;
        };

        auto flushReaction = [&]()
        {
            if (!inRxn) return;
            if (curRxn.name.empty())
                throw std::runtime_error(
                    "MechanismReader: reaction name must not be empty");
            if (!rxnPresence.A)
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required field 'A'");
            if (!rxnPresence.n)
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required field 'n'");
            if (!rxnPresence.Ea)
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required field 'Ea'");
            if (!rxnPresence.DH)
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required field 'DH'");
            if (!rxnPresence.concExp || curRxn.concExp.empty())
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required section 'conc_exp'");
            if (!rxnPresence.stoich || curRxn.stoich.empty())
                throw std::runtime_error(
                    "MechanismReader: reaction '" + curRxn.name
                    + "' is missing required section 'stoich'");
            mech.reactions.push_back(curRxn);
            curRxn = ReactionData{};
            rxnPresence = ReactionPresence{};
            inRxn = false;
        };

        std::string line;
        while (std::getline(f, line))
        {
            // Remove comments
            size_t cpos = line.find('#');
            if (cpos != std::string::npos) line = line.substr(0, cpos);
            if (trim(line).empty() || trim(line) == "---") continue;

            std::string tl = trim(line);

            // Section headers
            if (tl == "species:")
            {
                if (section != NONE || !hasMechanism || hasSpeciesSection)
                    throw std::runtime_error(
                        "MechanismReader: invalid or duplicate species section");
                hasSpeciesSection = true;
                section = SPECIES;
                continue;
            }
            if (tl == "reactions:")
            {
                if (section != SPECIES || !hasSpeciesSection
                    || hasReactionSection)
                    throw std::runtime_error(
                        "MechanismReader: invalid or duplicate reactions section");
                flushSpecies();
                hasReactionSection = true;
                section = REACTIONS;
                continue;
            }
            if (inRxn && section == REACTION_BODY
                && (tl == "conc_exp:" || tl == "concExp:"))
            {
                if (rxnPresence.concExp)
                    throw std::runtime_error(
                        "MechanismReader: duplicate conc_exp section in reaction '"
                        + curRxn.name + "'");
                rxnPresence.concExp = true;
                section = CONC_EXP;
                continue;
            }
            if (inRxn && (section == REACTION_BODY || section == CONC_EXP)
                && tl == "stoich:")
            {
                if (rxnPresence.stoich)
                    throw std::runtime_error(
                        "MechanismReader: duplicate stoich section in reaction '"
                        + curRxn.name + "'");
                rxnPresence.stoich = true;
                section = STOICH;
                continue;
            }

            // Top-level key: value
            if (section == NONE)
            {
                std::string k, v;
                if (!parseKeyVal(tl, k, v))
                    throw std::runtime_error(
                        "MechanismReader: invalid top-level entry '" + tl + "'");
                if (k != "mechanism")
                    throw std::runtime_error(
                        "MechanismReader: unknown top-level key '" + k + "'");
                if (hasMechanism)
                    throw std::runtime_error(
                        "MechanismReader: duplicate mechanism field");
                if (v.empty())
                    throw std::runtime_error(
                        "MechanismReader: mechanism name must not be empty");
                mech.name = v;
                hasMechanism = true;
                continue;
            }

            // --- species section ---
            if (section == SPECIES)
            {
                if (startsWith(tl, "- name:") || startsWith(tl, "-name:"))
                {
                    flushSpecies();
                    curSpec = SpeciesData{};
                    specPresence = SpeciesPresence{};
                    curSpec.name = trim(tl.substr(tl.find(':') + 1));
                    inSpec = true;
                }
                else
                {
                    if (!inSpec)
                        throw std::runtime_error(
                            "MechanismReader: species data appears before a species name");
                    std::string k, v;
                    if (!parseKeyVal(tl, k, v))
                        throw std::runtime_error(
                            "MechanismReader: invalid species entry '" + tl + "'");
                    if (k == "MW" || k == "mw")
                    {
                        if (specPresence.MW)
                            throw std::runtime_error(
                                "MechanismReader: duplicate MW for species '"
                                + curSpec.name + "'");
                        curSpec.MW = parseDouble(v);
                        specPresence.MW = true;
                    }
                    else if (k == "hf")
                    {
                        if (specPresence.hf)
                            throw std::runtime_error(
                                "MechanismReader: duplicate hf for species '"
                                + curSpec.name + "'");
                        curSpec.hf = parseDouble(v);
                        specPresence.hf = true;
                    }
                    else
                        throw std::runtime_error(
                            "MechanismReader: unknown species key '" + k + "'");
                }
                continue;
            }

            // --- reactions section (outer list) ---
            if (section == REACTIONS)
            {
                if (startsWith(tl, "- name:") || startsWith(tl, "-name:"))
                {
                    flushReaction();
                    curRxn = ReactionData{};
                    rxnPresence = ReactionPresence{};
                    curRxn.name = trim(tl.substr(tl.find(':') + 1));
                    inRxn = true;
                    section = REACTION_BODY;
                }
                else
                    throw std::runtime_error(
                        "MechanismReader: reaction data appears before a reaction name");
                continue;
            }

            // --- reaction body: A, n, Ea, DH, reversible ---
            if (section == REACTION_BODY)
            {
                if (startsWith(tl, "- name:") || startsWith(tl, "-name:"))
                {
                    flushReaction();
                    curRxn = ReactionData{};
                    rxnPresence = ReactionPresence{};
                    curRxn.name = trim(tl.substr(tl.find(':') + 1));
                    inRxn = true;
                    continue;
                }
                std::string k, v;
                if (!parseKeyVal(tl, k, v))
                    throw std::runtime_error(
                        "MechanismReader: invalid reaction entry '" + tl + "'");
                if (k == "A")
                {
                    if (rxnPresence.A)
                        throw std::runtime_error(
                            "MechanismReader: duplicate A in reaction '"
                            + curRxn.name + "'");
                    curRxn.A = parseDouble(v);
                    rxnPresence.A = true;
                }
                else if (k == "n")
                {
                    if (rxnPresence.n)
                        throw std::runtime_error(
                            "MechanismReader: duplicate n in reaction '"
                            + curRxn.name + "'");
                    curRxn.n = parseDouble(v);
                    rxnPresence.n = true;
                }
                else if (k == "Ea")
                {
                    if (rxnPresence.Ea)
                        throw std::runtime_error(
                            "MechanismReader: duplicate Ea in reaction '"
                            + curRxn.name + "'");
                    curRxn.Ea = parseDouble(v);
                    rxnPresence.Ea = true;
                }
                else if (k == "DH")
                {
                    if (rxnPresence.DH)
                        throw std::runtime_error(
                            "MechanismReader: duplicate DH in reaction '"
                            + curRxn.name + "'");
                    curRxn.DH = parseDouble(v);
                    rxnPresence.DH = true;
                }
                else if (k == "reversible")
                {
                    if (rxnPresence.reversible)
                        throw std::runtime_error(
                            "MechanismReader: duplicate reversible in reaction '"
                            + curRxn.name + "'");
                    curRxn.reversible = parseBool(v);
                    rxnPresence.reversible = true;
                }
                else
                    throw std::runtime_error(
                        "MechanismReader: unknown reaction key '" + k + "'");
                continue;
            }

            // --- concentration exponents ---
            if (section == CONC_EXP)
            {
                if (tl == "stoich:") { section = STOICH; continue; }
                // detect if a new reaction starts
                if (startsWith(tl, "- name:") || startsWith(tl, "-name:"))
                {
                    flushReaction();
                    curRxn = ReactionData{};
                    rxnPresence = ReactionPresence{};
                    curRxn.name = trim(tl.substr(tl.find(':') + 1));
                    inRxn = true; section = REACTION_BODY; continue;
                }
                std::string k, v;
                if (!parseKeyVal(tl, k, v) || k.empty())
                    throw std::runtime_error(
                        "MechanismReader: invalid concentration exponent entry '"
                        + tl + "'");
                if (!curRxn.concExp.emplace(k, parseDouble(v)).second)
                    throw std::runtime_error(
                        "MechanismReader: duplicate concentration exponent for species '"
                        + k + "' in reaction '" + curRxn.name + "'");
                continue;
            }

            // --- net stoichiometry ---
            if (section == STOICH)
            {
                if (startsWith(tl, "- name:") || startsWith(tl, "-name:"))
                {
                    flushReaction();
                    curRxn = ReactionData{};
                    rxnPresence = ReactionPresence{};
                    curRxn.name = trim(tl.substr(tl.find(':') + 1));
                    inRxn = true; section = REACTION_BODY; continue;
                }
                std::string k, v;
                if (!parseKeyVal(tl, k, v) || k.empty())
                    throw std::runtime_error(
                        "MechanismReader: invalid stoichiometry entry '"
                        + tl + "'");
                if (!curRxn.stoich.emplace(k, parseDouble(v)).second)
                    throw std::runtime_error(
                        "MechanismReader: duplicate stoichiometry for species '"
                        + k + "' in reaction '" + curRxn.name + "'");
                continue;
            }

            throw std::runtime_error(
                "MechanismReader: unsupported entry '" + tl + "'");
        }

        // flush last items
        flushSpecies();
        flushReaction();

        if (!hasMechanism)
            throw std::runtime_error(
                "MechanismReader: missing required mechanism field in '"
                + filename + "'");
        if (mech.species.empty())
            throw std::runtime_error("MechanismReader: no species found in '" + filename + "'");
        if (mech.reactions.empty())
            throw std::runtime_error("MechanismReader: no reactions found in '" + filename + "'");

        return mech;
    }
};

} // End namespace Foam

#endif
