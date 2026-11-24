#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "PhaseEquilibrium/TabulatedIdealCaloric.H"

using namespace Foam::phaseEquilibrium;

namespace
{

constexpr const char* expectedSchema =
    "rfam.ideal_caloric_coolprop72/v1";
constexpr const char* expectedNotes =
    "ideal_caloric_only;not_liquid_or_two_phase_transport;not_HEOS_residual";

int failures = 0;

void check(const std::string& name, bool condition)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] " << name
              << '\n';
    if (!condition) ++failures;
}

void checkClose(
    const std::string& name,
    double actual,
    double expected,
    double relativeTolerance,
    double absoluteTolerance)
{
    const double error = std::abs(actual - expected);
    const double tolerance = absoluteTolerance
        + relativeTolerance*std::max(std::abs(actual), std::abs(expected));
    const bool condition = std::isfinite(actual) && error <= tolerance;
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] " << name;
    if (!condition)
    {
        std::cout << " got=" << std::setprecision(17) << actual
                  << " expected=" << expected
                  << " error=" << error
                  << " tolerance=" << tolerance;
        ++failures;
    }
    std::cout << '\n';
}

std::vector<std::string> split(const std::string& line)
{
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

struct FrozenRow
{
    std::string species;
    double temperature;
    double sensibleEnthalpy;
    double cp;
};

std::filesystem::path locateReference(int argc, char** argv)
{
    if (argc > 1) return argv[1];

    const std::filesystem::path sourcePath(__FILE__);
    const std::vector<std::filesystem::path> candidates =
    {
        "tests/reference/ideal_caloric_coolprop72.csv",
        "../tests/reference/ideal_caloric_coolprop72.csv",
        "../../tests/reference/ideal_caloric_coolprop72.csv",
        sourcePath.parent_path().parent_path()
            / "reference/ideal_caloric_coolprop72.csv"
    };
    for (const auto& candidate : candidates)
    {
        if (std::filesystem::exists(candidate)) return candidate;
    }
    return candidates.front();
}

std::vector<FrozenRow> loadReference(const std::filesystem::path& path)
{
    std::ifstream file(path);
    if (!file)
    {
        check("frozen reference opens: " + path.string(), false);
        return {};
    }

    std::string line;
    std::getline(file, line);
    const std::string expectedHeader =
        "schema_version,species,coolprop_fluid,backend,coolprop_version,"
        "source_api,enthalpy_property,cp_property,reference_temperature_k,"
        "query_pressure_pa,temperature_k,hsens_ideal_j_kg,cp0_j_kg_k,notes";
    check("frozen reference schema columns are exact", line == expectedHeader);

    std::vector<FrozenRow> rows;
    std::size_t lineNumber = 1;
    while (std::getline(file, line))
    {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split(line);
        const std::string prefix = "reference line " + std::to_string(lineNumber);
        if (fields.size() != 14)
        {
            check(prefix + " has 14 fields", false);
            continue;
        }
        check(prefix + " schema", fields[0] == expectedSchema);
        check(prefix + " backend", fields[3] == "HEOS");
        check(prefix + " CoolProp version", fields[4] == "7.2.0");
        check(prefix + " source API", fields[5] == "PropsSI");
        check(prefix + " enthalpy property",
            fields[6] == "HMASS_IDEALGAS");
        check(prefix + " cp property", fields[7] == "CP0MASS");
        check(prefix + " scope note", fields[13] == expectedNotes);

        try
        {
            const double referenceTemperature = std::stod(fields[8]);
            const double queryPressure = std::stod(fields[9]);
            checkClose(prefix + " reference temperature",
                referenceTemperature, 298.15, 0.0, 1.0e-13);
            checkClose(prefix + " query pressure",
                queryPressure, 1.0e5, 0.0, 0.0);
            rows.push_back(
                {
                    fields[1],
                    std::stod(fields[10]),
                    std::stod(fields[11]),
                    std::stod(fields[12])
                });
        }
        catch (...)
        {
            check(prefix + " numeric fields parse", false);
        }
    }
    return rows;
}

void checkFailClosed(const TabulatedIdealCaloric& caloric)
{
    const auto unknown = caloric.evaluate("N2", 200.0);
    check("unknown species is rejected",
        unknown.status == TabulatedIdealCaloricStatus::UNKNOWN_SPECIES);
    check("unknown species exposes no state", !unknown.state.has_value());

    for (double temperature :
        {
            -1.0,
            0.0,
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN()
        })
    {
        const auto invalid = caloric.evaluate("O2", temperature);
        check("nonphysical/nonfinite temperature fails explicitly",
            invalid.status
                == TabulatedIdealCaloricStatus::INVALID_TEMPERATURE);
        check("invalid temperature exposes no state",
            !invalid.state.has_value());
    }

    for (const auto& range :
        std::vector<std::pair<std::string, std::pair<double, double>>>
        {
            {"O2", {80.0, 300.0}},
            {"C3H8", {80.0, 350.0}}
        })
    {
        const std::string& species = range.first;
        const double minimum = range.second.first;
        const double maximum = range.second.second;
        const auto atMinimum = caloric.evaluate(species, minimum);
        const auto atMaximum = caloric.evaluate(species, maximum);
        check(species + " minimum endpoint is supported", atMinimum.success());
        check(species + " maximum endpoint is supported", atMaximum.success());

        const auto below = caloric.evaluate(
            species, std::nextafter(minimum, -std::numeric_limits<double>::infinity()));
        const auto above = caloric.evaluate(
            species, std::nextafter(maximum, std::numeric_limits<double>::infinity()));
        check(species + " below endpoint fails without clamp",
            below.status == TabulatedIdealCaloricStatus::OUT_OF_RANGE);
        check(species + " below endpoint exposes no state",
            !below.state.has_value());
        check(species + " above endpoint fails without clamp",
            above.status == TabulatedIdealCaloricStatus::OUT_OF_RANGE);
        check(species + " above endpoint exposes no state",
            !above.state.has_value());
    }
}

void checkFrozenRows(
    const TabulatedIdealCaloric& caloric,
    const std::vector<FrozenRow>& rows)
{
    std::map<std::string, std::size_t> counts;
    for (const FrozenRow& row : rows)
    {
        ++counts[row.species];
        const auto result = caloric.evaluate(row.species, row.temperature);
        check(row.species + " frozen node evaluates", result.success());
        if (!result.state) continue;
        checkClose(row.species + " frozen h node",
            result.state->sensibleEnthalpy, row.sensibleEnthalpy,
            3.0e-13, 2.0e-10);
        checkClose(row.species + " frozen cp node",
            result.state->cp, row.cp, 3.0e-13, 2.0e-10);
    }
    check("O2 frozen node count", counts["O2"] == 46);
    check("C3H8 frozen node count", counts["C3H8"] == 70);

    const auto oxygenReference = caloric.evaluate("O2", 298.15);
    const auto propaneReference = caloric.evaluate("C3H8", 298.15);
    check("O2 contains exact 298.15 K reference",
        oxygenReference.success()
        && oxygenReference.state->sensibleEnthalpy == 0.0);
    check("C3H8 contains exact 298.15 K reference",
        propaneReference.success()
        && propaneReference.state->sensibleEnthalpy == 0.0);
    check("O2 entropy uses the exact 298.15 K reference",
        oxygenReference.success()
        && oxygenReference.state->sensibleEntropy == 0.0);
    check("C3H8 entropy uses the exact 298.15 K reference",
        propaneReference.success()
        && propaneReference.state->sensibleEntropy == 0.0);

    // Independently frozen CoolProp 7.2 quad(CP0MASS(T)/T, 298.15, T)
    // values.  SMASS_IDEALGAS on an HEOS state is density-dependent and is
    // not the standard-state entropy integral needed here.
    const auto oxygen130 = caloric.evaluate("O2", 130.0);
    const auto propane200 = caloric.evaluate("C3H8", 200.0);
    checkClose("O2 130 K frozen ideal entropy increment",
        oxygen130.state ? oxygen130.state->sensibleEntropy
                        : std::numeric_limits<double>::quiet_NaN(),
        -7.56420849674296000e+02, 2.0e-9, 2.0e-6);
    checkClose("C3H8 200 K frozen ideal entropy increment",
        propane200.state ? propane200.state->sensibleEntropy
                         : std::numeric_limits<double>::quiet_NaN(),
        -5.77151255702291300e+02, 2.0e-9, 2.0e-6);
}

void checkDenseThermodynamicClosure(
    const TabulatedIdealCaloric& caloric,
    const std::string& species,
    double minimum,
    double maximum)
{
    constexpr int samples = 20000;
    const double spacing = (maximum - minimum)/samples;
    auto previous = caloric.evaluate(species, minimum);
    bool allFinite = previous.success();
    bool allPositive = previous.success() && previous.state->cp > 0.0;
    bool strictlyMonotone = true;
    double maximumDerivativeRelativeError = 0.0;
    double maximumEntropyDerivativeRelativeError = 0.0;

    for (int i = 1; i <= samples; ++i)
    {
        const double temperature = minimum + spacing*i;
        const auto current = caloric.evaluate(species, temperature);
        allFinite = allFinite && current.success()
            && std::isfinite(current.state->sensibleEnthalpy)
            && std::isfinite(current.state->cp)
            && std::isfinite(current.state->sensibleEntropy);
        allPositive = allPositive && current.success() && current.state->cp > 0.0;
        if (!current.success() || !previous.success())
        {
            strictlyMonotone = false;
            previous = current;
            continue;
        }
        strictlyMonotone = strictlyMonotone
            && current.state->sensibleEnthalpy
                > previous.state->sensibleEnthalpy;

        if (i < samples)
        {
            const double delta = 1.0e-3;
            const auto lower = caloric.evaluate(species, temperature - delta);
            const auto upper = caloric.evaluate(species, temperature + delta);
            if (lower.success() && upper.success())
            {
                const double derivative =
                    (upper.state->sensibleEnthalpy
                   - lower.state->sensibleEnthalpy)/(2.0*delta);
                const double relativeError =
                    std::abs(derivative - current.state->cp)
                    /std::max(std::abs(current.state->cp), 1.0);
                maximumDerivativeRelativeError = std::max(
                    maximumDerivativeRelativeError, relativeError);
                const double entropyDerivative =
                    (upper.state->sensibleEntropy
                   - lower.state->sensibleEntropy)/(2.0*delta);
                const double expectedEntropyDerivative =
                    current.state->cp/temperature;
                maximumEntropyDerivativeRelativeError = std::max(
                    maximumEntropyDerivativeRelativeError,
                    std::abs(entropyDerivative - expectedEntropyDerivative)
                        /std::max(std::abs(expectedEntropyDerivative), 1.0));
            }
            else
            {
                allFinite = false;
            }
        }
        previous = current;
    }

    check(species + " dense scan remains finite", allFinite);
    check(species + " dense scan has cp>0", allPositive);
    check(species + " h(T) is strictly monotone", strictlyMonotone);
    check(species + " analytic cp equals dense numerical dh/dT",
        maximumDerivativeRelativeError < 2.0e-8);
    check(species + " entropy derivative equals cp/T",
        maximumEntropyDerivativeRelativeError < 2.0e-8);
    std::cout << "    " << species
              << " max_derivative_relative_error="
              << std::setprecision(8) << maximumDerivativeRelativeError
              << " max_entropy_derivative_relative_error="
              << maximumEntropyDerivativeRelativeError << '\n';
}

void checkKnotContinuity(
    const TabulatedIdealCaloric& caloric,
    const std::vector<FrozenRow>& rows)
{
    std::map<std::string, std::vector<FrozenRow>> bySpecies;
    for (const FrozenRow& row : rows) bySpecies[row.species].push_back(row);

    for (const auto& entry : bySpecies)
    {
        const std::string& species = entry.first;
        const auto& speciesRows = entry.second;
        double maximumSlopeRelativeError = 0.0;
        double maximumCpJumpRelative = 0.0;
        for (std::size_t i = 1; i + 1 < speciesRows.size(); ++i)
        {
            const double temperature = speciesRows[i].temperature;
            const double adjacentSpacing = std::min(
                temperature - speciesRows[i - 1].temperature,
                speciesRows[i + 1].temperature - temperature);
            const double delta = adjacentSpacing*1.0e-6;
            const auto left = caloric.evaluate(species, temperature - delta);
            const auto knot = caloric.evaluate(species, temperature);
            const auto right = caloric.evaluate(species, temperature + delta);
            if (!left.success() || !knot.success() || !right.success())
            {
                check(species + " interior knot neighborhood evaluates", false);
                continue;
            }
            const double leftSlope =
                (knot.state->sensibleEnthalpy
               - left.state->sensibleEnthalpy)/delta;
            const double rightSlope =
                (right.state->sensibleEnthalpy
               - knot.state->sensibleEnthalpy)/delta;
            const double scale = std::max(knot.state->cp, 1.0);
            maximumSlopeRelativeError = std::max(
                maximumSlopeRelativeError,
                std::max(
                    std::abs(leftSlope - knot.state->cp)/scale,
                    std::abs(rightSlope - knot.state->cp)/scale));
            maximumCpJumpRelative = std::max(
                maximumCpJumpRelative,
                std::abs(right.state->cp - left.state->cp)/scale);
        }
        check(species + " h and dh/dT are continuous at every knot",
            maximumSlopeRelativeError < 2.0e-7
            && maximumCpJumpRelative < 2.0e-7);
    }
}

} // End anonymous namespace

int main(int argc, char** argv)
{
    const auto& caloric = tabulatedIdealCaloric();
    check("O2 table is available", caloric.has("O2"));
    check("C3H8 table is available", caloric.has("C3H8"));
    check("unsupported N2 table is absent", !caloric.has("N2"));

    double minimum = 0.0;
    double maximum = 0.0;
    check("O2 range is reported",
        caloric.temperatureRange("O2", minimum, maximum)
        && minimum == 80.0 && maximum == 300.0);
    check("C3H8 range is reported",
        caloric.temperatureRange("C3H8", minimum, maximum)
        && minimum == 80.0 && maximum == 350.0);
    check("unknown range fails",
        !caloric.temperatureRange("N2", minimum, maximum));

    const auto rows = loadReference(locateReference(argc, argv));
    check("frozen reference has 116 data rows", rows.size() == 116);
    checkFailClosed(caloric);
    checkFrozenRows(caloric, rows);
    checkDenseThermodynamicClosure(caloric, "O2", 80.0, 300.0);
    checkDenseThermodynamicClosure(caloric, "C3H8", 80.0, 350.0);
    checkKnotContinuity(caloric, rows);

    check("status names remain explicit",
        std::string(tabulatedIdealCaloricStatusName(
            TabulatedIdealCaloricStatus::OUT_OF_RANGE)) == "OUT_OF_RANGE"
        && std::string(tabulatedIdealCaloricStatusName(
            TabulatedIdealCaloricStatus::INVALID_FROZEN_TABLE))
            == "INVALID_FROZEN_TABLE");

    if (failures != 0)
    {
        std::cerr << "TABULATED_IDEAL_CALORIC_FAIL failures="
                  << failures << '\n';
        return 1;
    }

    std::cout << "TABULATED_IDEAL_CALORIC_PASS\n"
              << "NOTE: CoolProp 7.2.0 ideal calorics only; no liquid/two-phase "
              << "transport or HEOS residual properties.\n";
    return 0;
}
