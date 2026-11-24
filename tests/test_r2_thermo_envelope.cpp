#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "NASA7Thermo.H"

namespace
{

int failures = 0;
int checks = 0;

void check(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

bool contains(double temperature, double minimum, double maximum)
{
    return std::isfinite(temperature)
        && temperature >= minimum
        && temperature <= maximum;
}

} // namespace

int main()
{
    const Foam::NASA7Thermo& thermo = Foam::nasa7Thermo();
    const Foam::SpeciesRegistry& registry = thermo.speciesRegistry();
    const std::vector<std::string> species = {
        "C3H8", "O2", "CO", "CO2", "H2O", "N2"};

    double sharedMinimum = -std::numeric_limits<double>::infinity();
    double sharedMaximum = std::numeric_limits<double>::infinity();
    std::vector<std::string> lowerTailReferenceSpecies;
    for (const std::string& name : species)
    {
        double minimum = 0.0;
        double maximum = 0.0;
        check(
            thermo.temperatureRange(name, minimum, maximum),
            "R2 species has NASA7 data: " + name);
        check(
            minimum > 0.0 && minimum < maximum,
            "R2 species NASA7 range is ordered: " + name);
        sharedMinimum = std::max(sharedMinimum, minimum);
        sharedMaximum = std::min(sharedMaximum, maximum);

        const double reference = thermo.h_mole(
            name, Foam::NASA7Thermo::sensibleReferenceTemperature);
        double expectedReference = reference;
        if (Foam::NASA7Thermo::sensibleReferenceTemperature < minimum)
        {
            lowerTailReferenceSpecies.push_back(name);
            expectedReference = thermo.h_mole(name, minimum)
                + thermo.cp_mole(name, minimum)
                *(Foam::NASA7Thermo::sensibleReferenceTemperature - minimum);
        }
        check(
            std::abs(reference - expectedReference) <= 1.0e-12,
            "RFAM reference enthalpy follows the lower-tail policy: " + name);
        const int global = registry.index(name);
        check(global >= 0, "R2 species is present in SpeciesRegistry: " + name);
        check(
            global >= 0 && std::isfinite(
                reference - registry.formationEnthalpy(global)),
            "RFAM-to-registry reference correction is finite: " + name);
    }

    constexpr double comparisonMaximum = 3400.0;
    constexpr double upperGuard = 100.0;
    check(sharedMinimum == 300.0, "shared R2 NASA7 minimum is 300 K");
    check(sharedMaximum == 3500.0, "shared R2 NASA7 maximum is 3500 K");
    check(
        comparisonMaximum + upperGuard == sharedMaximum,
        "trajectory comparison keeps a 100 K upper guard");
    check(
        lowerTailReferenceSpecies
            == std::vector<std::string>({"C3H8", "N2"}),
        "only C3H8 and N2 use the 298.15 K lower linear tail");
    check(
        contains(comparisonMaximum, sharedMinimum, sharedMaximum),
        "3400 K comparison boundary is inside the shared range");
    check(
        contains(sharedMaximum, sharedMinimum, sharedMaximum),
        "3500 K polynomial boundary is supported");
    check(
        !contains(3500.001, sharedMinimum, sharedMaximum),
        "temperature above 3500 K is outside evidence scope");
    check(
        !contains(std::numeric_limits<double>::quiet_NaN(),
            sharedMinimum, sharedMaximum),
        "non-finite temperature is outside evidence scope");

    double minimum = 0.0;
    double maximum = 0.0;
    check(
        !thermo.temperatureRange("OH", minimum, maximum),
        "unregistered radical thermo is not silently synthesized");

    if (failures != 0)
    {
        std::cerr << "R2_THERMO_ENVELOPE_TEST_FAIL failures=" << failures
                  << " checks=" << checks << '\n';
        return 1;
    }

    std::cout << "R2_THERMO_ENVELOPE_TEST_PASS checks=" << checks
              << " shared_min_k=" << sharedMinimum
              << " shared_max_k=" << sharedMaximum
              << " comparison_max_k=" << comparisonMaximum
              << " production_ready=0\n";
    return 0;
}
