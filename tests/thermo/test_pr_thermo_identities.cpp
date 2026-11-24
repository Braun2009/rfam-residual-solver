#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "PengRobinsonEOS.H"
#include "SpeciesCriticalDB.H"
#include "VdWOneFluidMixing.H"

using namespace Foam;

namespace
{

constexpr double R = 8.314462618;
int failures = 0;

void check(
    const std::string& name, bool ok,
    double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

double relativeError(double got, double expected, double floor = 1.0e-14)
{
    return std::abs(got - expected)
        / std::max({std::abs(got), std::abs(expected), floor});
}

double firstDerivative5(
    const std::function<double(double)>& f, double x, double h)
{
    return (f(x - 2.0*h) - 8.0*f(x - h)
          + 8.0*f(x + h) - f(x + 2.0*h))/(12.0*h);
}

double secondDerivative5(
    const std::function<double(double)>& f, double x, double h)
{
    return (-f(x + 2.0*h) + 16.0*f(x + h) - 30.0*f(x)
          + 16.0*f(x - h) - f(x - 2.0*h))/(12.0*h*h);
}

double residualGibbsMolar(
    const PengRobinsonEOS& eos,
    const CriticalProperties& props,
    double p, double T, const std::vector<double>& Y)
{
    const auto root = eos.solveZ(p, T, Y);
    if (!root.valid) return std::numeric_limits<double>::quiet_NaN();

    const double kappa = 0.37464 + 1.54226*props.omega
        - 0.26992*props.omega*props.omega;
    const double alpha = std::pow(
        1.0 + kappa*(1.0 - std::sqrt(T/props.Tc)), 2);
    const double a = 0.45724*R*R*props.Tc*props.Tc/props.Pc*alpha;
    const double b = 0.07780*R*props.Tc/props.Pc;
    const double A = a*p/(R*R*T*T);
    const double B = b*p/(R*T);
    const double sqrt2 = std::sqrt(2.0);
    const double term =
        (root.Z + (1.0 + sqrt2)*B)
        /(root.Z - (sqrt2 - 1.0)*B);
    return R*T*(root.Z - 1.0 - std::log(root.Z - B)
        - A/(2.0*sqrt2*B)*std::log(term));
}

void checkAlphaDerivatives(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& label,
    const std::vector<std::string>& names,
    const std::vector<double>& Y,
    double T)
{
    VdWOneFluidMixing mixing(critical, kij);
    mixing.setSpeciesNames(names);
    mixing.enableMultiSpeciesMixing();

    const double h = 1.0e-3*T;
    const auto a = [&](double temperature)
    {
        return mixing.a_quad(Y, temperature);
    };
    const double daNumeric = firstDerivative5(a, T, h);
    const double d2aNumeric = secondDerivative5(a, T, h);
    const double da = mixing.dadT_quad(Y, T);
    const double d2a = mixing.d2adT2_quad(Y, T);

    check(label + " da/dT identity",
        relativeError(da, daNumeric) < 2.0e-6, da, daNumeric);
    check(label + " d2a/dT2 identity",
        relativeError(d2a, d2aNumeric) < 2.0e-5, d2a, d2aNumeric);
}

} // namespace

int main()
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();

    std::cout << "=== Peng-Robinson thermodynamic identities ===\n";

    checkAlphaDerivatives(
        critical, kij, "O2 130 K", {"O2"}, {1.0}, 130.0);
    checkAlphaDerivatives(
        critical, kij, "C3H8 450 K", {"C3H8"}, {1.0}, 450.0);
    checkAlphaDerivatives(
        critical, kij, "O2/C3H8 400 K",
        {"O2", "C3H8"}, {0.784, 0.216}, 400.0);
    checkAlphaDerivatives(
        critical, kij, "four-species 700 K",
        {"O2", "C3H8", "CO2", "H2O"},
        {0.20, 0.10, 0.40, 0.30}, 700.0);

    const double p = 4.0e6;
    const double T = 130.0;
    const CriticalProperties& oxygen = critical.get("O2");

    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames({"O2"});
    const std::vector<double> Y = {1.0};
    const auto root = eos.solveZ(p, T, Y);
    check("O2 stable root is valid", root.valid);

    if (root.valid)
    {
        const double hDeparture = eos.enthalpyDeparture(p, T, root.Z, Y);
        const double sDeparture = eos.entropyDeparture(p, T, root.Z, Y);
        const double cvDeparture = eos.cvDeparture(p, T, root.Z, Y);
        const double hStep = 1.0e-2;
        const auto gOverT = [&](double temperature)
        {
            return residualGibbsMolar(
                eos, oxygen, p, temperature, Y)/temperature;
        };
        const double hFromGibbs =
            -T*T*firstDerivative5(gOverT, T, hStep)/oxygen.MW;

        check("O2 h departure reference",
            std::abs(hDeparture - (-176268.71)) < 100.0,
            hDeparture, -176268.71);
        check("Gibbs-to-enthalpy identity",
            relativeError(hDeparture, hFromGibbs) < 2.0e-6,
            hDeparture, hFromGibbs);
        const auto gibbs = [&](double temperature)
        {
            return residualGibbsMolar(
                eos, oxygen, p, temperature, Y)/oxygen.MW;
        };
        const double sFromGibbs = -firstDerivative5(gibbs, T, hStep);
        const double sFromHandG = (hDeparture - root.G/oxygen.MW)/T;
        check("Gibbs-to-entropy derivative identity",
            relativeError(sDeparture, sFromGibbs) < 2.0e-6,
            sDeparture, sFromGibbs);
        check("h-g=T*s departure identity",
            relativeError(sDeparture, sFromHandG) < 2.0e-10,
            sDeparture, sFromHandG);
        check("O2 Cv departure reference",
            std::abs(cvDeparture - 204.304) < 1.0,
            cvDeparture, 204.304);

        const double rho = p*oxygen.MW/(root.Z*R*T);
        const auto residualInternalEnergy = [&](double temperature)
        {
            const double pressure = eos.pressure(rho, temperature, Y);
            const double Z = pressure*oxygen.MW/(rho*R*temperature);
            const double hResidual =
                eos.enthalpyDeparture(pressure, temperature, Z, Y);
            return hResidual - (Z - 1.0)*R*temperature/oxygen.MW;
        };
        const double cvFromEnergy =
            firstDerivative5(residualInternalEnergy, T, hStep);
        check("fixed-density Cv identity",
            relativeError(cvDeparture, cvFromEnergy) < 5.0e-6,
            cvDeparture, cvFromEnergy);
    }

    PengRobinsonEOS pseudoCriticalEos(critical, kij);
    const std::vector<double> legacyY = {1.0, 0.0};
    const auto legacyRoot = pseudoCriticalEos.solveZ(p, T, legacyY);
    check("pseudo-critical O2 root is valid", legacyRoot.valid);
    if (legacyRoot.valid)
    {
        const double hDeparture = pseudoCriticalEos.enthalpyDeparture(
            p, T, legacyRoot.Z, legacyY);
        const double cvDeparture = pseudoCriticalEos.cvDeparture(
            p, T, legacyRoot.Z, legacyY);
        check("pseudo-critical h uses corrected derivative/sign",
            std::abs(hDeparture - (-176268.71)) < 100.0,
            hDeparture, -176268.71);
        check("pseudo-critical Cv uses corrected derivative/sign",
            std::abs(cvDeparture - 204.304) < 1.0,
            cvDeparture, 204.304);
    }

    std::cout << (failures == 0
        ? "PR_THERMO_IDENTITIES_PASS\n"
        : "PR_THERMO_IDENTITIES_FAIL\n");
    return failures == 0 ? 0 : 1;
}
