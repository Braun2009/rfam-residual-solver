#ifdef RFAM_HAVE_TEQP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "NASA7Thermo.H"
#include "RealFluidState.H"
#include "TeqpMultifluidThermo.H"

using namespace Foam;

namespace
{

bool close(double actual, double expected, double relative, double absolute)
{
    return std::isfinite(actual) && std::isfinite(expected)
        && std::abs(actual-expected)
            <= absolute + relative*std::max(std::abs(actual),std::abs(expected));
}

int failures = 0;

void check(bool condition, const std::string& label)
{
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << label << '\n';
    }
}

} // namespace

int main()
{
    // A representative hot combustion-product state from the A406 iter=1000
    // sector field (cell 2,2,5).  Y is RFAM mass basis, whereas teqp's
    // multifluid derivative API takes mole fractions.
    const std::vector<double> massFractions{
        0.0697176519051,
        0.193696064,
        5.91804843829e-19,
        2.30017422924e-13,
        9.77962086206e-18,
        0.0186380925168,
        0.254318153846,
        0.00366280543697,
        0.459967232294};
    constexpr double rhoMass = 6.84183623359;
    constexpr double temperature = 2119.30972409;

    check(RealFluidState::teqpAvailable(),
        "teqp multifluid backend is constructible");
    if (failures != 0)
    {
        std::cout << "TEQP_MULTIFLUID_MASS_BASIS_FAIL\n";
        return 1;
    }

    const TeqpMultifluidThermo& thermo = RealFluidState::teqpInstance();
    std::vector<double> moleFractions;
    double mixtureMolecularWeight =
        std::numeric_limits<double>::quiet_NaN();
    check(TeqpMultifluidThermo::massToMoleFractions(
            massFractions,moleFractions,&mixtureMolecularWeight),
        "mass-to-mole conversion succeeds");
    check(moleFractions.size() == massFractions.size()
            && close(mixtureMolecularWeight,0.029226986102815624,
                2.0e-12,1.0e-15),
        "mass-basis molecular weight is exact");
    check(close(moleFractions[0],0.06367853930000902,1.0e-12,1.0e-14)
            && close(moleFractions[6],0.4125915971418405,1.0e-12,1.0e-14)
            && close(moleFractions[8],0.3054671356416093,1.0e-12,1.0e-14),
        "hot-product mole fractions are not treated as mass fractions");
    if (failures != 0)
    {
        std::cout << "TEQP_MULTIFLUID_MASS_BASIS_FAIL\n";
        return 1;
    }

    const double rhoMolar = rhoMass/mixtureMolecularWeight;
    const double directAr01 = teqp::TDXDerivatives<
        TeqpMultifluidThermo::ModelT,double,std::vector<double>>
            ::template get_Arxy<0,1>(
                thermo.model(),temperature,rhoMolar,moleFractions);
    const double directPressure = rhoMolar*thermo.model().R(moleFractions)
        *temperature*(1.0+directAr01);
    const double wrappedPressure = thermo.pressure(
        temperature,rhoMolar,massFractions);
    check(close(wrappedPressure,directPressure,2.0e-13,1.0e-5),
        "wrapped pressure passes mole fractions to teqp");

    const double directAr10 = teqp::TDXDerivatives<
        TeqpMultifluidThermo::ModelT,double,std::vector<double>>
            ::template get_Arxy<1,0>(
                thermo.model(),temperature,rhoMolar,moleFractions);
    double massSum = 0.0;
    for (double yi : massFractions) massSum += yi;
    std::vector<double> normalizedMassFractions = massFractions;
    for (double& yi : normalizedMassFractions) yi /= massSum;
    const double directEnergy = nasa7Thermo().hsens_mass_mix(
        TeqpMultifluidThermo::speciesNames(),normalizedMassFractions,
        temperature)
        - TeqpMultifluidThermo::R_UNIVERSAL*temperature
            /mixtureMolecularWeight
        + TeqpMultifluidThermo::R_UNIVERSAL*temperature*directAr10
            /mixtureMolecularWeight;
    const double wrappedEnergy = thermo.eSens(
        temperature,rhoMolar,massFractions);
    check(close(wrappedEnergy,directEnergy,2.0e-13,1.0e-5),
        "wrapped sensible energy uses mass ideal term and mole residual term");

    double recoveredTemperature = std::numeric_limits<double>::quiet_NaN();
    double recoveredPressure = std::numeric_limits<double>::quiet_NaN();
    check(thermo.recoverRhoEY(
            rhoMass,wrappedEnergy,massFractions,
            recoveredTemperature,recoveredPressure)
            && close(recoveredTemperature,temperature,2.0e-8,1.0e-6)
            && close(recoveredPressure,wrappedPressure,2.0e-8,1.0e-2),
        "multifluid rho-e recovery round-trips a hot mixed state");

    // The diagnostic chemistry wrapper intentionally does not own a GERG
    // residual-table pointer.  A mixture must still reach teqp on A406 rather
    // than silently falling through to PR solely because that optional table
    // is absent.
    RealFluidState state;
    state.setSpecies(TeqpMultifluidThermo::speciesNames());
    const RealFluidState::Status status = state.updateRhoEY(
        rhoMass,wrappedEnergy,massFractions,temperature);
    check(status == RealFluidState::Status::SUCCESS
            && close(state.T(),temperature,2.0e-8,1.0e-6)
            && close(state.p(),wrappedPressure,2.0e-8,1.0e-2),
        "no-table mixed rho-e recovery routes through teqp");

    std::cout << (failures == 0
        ? "TEQP_MULTIFLUID_MASS_BASIS_PASS\n"
        : "TEQP_MULTIFLUID_MASS_BASIS_FAIL\n");
    return failures == 0 ? 0 : 1;
}

#else

int main()
{
    return 0;
}

#endif
