/*---------------------------------------------------------------------------*\
    test_cantera_compare - RFAM kinetics production-readiness gate.

    Compares RFAM TwoStepCO (2-step mechanism, C3H8/O2) against
    Cantera GRI-Mech 3.0 reference values from scripts/cantera_reference.py.

    Metrics:
      1. Ignition delay tau_ign(1200K, 4MPa, phi=1): RFAM legacy fixed-Cv
         diagnostic vs Cantera 144.53 us constant-pressure reference.
      2. Post-combustion T: RFAM legacy fixed-Cv diagnostic vs Cantera
         Tad=3615.6K constant-pressure equilibrium at 4MPa. The closures are
         not directly comparable, and the gap is recorded honestly.
      3. Species plausibility: fuel consumed, CO2/H2O produced, element C/H/O
         conserved.

    HONEST BOUNDARY: TwoStepCO is a global 2-step mechanism that is not
    calibrated for this ignition-delay or adiabatic-state reference. The RFAM
    fixed-Cv diagnostic is not thermodynamically equivalent to the Cantera
    constant-pressure calculation. Su comparison is NOT done here because RFAM
    lacks diffusion-coupled species transport for a
    proper 1D flame (that's Stage 4 work after 2D domain is established).

    Build (requires SUNDIALS):
      see local_test.sh: run cantera_compare
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "src/thermo/ThermoDataProvenance.H"
#include "src/thermo/SpeciesRegistry.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/combustion/FiniteRateReactor.H"

#ifndef RFAM_CANTERA_GENERATOR_SHA256
#error "RFAM_CANTERA_GENERATOR_SHA256 must be bound at build time"
#endif

#ifndef RFAM_CANTERA_REFERENCE_SHA256
#error "RFAM_CANTERA_REFERENCE_SHA256 must be bound at build time"
#endif

#ifndef RFAM_CANTERA_MECHANISM_SHA256
#error "RFAM_CANTERA_MECHANISM_SHA256 must be bound at build time"
#endif

using namespace Foam;
using namespace std;

namespace
{

struct CanteraReference
{
    string canteraVersion;
    string generatorSha256;
    string mechanismName;
    string mechanismSha256;
    string fuel;
    string oxidizer;
    double ignitionInitialTemperature = 0.0;
    double ignitionPressure = 0.0;
    double ignitionPhi = 0.0;
    double ignitionDeltaTemperature = 0.0;
    double ignitionDelayUs = 0.0;
    double adiabaticInitialTemperature = 0.0;
    double adiabaticPressure = 0.0;
    double adiabaticPhi = 0.0;
    double adiabaticTemperature = 0.0;
    double flameSpeedInitialTemperature = 0.0;
    double flameSpeedPressure = 0.0;
    double flameSpeedPhi = 0.0;
    double laminarFlameSpeed = 0.0;
};

string readBinaryFile(const string& path)
{
    ifstream input(path, ios::binary);
    if (!input.good())
        throw runtime_error("provenance-open-failed:" + path);
    return string(
        istreambuf_iterator<char>(input), istreambuf_iterator<char>());
}

bool validSha256(const string& value)
{
    if (value.size() != 64) return false;
    for (unsigned char character : value)
        if (!isdigit(character)
            && !(character >= 'a' && character <= 'f'))
            return false;
    return true;
}

double parseFiniteDouble(const string& value, const string& key)
{
    size_t consumed = 0;
    double parsed = 0.0;
    try
    {
        parsed = stod(value, &consumed);
    }
    catch (const exception&)
    {
        throw runtime_error("reference-invalid-number:" + key);
    }
    if (consumed != value.size() || !isfinite(parsed))
        throw runtime_error("reference-invalid-number:" + key);
    return parsed;
}

CanteraReference parseReference(const string& contents)
{
    map<string, string> values;
    istringstream input(contents);
    string line;
    size_t lineNumber = 0;
    while (getline(input, line))
    {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty())
            throw runtime_error(
                "reference-empty-line:" + to_string(lineNumber));
        const size_t separator = line.find('=');
        if (separator == string::npos || separator == 0
            || separator + 1 >= line.size()
            || line.find('=', separator + 1) != string::npos)
            throw runtime_error(
                "reference-invalid-line:" + to_string(lineNumber));
        const string key = line.substr(0, separator);
        const string value = line.substr(separator + 1);
        if (!values.emplace(key, value).second)
            throw runtime_error("reference-duplicate-key:" + key);
    }

    static const vector<string> requiredKeys = {
        "schema", "cantera_version", "generator_sha256",
        "mechanism_name", "mechanism_sha256", "fuel", "oxidizer",
        "ignition_initial_temperature_k", "ignition_pressure_pa",
        "ignition_phi", "ignition_delta_temperature_k",
        "ignition_delay_us", "adiabatic_initial_temperature_k",
        "adiabatic_pressure_pa", "adiabatic_phi",
        "adiabatic_temperature_k", "flame_speed_initial_temperature_k",
        "flame_speed_pressure_pa", "flame_speed_phi",
        "laminar_flame_speed_m_per_s"};
    if (values.size() != requiredKeys.size())
        throw runtime_error("reference-key-count-mismatch");
    for (const string& key : requiredKeys)
        if (values.find(key) == values.end())
            throw runtime_error("reference-missing-key:" + key);

    if (values.at("schema") != "rfam.cantera-kinetics-reference/v1")
        throw runtime_error("reference-schema-mismatch");

    CanteraReference reference;
    reference.canteraVersion = values.at("cantera_version");
    reference.generatorSha256 = values.at("generator_sha256");
    reference.mechanismName = values.at("mechanism_name");
    reference.mechanismSha256 = values.at("mechanism_sha256");
    reference.fuel = values.at("fuel");
    reference.oxidizer = values.at("oxidizer");
    reference.ignitionInitialTemperature = parseFiniteDouble(
        values.at("ignition_initial_temperature_k"),
        "ignition_initial_temperature_k");
    reference.ignitionPressure = parseFiniteDouble(
        values.at("ignition_pressure_pa"), "ignition_pressure_pa");
    reference.ignitionPhi = parseFiniteDouble(
        values.at("ignition_phi"), "ignition_phi");
    reference.ignitionDeltaTemperature = parseFiniteDouble(
        values.at("ignition_delta_temperature_k"),
        "ignition_delta_temperature_k");
    reference.ignitionDelayUs = parseFiniteDouble(
        values.at("ignition_delay_us"), "ignition_delay_us");
    reference.adiabaticInitialTemperature = parseFiniteDouble(
        values.at("adiabatic_initial_temperature_k"),
        "adiabatic_initial_temperature_k");
    reference.adiabaticPressure = parseFiniteDouble(
        values.at("adiabatic_pressure_pa"), "adiabatic_pressure_pa");
    reference.adiabaticPhi = parseFiniteDouble(
        values.at("adiabatic_phi"), "adiabatic_phi");
    reference.adiabaticTemperature = parseFiniteDouble(
        values.at("adiabatic_temperature_k"),
        "adiabatic_temperature_k");
    reference.flameSpeedInitialTemperature = parseFiniteDouble(
        values.at("flame_speed_initial_temperature_k"),
        "flame_speed_initial_temperature_k");
    reference.flameSpeedPressure = parseFiniteDouble(
        values.at("flame_speed_pressure_pa"), "flame_speed_pressure_pa");
    reference.flameSpeedPhi = parseFiniteDouble(
        values.at("flame_speed_phi"), "flame_speed_phi");
    reference.laminarFlameSpeed = parseFiniteDouble(
        values.at("laminar_flame_speed_m_per_s"),
        "laminar_flame_speed_m_per_s");
    return reference;
}

bool closeReferenceValue(double value, double expected)
{
    return abs(value - expected)
        <= 1e-12 * max(1.0, abs(expected));
}

void validateReferenceContract(const CanteraReference& reference)
{
    static const string pinnedMechanismSha256 =
        "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345";
    if (reference.canteraVersion != "3.1.0"
        || reference.mechanismName != "gri30.yaml"
        || reference.mechanismSha256 != pinnedMechanismSha256
        || !validSha256(reference.mechanismSha256)
        || reference.fuel != "C3H8" || reference.oxidizer != "O2"
        || !closeReferenceValue(reference.ignitionInitialTemperature, 1200.0)
        || !closeReferenceValue(reference.ignitionPressure, 4.0e6)
        || !closeReferenceValue(reference.ignitionPhi, 1.0)
        || !closeReferenceValue(reference.ignitionDeltaTemperature, 400.0)
        || !closeReferenceValue(reference.adiabaticInitialTemperature, 300.0)
        || !closeReferenceValue(reference.adiabaticPressure, 4.0e6)
        || !closeReferenceValue(reference.adiabaticPhi, 1.0)
        || !closeReferenceValue(reference.flameSpeedInitialTemperature, 300.0)
        || !closeReferenceValue(reference.flameSpeedPressure, 101325.0)
        || !closeReferenceValue(reference.flameSpeedPhi, 1.0)
        || reference.ignitionDelayUs <= 0.0
        || reference.adiabaticTemperature <= 0.0
        || reference.laminarFlameSpeed <= 0.0)
        throw runtime_error("reference-contract-mismatch");
}

} // namespace

static int g_pass = 0, g_fail = 0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if (!ok) cout << "   (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

static constexpr double MIN_IGNITION_DELAY_RATIO = 0.1;
static constexpr double MAX_IGNITION_DELAY_RATIO = 10.0;
static constexpr double MAX_TEMPERATURE_RELATIVE_ERROR = 0.10;

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        cerr << "CANTERA_KINETICS_GATE_FAIL provenance-arguments-required\n";
        return 2;
    }

    CanteraReference reference;
    string generatorSha256;
    string referenceSha256;
    string mechanismSha256;
    try
    {
        const string generatorContents = readBinaryFile(argv[1]);
        const string referenceContents = readBinaryFile(argv[2]);
        const string mechanismContents = readBinaryFile(argv[3]);
        generatorSha256 = thermoProvenance::sha256Hex(generatorContents);
        referenceSha256 = thermoProvenance::sha256Hex(referenceContents);
        mechanismSha256 = thermoProvenance::sha256Hex(mechanismContents);
        if (!validSha256(RFAM_CANTERA_GENERATOR_SHA256)
            || generatorSha256 != RFAM_CANTERA_GENERATOR_SHA256)
            throw runtime_error("generator-sha256-mismatch");
        if (!validSha256(RFAM_CANTERA_REFERENCE_SHA256)
            || referenceSha256 != RFAM_CANTERA_REFERENCE_SHA256)
            throw runtime_error("reference-sha256-mismatch");
        if (!validSha256(RFAM_CANTERA_MECHANISM_SHA256)
            || mechanismSha256 != RFAM_CANTERA_MECHANISM_SHA256)
            throw runtime_error("mechanism-sha256-mismatch");
        reference = parseReference(referenceContents);
        if (!validSha256(reference.generatorSha256)
            || reference.generatorSha256 != generatorSha256)
            throw runtime_error("reference-generator-sha256-mismatch");
        validateReferenceContract(reference);
        if (reference.mechanismSha256 != mechanismSha256)
            throw runtime_error("reference-mechanism-sha256-mismatch");
    }
    catch (const exception& error)
    {
        cerr << "CANTERA_KINETICS_GATE_FAIL " << error.what() << "\n";
        return 2;
    }

    cout << "=== RFAM vs Cantera reference comparison ===\n";
    cout << "    Cantera: " << reference.canteraVersion
         << " / " << reference.mechanismName << "\n";
    cout << "    RFAM:    TwoStepCO (2-step global, C3H8+3.5O2->3CO+4H2O, CO+0.5O2->CO2)\n\n";
    cout << "    generator path: " << argv[1] << "\n";
    cout << "    reference path: " << argv[2] << "\n";
    cout << "    mechanism path: " << argv[3] << "\n";
    cout << "    generator SHA256: " << generatorSha256 << "\n";
    cout << "    reference SHA256: " << referenceSha256 << "\n";
    cout << "    mechanism SHA256: " << mechanismSha256 << "\n\n";

    SpeciesRegistry reg;
    const int nS = reg.size();

    // ---- IC: stoich C3H8/O2 at 1200K, 4MPa (= RFAM A-box IC) ----
    const double p0 = reference.ignitionPressure;
    const double T0 = reference.ignitionInitialTemperature;
    vector<double> Y0(nS, 0.0);
    Y0[reg.index("C3H8")] = 0.216;
    Y0[reg.index("O2")]   = 0.784;
    // rho from ideal gas approximation (good enough for 0-D IC)
    const double R_mix = 8.314462618 * (Y0[reg.index("C3H8")] / 0.0440956
                                       + Y0[reg.index("O2")]   / 0.0319988);
    const double rho0 = p0 / (R_mix * T0);

    // ---- Run RFAM 0-D const-p reactor ----
    TwoStepCO chem(1e10, 0.0, 125520.0, 2e12, 0.0, 167360.0);
    FiniteRateReactor reactor(
        &chem,reg,FiniteRateReactor::LegacyFixedCvDiagnostic(1500.0));
    // Explicit legacy diagnostic retained only for comparison with the frozen
    // historical Cantera audit; CFD operators do not use this closure.

    double T = T0;
    const double dT_rise = reference.ignitionDeltaTemperature;
    vector<double> Y = Y0;
    double dE = 0.0;
    double tau_ign = -1.0;

    const double dt = 1.0e-7;   // 100 ns explicit steps
    const int maxSteps = 100000;
    double t = 0.0;

    for (int step = 0; step < maxSteps; ++step)
    {
        reactor.reactLegacyFixedCv(rho0, p0, T, Y, dt, dE);
        t += dt;
        if (tau_ign < 0.0 && T - T0 > dT_rise)
            tau_ign = t;
        if (Y[reg.index("C3H8")] < 1e-6 && t > 1e-6) break;  // fully burnt
    }

    double Tfinal = T;
    double tau_us = tau_ign * 1e6;

    cout << fixed << setprecision(1);
    cout << "  IC: T=" << T0 << "K  p=" << p0/1e6 << "MPa  phi=1  C3H8/O2\n";
    cout << "  RFAM  tau_ign = " << setprecision(2) << tau_us << " µs"
         << "  (Cantera=" << reference.ignitionDelayUs << " µs)\n";
    cout << "  RFAM  T_final = " << setprecision(1) << Tfinal << "K"
         << "  (Cantera Tad=" << reference.adiabaticTemperature
         << "K const-p)\n\n";

    // ---- Metric 1: tau_ign detection ----
    // TwoStepCO is a GLOBAL 2-step mechanism with no ignition-delay calibration.
    // At 1200K the Arrhenius rate is so large that ignition fires within the
    // first time-step (<<1 µs). GRI-Mech gives 144.53 µs because it resolves
    // the low-T chain-branching induction chemistry absent in the global mech.
    // CONSEQUENCE: tau_ign comparison is QUALITATIVE only. We only verify that
    // (a) ignition was detected at all, and (b) T rose significantly.
    // Quantitative ignition delay requires a skeletal/detailed mechanism.
    // [[technical-debt]]: add n-heptane or USC-Mech II skeletal for tau_ign validation.
    double ignitionRatio = 0.0;
    bool ignitionWithinTolerance = false;
    if (tau_ign > 0.0)
    {
        ignitionRatio = tau_us / reference.ignitionDelayUs;
        ignitionWithinTolerance =
            ignitionRatio >= MIN_IGNITION_DELAY_RATIO
            && ignitionRatio <= MAX_IGNITION_DELAY_RATIO;
        cout << "  ratio = " << setprecision(4) << ignitionRatio << "x Cantera"
             << "  [NOTE: TwoStepCO is not calibrated for ignition delay —"
             << " ratio << 1 is expected]\n";
        check("tau_ign: ignition detected (T rose > 400K)", tau_ign > 0.0);
        check("production ignition-delay tolerance rejects current candidate",
              !ignitionWithinTolerance, ignitionRatio, MIN_IGNITION_DELAY_RATIO);
    }
    else
    {
        check("tau_ign: ignition detected", false);
    }

    // ---- Metric 2: T_final vs Cantera Tad ----
    // RFAM (ideal-Cv, no dissociation): T_final >> Cantera Tad because:
    //   (a) ideal Cv=1500 instead of real mixture Cv
    //   (b) no dissociation (A+B = reactants can't recombine endothermically)
    // Expected: T_final 1.5-2.5x Cantera Tad. We just document the gap.
    double T_ratio = Tfinal / reference.adiabaticTemperature;
    const double temperatureRelativeError =
        fabs(Tfinal - reference.adiabaticTemperature)
        / reference.adiabaticTemperature;
    const bool temperatureWithinTolerance =
        temperatureRelativeError <= MAX_TEMPERATURE_RELATIVE_ERROR;
    cout << "  T_final/Cantera_Tad = " << setprecision(2) << T_ratio
         << "x  [expected ~2x: no dissociation + ideal-Cv closure]\n";
    check("T_final > 0 and finite", Tfinal > 0.0 && isfinite(Tfinal), Tfinal);
    check("T_final > T0 (ignition occurred)", Tfinal > T0 * 1.5, Tfinal, T0);
    check("production adiabatic-temperature tolerance rejects current candidate",
          !temperatureWithinTolerance,
          temperatureRelativeError, MAX_TEMPERATURE_RELATIVE_ERROR);

    // ---- Metric 3: species plausibility ----
    check("C3H8 consumed  (Y < 1e-4)", Y[reg.index("C3H8")] < 1e-4,
          Y[reg.index("C3H8")]);
    // With reversible R2, CO2 partially dissociates back to CO at high T.
    // Check total carbon-containing products (CO+CO2) instead of CO2 alone.
    // Stoich C3H8/O2: 3 mol CO per mol C3H8 → mass fraction ~3×28/44×0.216 ≈ 0.41
    double Y_CO_CO2 = Y[reg.index("CO")] + Y[reg.index("CO2")];
    check("CO+CO2 carbon products > 0.3",  Y_CO_CO2  > 0.3, Y_CO_CO2);
    check("H2O produced   (Y > 0.1)",  Y[reg.index("H2O")]  > 0.1,
          Y[reg.index("H2O")]);
    // element conservation (total atoms)
    double C0=0, C1=0, H0=0, H1=0, O0=0, O1=0;
    for (int g=0;g<nS;++g) {
        C0 += Y0[g]/reg.molecularWeight(g)*reg.nC(g);
        C1 += Y[g] /reg.molecularWeight(g)*reg.nC(g);
        H0 += Y0[g]/reg.molecularWeight(g)*reg.nH(g);
        H1 += Y[g] /reg.molecularWeight(g)*reg.nH(g);
        O0 += Y0[g]/reg.molecularWeight(g)*reg.nO(g);
        O1 += Y[g] /reg.molecularWeight(g)*reg.nO(g);
    }
    double ref_atom = C0+H0+O0;
    check("element C conserved (< 1e-4)", fabs(C1-C0)/(ref_atom+1e-30) < 1e-4);
    check("element H conserved (< 1e-4)", fabs(H1-H0)/(ref_atom+1e-30) < 1e-4);
    check("element O conserved (< 1e-4)", fabs(O1-O0)/(ref_atom+1e-30) < 1e-4);

    // Production readiness requires comparable energy closures and both
    // quantitative reference checks. The present legacy fixed-Cv diagnostic is
    // intentionally not comparable to the constant-pressure reference.
    const bool referenceContractsComparable = false;
    const bool productionKineticsReady = referenceContractsComparable
        && ignitionWithinTolerance && temperatureWithinTolerance;
    cout << scientific << setprecision(8)
         << "  production_ready=" << productionKineticsReady
         << " reference_contracts_comparable="
         << referenceContractsComparable
         << " ignition_ratio=" << ignitionRatio
         << " temperature_relative_error=" << temperatureRelativeError
         << "\n" << defaultfloat;
    check("production kinetics remains fail-closed", !productionKineticsReady);

    cout << "\n  NOTE: Su (laminar flame speed) comparison deferred to Stage 4.\n";
    cout << "  Cantera Su(1atm,300K,phi=1,C3H8/O2) = "
         << reference.laminarFlameSpeed << " m/s  [stored for future]\n";
    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 && !productionKineticsReady
        ? "CANTERA_KINETICS_GATE_BLOCKED"
        : "CANTERA_KINETICS_GATE_FAIL") << "\n";
    return g_fail==0 && !productionKineticsReady ? 0 : 1;
}
