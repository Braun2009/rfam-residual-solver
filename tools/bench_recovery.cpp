// Micro-benchmark of the mixture conserved-state recovery cost per call, with
// and without the table-assisted single-phase fast path.
#include "CoolPropThermoModel.H"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using Foam::thermo::EquilibriumSeed;
using Foam::thermo::EquilibriumRegion;

int main(int argc, char** argv)
{
    std::string table, states;
    int repeat = 200, limit = 400;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--table" && i + 1 < argc) table = argv[++i];
        else if (a == "--states" && i + 1 < argc) states = argv[++i];
        else if (a == "--repeat" && i + 1 < argc) repeat = std::atoi(argv[++i]);
        else if (a == "--limit" && i + 1 < argc) limit = std::atoi(argv[++i]);
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (states.empty()) { std::fprintf(stderr, "need --states\n"); return 2; }
    if (!table.empty()) setenv("RFAM_TETRA_TABLE", table.c_str(), 1);
    auto model = Foam::thermo::makeCoolPropThermoModel();
    std::ifstream in(states);
    std::string line;
    int used = 0;
    double singleSeconds = 0., twoPhaseSeconds = 0.;
    int singleCount = 0, twoPhaseCount = 0, failures = 0;
    while (std::getline(in, line) && used < limit)
    {
        if (line.empty() || line[0] == '#') continue;
        double y, t, p;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf", &y, &t, &p) != 3) continue;
        const auto st = model->initializePTY(p, t, {y, 1. - y});
        if (!st.success()) continue;
        const bool single = st.state->region != EquilibriumRegion::TWO_PHASE;
        const std::vector<double> comp{y, 1. - y};
        const double rho0 = st.state->rho, e0 = st.state->e;
        EquilibriumSeed seed;
        seed.pressure = st.state->p; seed.temperature = st.state->T;
        seed.previousTemperature = st.state->T;
        seed.previousRegion = st.state->region;
        const auto start = std::chrono::steady_clock::now();
        for (int r = 0; r < repeat; ++r)
        {
            const double rho = rho0 * (1. + 1e-7 * r);
            const double e = e0 * (1. + 1e-9 * r);
            if (!model->recoverRhoE(rho, e, comp, &seed).success()) ++failures;
        }
        const double sec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        if (single) { singleSeconds += sec; ++singleCount; }
        else { twoPhaseSeconds += sec; ++twoPhaseCount; }
        ++used;
    }
    const auto ms = [&](double s, int n)
    { return n ? 1000. * s / (double(n) * repeat) : 0.; };
    std::printf("BENCH states=%d single=%d twoPhase=%d failures=%d repeat=%d\n",
                used, singleCount, twoPhaseCount, failures, repeat);
    std::printf("BENCH_PER_CALL_MS single=%.3f twoPhase=%.3f\n",
                ms(singleSeconds, singleCount), ms(twoPhaseSeconds, twoPhaseCount));
    return 0;
}
