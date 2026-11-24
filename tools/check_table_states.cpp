// Compare the table-interpolated equilibrium sound speed with the direct EOS
// isentrope at a list of (Y_O2, T, p) states (CSV). Verification tool.
#include "CoolPropThermoModel.H"
#include "TableAcoustic.H"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

int main(int argc, char** argv)
{
    std::string table;
    std::string states;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--table" && i + 1 < argc) table = argv[++i];
        else if (a == "--states" && i + 1 < argc) states = argv[++i];
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (table.empty() || states.empty())
    {
        std::fprintf(stderr, "usage: check_table_states --table T.bin --states states.csv\n");
        return 2;
    }
    setenv("RFAM_TETRA_TABLE", table.c_str(), 1);
    const auto* tab = Foam::thermo::tableAcoustic();
    if (!tab) { std::fprintf(stderr, "table load failed: %s\n", table.c_str()); return 3; }
    auto model = Foam::thermo::makeCoolPropThermoModel();
    std::ifstream in(states);
    std::string line;
    int count = 0, directFail = 0, tableMiss = 0;
    double worst = 0., sum = 0., worstAbs = 0.;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#') continue;
        double y, t, p;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf", &y, &t, &p) != 3) continue;
        double aTable = -1.;
        const bool haveTable = tab->lookup(y, t, p, aTable);
        const auto r = model->initializePTY(p, t, {y, 1. - y});
        ++count;
        if (!r.success() || !r.state) { ++directFail; continue; }
        const double aDirect = r.state->acoustic.soundSpeed();
        if (!haveTable) { ++tableMiss; continue; }
        const double rel = std::fabs(std::sqrt(aTable) - aDirect) / aDirect;
        if (rel > worst) { worst = rel; worstAbs = std::sqrt(aTable) - aDirect; }
        sum += rel;
        std::printf("STATE y=%.6g T=%.4g p=%.6g direct=%.6g table=%.6g rel=%.3e\n",
                    y, t, p, aDirect, std::sqrt(aTable), rel);
    }
    std::printf("SUMMARY n=%d directFail=%d tableMiss=%d meanRel=%.3e worstRel=%.3e worstAbs=%.3e\n",
                count, directFail, tableMiss, count ? sum / count : 0., worst, worstAbs);
    return 0;
}
