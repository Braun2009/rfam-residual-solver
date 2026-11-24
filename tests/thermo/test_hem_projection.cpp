// Domain-level two-phase projection: (rho, e) inside the density dome but
// with the energy outside the saturated band (numerical mixing noise) must
// recover via the projection branch: keep rho, project e onto the
// saturated endpoints, flag the projection audibly.
//
// Output marker: HEM_PROJECTION_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(const std::string& name, bool ok,
           double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

int run()
{
    static GERGResidualTable gerg;
    static GERGSaturationTable sat;
    if (!gerg.load("tests/reference/gerg_residual_table_v3.csv")
        || !sat.load("tests/reference/gerg_saturation_table_v4.csv"))
    {
        check("tables load", false);
        std::cout << "HEM_PROJECTION_FAIL\n";
        return 1;
    }
    SpeciesRegistry registry;
    const auto names = registry.names();
    std::vector<double> YO2(names.size(), 0.0);
    YO2[registry.index("O2")] = 1.0;

    // Reference two-phase state at 100 K, q=0.4.
    double rhoV = 0.0, rhoL = 0.0, eV = 0.0, eL = 0.0;
    if (!sat.saturationEndpoints(100.0, 1.0, rhoV, rhoL, &eV, &eL))
    {
        check("endpoints", false);
        std::cout << "HEM_PROJECTION_FAIL\n";
        return 1;
    }
    const double rhoMid = 0.5 * (rhoV + rhoL);   // inside the density dome
    // Energy ABOVE the saturated-vapour band (super-heated vapour noise):
    // the exact HEM cannot match rho (x clamps to 1 -> rhoV), so the
    // projection branch must keep rho and project e.
    const double eNoisy = eV * 1.2;

    RealFluidState st;
    st.setSpecies(names);
    st.setGERGTable(&gerg);
    st.setSaturationTable(&sat);
    const auto status = st.updateRhoEY(rhoMid, eNoisy, YO2);
    check("projection recovery SUCCESS + two-phase",
        status == RealFluidState::Status::SUCCESS && st.isTwoPhase());
    if (status == RealFluidState::Status::SUCCESS)
    {
        check("projection flagged (energyProjected)",
            st.energyProjected());
        check("projection keeps rho (mass conserved)",
            std::abs(st.rho() - rhoMid) / rhoMid < 1e-6,
            st.rho(), rhoMid);
        check("projection e within band",
            st.e() >= eL && st.e() <= eV, st.e(), eV);
        check("projection delta finite",
            std::isfinite(st.energyProjectionDelta()),
            st.energyProjectionDelta());
        std::cout << "  [INFO] rho=" << st.rho()
                  << " e_proj=" << st.e()
                  << " e_noisy=" << eNoisy
                  << " delta=" << st.energyProjectionDelta()
                  << " q=" << st.quality() << "\n";
    }

    // Control: exact two-phase recovery must NOT flag a projection.
    {
        const double x = 0.4;
        const double rhoMix = 1.0 / (x / rhoV + (1.0 - x) / rhoL);
        const double eMix = x * eV + (1.0 - x) * eL;
        RealFluidState ctl;
        ctl.setSpecies(names);
        ctl.setGERGTable(&gerg);
        ctl.setSaturationTable(&sat);
        const auto s = ctl.updateRhoEY(rhoMix, eMix, YO2);
        check("control exact HEM: no projection flag",
            s == RealFluidState::Status::SUCCESS && !ctl.energyProjected());
    }

    std::cout << (failures == 0 ? "HEM_PROJECTION_PASS\n"
                                : "HEM_PROJECTION_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
