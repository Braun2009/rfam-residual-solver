// Unit test for TetraTurbulenceOperator on the B benchmark mesh.
// Verifies: wall-distance field, initialization, one advance step,
// eddy-viscosity sanity, and transport hook integration.
// Does NOT require CoolProp flash per face -- uses a fixed lambda transport.
#include "TetraFiniteVolume.H"
#include "TetraTurbulenceOperator.H"
#include "CoolPropThermoModel.H"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace Foam;

int main(int argc, char** argv)
{
    if (argc != 2) { std::cerr << "usage: test_tetra_turbulence <mesh.rfamtet>\n"; return 2; }
    if (!thermo::coolPropThermoAvailable()) return 77;
    try
    {
        TetraMesh mesh; mesh.read(argv[1]);
        const std::size_t n = mesh.cellCount();
        std::cout << "mesh cells=" << n << " faces=" << mesh.faceCount() << "\n";

        // ---- build FV with fixed molecular transport (no CoolProp per face) ---
        auto model = thermo::makeCoolPropThermoModel();
        PrimitiveVars w(2); w.p = 4e6; w.T = 400.; w.Y = {0.784, 0.216};
        w.U = {20., 0., 0.};
        const auto q = model->initializePTY(w.p, w.T, w.Y);
        if (!q.success()) throw std::runtime_error("EOS init failed");

        TetraFiniteVolume fv(mesh, {"O2","C3H8"},
            std::unique_ptr<FluxScheme>(new AUSMFluxAdapter({"O2","C3H8"})));
        for (std::size_t c = 0; c < n; ++c) fv.setPrimitive(c, w, *q.state);

        // Fixed molecular transport: mu=2e-5, lambda=0.04
        fv.setMolecularTransportProvider([](const PrimitiveVars&,
                                            const thermo::EquilibriumState&)
        { return std::make_pair(2e-5, 0.04); });

        // Prime lastMu/lastLambda with one residual evaluation
        const TetraBoundaryProvider boundary = [&](std::size_t, std::size_t,
            PrimitiveVars& gw, thermo::EquilibriumState& gq)
        { gw = w; gq = *q.state; };
        std::vector<ConservedVars> res;
        fv.computeResidual(boundary, res);

        // ---- construct turbulence operator ----------------------------
        TetraTurbulenceOperator::Settings ts;
        ts.intensity    = 0.05;
        ts.lengthScale  = 1e-3;
        ts.wallZone     = 0;
        ts.prandtlT     = 0.85;
        ts.wallOmegaFactor = 60.0;

        TetraTurbulenceOperator turbulence(fv, ts);

        // Check wall-distance field
        const auto& y = turbulence.wallDistanceField();
        if (y.size() != n) throw std::runtime_error("wall distance size mismatch");
        const double yMin = *std::min_element(y.begin(), y.end());
        const double yMax = *std::max_element(y.begin(), y.end());
        std::cout << "wall_distance: y_min=" << yMin << " y_max=" << yMax << " m\n";
        if (!(yMin > 0.))   throw std::runtime_error("wall distance has zero or negative entry");
        if (!(yMax < 1.0))  throw std::runtime_error("wall distance max > 1 m (unreasonable)");
        if (!(yMax > yMin)) throw std::runtime_error("wall distance field is constant");

        // Check all k/omega finite after construction
        const auto& kf = turbulence.kField();
        const auto& of = turbulence.omegaField();
        for (std::size_t c = 0; c < n; ++c)
        {
            if (!std::isfinite(kf[c]) || kf[c] <= 0.)
                throw std::runtime_error("k not finite/positive at cell " + std::to_string(c));
            if (!std::isfinite(of[c]) || of[c] <= 0.)
                throw std::runtime_error("omega not finite/positive at cell " + std::to_string(c));
        }

        // ---- initialize from flow field ---
        turbulence.initialize();

        // ---- wire transport hook and advance one step ------------------
        TetraTurbulenceOperator* raw = &turbulence;
        fv.setTurbulentTransportProvider(
            [raw](std::size_t c, const PrimitiveVars& wc, const thermo::EquilibriumState& qc)
            { return raw->transportAddition(c, wc, qc); });
        if (!fv.hasTurbulentTransport())
            throw std::runtime_error("turbulent transport hook not registered");

        // Compute dt vector (simple fixed dt)
        std::vector<double> dtCell(n, 1e-7);

        turbulence.advance(dtCell);

        // Check results
        const double nutMax = turbulence.maxNut();
        const std::size_t wallCells = turbulence.wallAdjacentCells();
        std::cout << "nut_max=" << nutMax << " turb_wall_cells=" << wallCells << "\n";

        if (!std::isfinite(nutMax))    throw std::runtime_error("nut_max is not finite");
        if (!(nutMax > 0.))            throw std::runtime_error("nut_max is zero (no turbulence produced)");
        if (!(nutMax < 1.0))           throw std::runtime_error("nut_max > 1 m2/s (unreasonable)");
        if (wallCells == 0)            throw std::runtime_error("no wall-adjacent cells found");
        if (wallCells >= n)            throw std::runtime_error("all cells flagged as wall-adjacent");

        // Wall cells must have k ~ k_min
        for (std::size_t c = 0; c < n; ++c)
        {
            if (!std::isfinite(kf[c]) || kf[c] <= 0.)
                throw std::runtime_error("k non-positive after advance at cell " + std::to_string(c));
            if (!std::isfinite(of[c]) || of[c] <= 0.)
                throw std::runtime_error("omega non-positive after advance at cell " + std::to_string(c));
        }

        // transportAddition should return positive mu_t, lambda_t for interior cells
        std::size_t positiveMuT = 0;
        for (std::size_t c = 0; c < n; ++c)
        {
            const auto [muT, lamT] = turbulence.transportAddition(c, w, *q.state);
            if (muT > 0. && lamT > 0.) ++positiveMuT;
            if (muT < 0. || lamT < 0.)
                throw std::runtime_error("negative turbulent transport at cell " + std::to_string(c));
        }
        if (positiveMuT == 0)
            throw std::runtime_error("all turbulent mu_t are zero");

        std::cout << "TETRA_TURBULENCE_PASS"
                  << " cells=" << n
                  << " wall_cells=" << wallCells
                  << " y_min=" << yMin
                  << " y_max=" << yMax
                  << " nut_max=" << nutMax
                  << " positive_muT_cells=" << positiveMuT
                  << "\n";
        return 0;
    }
    catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << "\n"; return 1; }
}
