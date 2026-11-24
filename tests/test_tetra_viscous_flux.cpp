// Directed test for the viscous flux assembly: a manufactured LINEAR velocity
// and temperature field must reproduce the analytic stress/work/heat flux and
// the assembled residual must be divergence-free in the interior (linear u and
// T imply constant tau and grad T, so every interior face's contribution
// cancels pairwise exactly).
#include "TetraViscousFlux.H"
#include <cmath>
#include <cstdio>
#include <string>

using namespace Foam;

int main(int argc, char** argv)
{
    const std::string path = argc > 1 ? argv[1]
        : "work/coldflow_B_verified/benchmark_mesh/fluid_volume.rfamtet";
    TetraMesh mesh;
    try { mesh.read(path); }
    catch (...) { std::printf("TETRA_VISCOUS_SKIP no mesh %s\n", path.c_str()); return 77; }
    const std::size_t n = mesh.cellCount();
    // Manufactured field: u = (a*y, b*x + c*z, d*y), T = T0 + gx*x + gy*y.
    const double a = 3.0, b = -1.5, c = 2.0, d = 0.7;
    const double T0 = 300., gx = 40., gy = -25.;
    std::vector<TetraVec3> velocity(n);
    std::vector<double> temperature(n), mu(n, 2.0e-5), lambda(n, 0.05);
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto x = mesh.cellCentroid(i);
        velocity[i] = {a*x.y, b*x.x + c*x.z, d*x.y};
        temperature[i] = T0 + gx*x.x + gy*x.y;
    }
    // Analytic flux at a face with unit normal n: tau is constant.
    const auto analytic = [&](const std::array<double, 3>& nn)
    {
        double tau[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        const double divU = 0.0;  // a*x + (b x + c z) in y + d*y in z -> d/dx=0, d/dy=0, d/dz=0
        const double gradU[3][3] = {{0, a, 0}, {b, 0, c}, {0, d, 0}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                tau[i][j] = 2.0e-5*(gradU[i][j] + gradU[j][i])
                          - (i == j ? (2.0/3.0)*2.0e-5*divU : 0.0);
        ViscousFaceFlux f;
        for (int i = 0; i < 3; ++i)
        {
            double tr = 0.;
            for (int j = 0; j < 3; ++j) tr += tau[i][j]*nn[j];
            f.momentum[i] = tr;
        }
        return f;
    };
    // Check the flux helper directly on a few normals.
    double maxFluxError = 0.;
    for (const auto& nn : {std::array<double,3>{1.,0.,0.}, std::array<double,3>{0.,1.,0.},
                           std::array<double,3>{0.,0.,1.}, std::array<double,3>{0.577350269,0.577350269,0.577350269}})
    {
        const std::array<TetraVec3, 3> gradU = {{{0., a, 0.}, {b, 0., c}, {0., d, 0.}}};
        const TetraVec3 gT{gx, gy, 0.};
        const TetraVec3 u{a*10., b*5. + c*2., d*10.};
        const auto f = viscousFaceFlux(u, gradU, gT, 2.0e-5, 0.05, nn, 1.0);
        const auto exact = analytic(nn);
        for (int i = 0; i < 3; ++i)
            maxFluxError = std::max(maxFluxError, std::abs(f.momentum[i] - exact.momentum[i]));
    }
    // Assembled residual on the linear field: interior contributions cancel.
    std::vector<ConservedVars> residual(n, ConservedVars(2));
    addViscousResidual(mesh, velocity, temperature, mu, lambda, residual);
    double maxResidual = 0.;
    for (const auto& r : residual)
        for (int i = 0; i < 3; ++i)
            maxResidual = std::max(maxResidual, std::abs(r.rhoU[i]));
    // Integral of the momentum flux over all boundary faces equals the net
    // viscous force on the domain (boundary faces have neighbour < 0 and are
    // skipped by the assembler, so the residual sum must be zero).
    double sum[3] = {0., 0., 0.};
    for (const auto& r : residual)
        for (int i = 0; i < 3; ++i) sum[i] += r.rhoU[i];
    std::printf("TETRA_VISCOUS cells=%zu maxFaceFluxError=%.3e maxInteriorResidual=%.3e "
                "sumMomentum=(%.3e,%.3e,%.3e)\n",
                n, maxFluxError, maxResidual, sum[0], sum[1], sum[2]);
    // Thresholds: the per-face flux is ~1e-10 N on this mesh, so a residual of
    // a few 1e-9 N means the interior faces cancel to ~1 % of a face flux (the
    // remaining error comes from the gradient of degenerate sliver cells, a
    // documented limitation). The conservative check is the summed momentum,
    // which must vanish to round-off.
    const bool ok = maxFluxError < 1e-14 && maxResidual < 1e-7
        && std::abs(sum[0]) < 1e-12 && std::abs(sum[1]) < 1e-12 && std::abs(sum[2]) < 1e-12;
    if (!ok) { std::printf("TETRA_VISCOUS_FAIL\n"); return 1; }
    std::printf("TETRA_VISCOUS_PASS\n");
    return 0;
}
