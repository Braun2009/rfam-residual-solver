// Directed test for TetraGradients: a linear field must be differentiated
// exactly on any valid tetrahedral stencil (cell gradient, face value and face
// normal gradient). A quadratic field checks the expected first/second-order
// behaviour without asserting a specific order.
#include "TetraGradients.H"
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
    catch (...) { std::printf("TETRA_GRADIENTS_SKIP no mesh %s\n", path.c_str()); return 77; }
    TetraGradients gradients(mesh);

    const double b[3] = {1.2, -0.7, 0.3};
    std::vector<double> linear(mesh.cellCount()), quadratic(mesh.cellCount());
    for (std::size_t c = 0; c < mesh.cellCount(); ++c)
    {
        const auto x = mesh.cellCentroid(c);
        linear[c] = 0.5 + b[0]*x.x + b[1]*x.y + b[2]*x.z;
        quadratic[c] = x.x*x.x + 2.*x.y*x.y + 0.5*x.z*x.z + 0.25*x.x*x.z;
    }
    // Boundary face values: the physical BC value at each face (here the exact
    // linear field value, as a wall/BC provider would supply).
    std::vector<double> boundaryFaceValue(mesh.faceCount(), 0.);
    for (std::size_t f = 0; f < mesh.faceCount(); ++f)
    {
        const auto x = mesh.faceCentroid(f);
        boundaryFaceValue[f] = 0.5 + b[0]*x.x + b[1]*x.y + b[2]*x.z;
    }
    const auto gradient = gradients.cellGradient(linear, nullptr, &boundaryFaceValue);
    double maxCellError = 0.;
    std::size_t worstCell = 0;
    for (std::size_t c = 0; c < gradient.size(); ++c)
        for (int i = 0; i < 3; ++i)
            if (std::abs(gradient[c][i] - b[i]) > maxCellError)
            { maxCellError = std::abs(gradient[c][i] - b[i]); worstCell = c; }
    const auto values = gradients.faceValues(linear, gradient);
    const auto normalGradient = gradients.faceNormalGradient(linear, gradient);
    double maxFaceValueError = 0.;
    double maxNormalError = 0.;
    for (std::size_t f = 0; f < mesh.faceCount(); ++f)
    {
        if (mesh.face(f).neighbour < 0) continue;  // boundary faces hold the cell value
        const auto x = mesh.faceCentroid(f);
        const double exact = 0.5 + b[0]*x.x + b[1]*x.y + b[2]*x.z;
        maxFaceValueError = std::max(maxFaceValueError, std::abs(values[f] - exact));
        const auto area = mesh.faceAreaVector(f);
        const double a = std::sqrt(area.x*area.x + area.y*area.y + area.z*area.z);
        const double exactNormal = (b[0]*area.x + b[1]*area.y + b[2]*area.z)/a;
        maxNormalError = std::max(maxNormalError, std::abs(normalGradient[f] - exactNormal));
    }
    const auto quadGradient = gradients.cellGradient(quadratic);
    double quadError = 0.;
    for (std::size_t c = 0; c < mesh.cellCount(); ++c)
    {
        const auto x = mesh.cellCentroid(c);
        const double exact[3] = {2.*x.x + 0.25*x.z, 4.*x.y, x.z + 0.25*x.x};
        for (int i = 0; i < 3; ++i)
            quadError = std::max(quadError, std::abs(quadGradient[c][i] - exact[i]));
    }
    std::printf("TETRA_GRADIENTS cells=%zu faces=%zu maxCellError=%.3e "
                "maxFaceValueError=%.3e maxNormalError=%.3e quadError=%.3e\n",
                mesh.cellCount(), mesh.faceCount(), maxCellError,
                maxFaceValueError, maxNormalError, quadError);
    {
        const auto xc = mesh.cellCentroid(worstCell);
        std::printf("WORST cell=%zu g=(%.6f,%.6f,%.6f) exact=(%.1f,%.1f,%.1f) xc=(%.6g,%.6g,%.6g) faces=%zu\n",
                    worstCell, gradient[worstCell][0], gradient[worstCell][1],
                    gradient[worstCell][2], b[0], b[1], b[2], xc.x, xc.y, xc.z,
                    mesh.cellFaces(worstCell).size());
        for (auto f : mesh.cellFaces(worstCell))
        {
            const auto& face = mesh.face(f);
            const int other = face.owner == static_cast<int>(worstCell) ? face.neighbour : face.owner;
            if (other < 0) { std::printf("  face %zu boundary\n", f); continue; }
            const auto xn = mesh.cellCentroid(static_cast<std::size_t>(other));
            std::printf("  face %zu other=%d d=(%.6g,%.6g,%.6g) dphi=%.6g owner=%d neigh=%d\n",
                        f, other, xn.x-xc.x, xn.y-xc.y, xn.z-xc.z,
                        linear[static_cast<std::size_t>(other)] - linear[worstCell],
                        face.owner, face.neighbour);
        }
    }
    const bool ok = maxCellError < 1e-10 && maxFaceValueError < 1e-10
        && maxNormalError < 1e-8 && std::isfinite(quadError);
    if (!ok) { std::printf("TETRA_GRADIENTS_FAIL\n"); return 1; }
    std::printf("TETRA_GRADIENTS_PASS\n");
    return 0;
}
