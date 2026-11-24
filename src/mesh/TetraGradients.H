#ifndef RFAM_TETRA_GRADIENTS_H
#define RFAM_TETRA_GRADIENTS_H
// Cell-centred weighted least-squares gradients and face gradients on the
// tetra mesh. This is the first piece of the viscous/turbulence infrastructure
// (Phase 2): viscous fluxes, SST production and heat conduction all consume
// these gradients. No thermodynamics and no model assumptions live here.
#include "TetraMesh.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{
using TetraVec3 = std::array<double, 3>;

class TetraGradients
{
    const TetraMesh& mesh_;
    double regularization_;
public:
    // regularization: relative Tikhonov term for degenerate stencils.
    explicit TetraGradients(const TetraMesh& mesh, double regularization = 1e-12)
      : mesh_(mesh), regularization_(regularization) {}

    // Weighted least-squares cell gradient of a scalar cell field.
    // boundaryFaceValue (optional): values at boundary faces (wall/BC values).
    // Boundary cells whose interior neighbours are degenerate (e.g. flattened
    // plate-adjacent tetrahedra) need these to span the normal direction; when
    // absent such directions stay unconstrained and the regularization keeps
    // the corresponding gradient component at zero (the zero-gradient BC).
    std::vector<TetraVec3> cellGradient(const std::vector<double>& phi,
                                        const std::vector<double>* weight = nullptr,
                                        const std::vector<double>* boundaryFaceValue = nullptr) const
    {
        if (phi.size() != mesh_.cellCount())
            throw std::invalid_argument("gradient field size mismatch");
        if (boundaryFaceValue && boundaryFaceValue->size() != mesh_.faceCount())
            throw std::invalid_argument("gradient boundary field size mismatch");
        std::vector<TetraVec3> gradient(mesh_.cellCount(), TetraVec3{0., 0., 0.});
        for (std::size_t c = 0; c < mesh_.cellCount(); ++c)
        {
            const auto xc = mesh_.cellCentroid(c);
            double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            double b[3] = {0, 0, 0};
            for (auto f : mesh_.cellFaces(c))
            {
                const auto& face = mesh_.face(f);
                const int other = face.owner == static_cast<int>(c) ? face.neighbour : face.owner;
                if (other < 0 && !boundaryFaceValue) continue;
                const auto xn = other < 0 ? mesh_.faceCentroid(f)
                                          : mesh_.cellCentroid(static_cast<std::size_t>(other));
                const double d[3] = {xn.x - xc.x, xn.y - xc.y, xn.z - xc.z};
                const double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                if (!(d2 > 0.)) continue;
                const double w = (weight && other >= 0 ? (*weight)[static_cast<std::size_t>(other)] : 1.0) / d2;
                const double dphi = other < 0 ? (*boundaryFaceValue)[f] - phi[c]
                                              : phi[static_cast<std::size_t>(other)] - phi[c];
                for (int i = 0; i < 3; ++i)
                {
                    b[i] += w * d[i] * dphi;
                    for (int j = 0; j < 3; ++j) A[i][j] += w * d[i] * d[j];
                }
            }
            double trace = 0.;
            for (int i = 0; i < 3; ++i) trace += A[i][i];
            // Adaptive Tikhonov term: well-conditioned stencils stay exact for
            // linear fields (relative perturbation 1e-12); only degenerate
            // stencils (e.g. flattened plate-adjacent tetrahedra whose interior
            // neighbours are coplanar) are regularized, where the normal
            // direction is otherwise unconstrained and the zero-gradient BC
            // keeps that component at zero.
            const double detA =
                A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
              - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
              + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            const double scale = std::cbrt(std::max(trace, 0.)/3.);
            const double eps = std::abs(detA) > 1e-14*std::pow(scale, 3.)
                ? regularization_*std::max(trace, 1e-300)
                : 1e-9*std::max(trace, 1e-300);
            for (int i = 0; i < 3; ++i) A[i][i] += eps;
            solve3(A, b, gradient[c]);
        }
        return gradient;
    }

    // Face value interpolation plus non-orthogonal correction:
    //   phi_f = w_P phi_P + w_N phi_N + (g_avg . (x_f - x_w)) correction
    // The plain linear interpolation is returned when the correction weights
    // are degenerate.
    std::vector<double> faceValues(const std::vector<double>& phi,
                                   const std::vector<TetraVec3>& cellGradient) const
    {
        if (phi.size() != mesh_.cellCount() || cellGradient.size() != mesh_.cellCount())
            throw std::invalid_argument("faceValues field size mismatch");
        std::vector<double> values(mesh_.faceCount(), 0.);
        for (std::size_t f = 0; f < mesh_.faceCount(); ++f)
        {
            const auto& face = mesh_.face(f);
            const std::size_t P = static_cast<std::size_t>(face.owner);
            if (face.neighbour < 0) { values[f] = phi[P]; continue; }
            const std::size_t N = static_cast<std::size_t>(face.neighbour);
            const auto xP = mesh_.cellCentroid(P), xN = mesh_.cellCentroid(N);
            const auto xf = mesh_.faceCentroid(f);
            const double dx[3] = {xN.x - xP.x, xN.y - xP.y, xN.z - xP.z};
            const double df[3] = {xf.x - xP.x, xf.y - xP.y, xf.z - xP.z};
            const double d2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            const double wN = d2 > 0. ? (df[0]*dx[0] + df[1]*dx[1] + df[2]*dx[2]) / d2 : 0.5;
            const double wP = 1. - wN;
            const double avg[3] = {
                0.5*(cellGradient[P][0] + cellGradient[N][0]),
                0.5*(cellGradient[P][1] + cellGradient[N][1]),
                0.5*(cellGradient[P][2] + cellGradient[N][2])};
            // Non-orthogonal correction: the component of the gradient along
            // the vector from the interpolation point to the face centre.
            const double t[3] = {df[0] - wN*dx[0], df[1] - wN*dx[1], df[2] - wN*dx[2]};
            values[f] = wP*phi[P] + wN*phi[N] + avg[0]*t[0] + avg[1]*t[1] + avg[2]*t[2];
        }
        return values;
    }

    // Face-normal gradient for a scalar zero-gradient-free field:
    //   (grad phi . n)_f = (phi_N - phi_P)/|d| * (n . dhat) + correction
    std::vector<double> faceNormalGradient(const std::vector<double>& phi,
                                           const std::vector<TetraVec3>& cellGradient) const
    {
        if (phi.size() != mesh_.cellCount() || cellGradient.size() != mesh_.cellCount())
            throw std::invalid_argument("faceNormalGradient field size mismatch");
        std::vector<double> values(mesh_.faceCount(), 0.);
        for (std::size_t f = 0; f < mesh_.faceCount(); ++f)
        {
            const auto& face = mesh_.face(f);
            if (face.neighbour < 0) continue;   // wall/boundary value supplied by the caller
            const std::size_t P = static_cast<std::size_t>(face.owner);
            const std::size_t N = static_cast<std::size_t>(face.neighbour);
            const auto xP = mesh_.cellCentroid(P), xN = mesh_.cellCentroid(N);
            const double dx[3] = {xN.x - xP.x, xN.y - xP.y, xN.z - xP.z};
            const double d = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
            if (!(d > 0.)) continue;
            const auto area = mesh_.faceAreaVector(f);
            const double a = std::sqrt(area.x*area.x + area.y*area.y + area.z*area.z);
            if (!(a > 0.)) continue;
            const double n[3] = {area.x/a, area.y/a, area.z/a};
            const double sn = (dx[0]*n[0] + dx[1]*n[1] + dx[2]*n[2]);
            const double dn = std::abs(sn) > 1e-12 ? sn : 1e-12;
            const double orthogonal = (phi[N] - phi[P]) / dn;
            const double avg[3] = {
                0.5*(cellGradient[P][0] + cellGradient[N][0]),
                0.5*(cellGradient[P][1] + cellGradient[N][1]),
                0.5*(cellGradient[P][2] + cellGradient[N][2])};
            // Remove the component of the averaged gradient along the face
            // direction from the raw face difference (De Villiers & O'Brien
            // style correction for skewed tetrahedral stencils).
            const double correction = avg[0]*n[0] + avg[1]*n[1] + avg[2]*n[2];
            const double along = (avg[0]*dx[0] + avg[1]*dx[1] + avg[2]*dx[2]) / dn;
            values[f] = orthogonal + (correction - along)*1.0;
        }
        return values;
    }

private:
    static void solve3(double A[3][3], const double b[3], TetraVec3& x)
    {
        // Cramer's rule with a partial-pivot fallback; the system is symmetric
        // positive definite for non-degenerate tetrahedral stencils.
        const double detA =
            A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
          - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
          + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
        if (!(std::abs(detA) > 0.))
        {
            x = {0., 0., 0.};
            return;
        }
        const auto det = [](const double a[3][3])
        {
            return a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
                 - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
                 + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
        };
        for (int column = 0; column < 3; ++column)
        {
            double M[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    M[i][j] = j == column ? b[i] : A[i][j];
            x[column] = det(M) / detA;
        }
    }
};
}
#endif
