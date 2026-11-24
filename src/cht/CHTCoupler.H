#ifndef CHTCoupler_H
#define CHTCoupler_H

/*---------------------------------------------------------------------------*\
    CHTCoupler - monolithic 1D conjugate heat transfer (Stage 4 CHT)

    Strong (monolithic) coupling of two conduction domains A and B sharing an
    interface (e.g. fluid-side near-wall layer / solid wall, or two solid
    materials). A single global tridiagonal system is assembled across BOTH
    domains, so at the A|B interface:
        - temperature continuity, and
        - heat-flux balance
    are satisfied implicitly in one solve (no partitioned sub-iteration).

    Interface face conductance = series (harmonic) of the two half-cells:
        g_int = 1 / ( 0.5*dxA/kA + 0.5*dxB/kB )

    Implicit Backward-Euler (transient) or steady. Thomas tridiagonal solve.
    Pure std C++; reuses SolidBC from SolidEnergyEqn.H.
\*---------------------------------------------------------------------------*/

// ---------------------------------------------------------------------------
// Integration notes (read before wiring this into a solver).
//
// 1. SCOPE. This couples CONDUCTION domains only. The "fluid side" of a CHT
//    problem must first be reduced to a wall heat flux or a fixed wall
//    temperature before it reaches this class. It does not solve the fluid,
//    and it does not know about convection inside the fluid.
//
// 2. Monolithic here means one linear system, not one nonlinear residual. The
//    coupling is therefore exact only while k, rho*cp and the boundary data
//    are held fixed across the solve. A strongly nonlinear wall (cryogenic
//    LOX-side k varying with T, or a coolant that changes phase) violates that
//    assumption, and the correct treatment is an outer Newton iteration that
//    re-assembles this system per iteration -- not a sub-iteration on the
//    interface flux, which would silently return to the partitioned scheme
//    this class exists to avoid.
//
// 3. The 1D assumption is a real restriction. A rocket chamber wall is a
//    solid with a coolant channel; representing it as a chain of 1D cells is
//    defensible along a wall normal through a thin plate, and NOT defensible
//    where the heat flow is two-dimensional (near a coolant-channel corner, or
//    at the throat where axial conduction matters). Document which one each
//    case is before quoting a wall temperature.
//
// 4. Conservation is the property to test, not the code path: total energy
//    leaving A must enter B to solver tolerance, and the interface residual
//    must shrink with a refinement of dx and dt. The existing
//    test_cht2d_conservation and test_cht_relaxation_vs_conservative targets
//    are the right shape; keep them as the gate for any change here.
//
// 5. The steady solve (dt <= 0) drops the capacity term entirely. Mixing a
//    steady solid solve with a transient fluid march is a modelling decision:
//    the solid is then assumed to have no thermal inertia, which is wrong for
//    a startup transient and right for a burn-time-averaged wall temperature.
// ---------------------------------------------------------------------------

#include "SolidEnergyEqn.H"   // SolidBC, SolidBCType
#include <vector>

namespace Foam
{

class CHTCoupler
{
private:
    int    NA_, NB_, M_;
    double dxA_, dxB_, kA_, kB_, rcA_, rcB_;   // rc = rho*cp
    std::vector<double> T_;

    static void thomas(std::vector<double> a, std::vector<double> b,
                       std::vector<double> c, std::vector<double> d,
                       std::vector<double>& x)
    {
        int n = (int)b.size(); x.assign(n, 0.0);
        for (int i = 1; i < n; ++i) { double m=a[i]/b[i-1]; b[i]-=m*c[i-1]; d[i]-=m*d[i-1]; }
        x[n-1] = d[n-1]/b[n-1];
        for (int i = n-2; i >= 0; --i) x[i] = (d[i]-c[i]*x[i+1])/b[i];
    }

    double kOf(int i)  const { return (i < NA_) ? kA_  : kB_;  }
    double dxOf(int i) const { return (i < NA_) ? dxA_ : dxB_; }
    double rcOf(int i) const { return (i < NA_) ? rcA_ : rcB_; }

    static double faceBC(const SolidBC& bc, double& Td, double& q, double dx, double k)
    {
        Td = 0.0; q = 0.0;
        switch (bc.type)
        {
            case SolidBCType::FixedT:     Td = bc.a; return k/(0.5*dx);
            case SolidBCType::Convection: Td = bc.b; return 1.0/(1.0/bc.a + 0.5*dx/k);
            case SolidBCType::FixedFlux:  q  = bc.a; return 0.0;
            default:                                 return 0.0;
        }
    }

    void assembleSolve(const SolidBC& left, const SolidBC& right, double dt)
    {
        std::vector<double> a(M_,0.0), b(M_,0.0), c(M_,0.0), d(M_,0.0), x;
        for (int i = 0; i < M_; ++i)
        {
            double dx  = dxOf(i);
            double cap = (dt > 0.0) ? rcOf(i)*dx/dt : 0.0;
            double diag = cap, rhs = (dt > 0.0) ? cap*T_[i] : 0.0;

            if (i > 0) {                       // west internal face (handles A|B interface too)
                double g = 1.0/(0.5*dxOf(i-1)/kOf(i-1) + 0.5*dx/kOf(i));
                a[i] = -g; diag += g;
            } else {                           // left boundary (domain A)
                double g,Td,q; g = faceBC(left, Td, q, dx, kOf(i)); diag += g; rhs += g*Td + q;
            }

            if (i < M_-1) {                    // east internal face
                double g = 1.0/(0.5*dx/kOf(i) + 0.5*dxOf(i+1)/kOf(i+1));
                c[i] = -g; diag += g;
            } else {                           // right boundary (domain B)
                double g,Td,q; g = faceBC(right, Td, q, dx, kOf(i)); diag += g; rhs += g*Td + q;
            }
            b[i] = diag; d[i] = rhs;
        }
        thomas(a, b, c, d, x);
        T_ = x;
    }

public:
    CHTCoupler(int NA, double LA, double kA, double rcA,
               int NB, double LB, double kB, double rcB, double T0 = 300.0)
        : NA_(NA), NB_(NB), M_(NA+NB),
          dxA_(LA/NA), dxB_(LB/NB), kA_(kA), kB_(kB), rcA_(rcA), rcB_(rcB),
          T_(NA+NB, T0) {}

    void step(double dt, const SolidBC& l, const SolidBC& r) { assembleSolve(l, r, dt); }
    void solveSteady(const SolidBC& l, const SolidBC& r)     { assembleSolve(l, r, -1.0); }

    int    M()  const { return M_; }
    int    NA() const { return NA_; }
    double LA() const { return NA_*dxA_; }
    double x(int i) const { return (i < NA_) ? (i+0.5)*dxA_ : NA_*dxA_ + (i-NA_+0.5)*dxB_; }
    double T(int i) const { return T_[i]; }

    // Heat flux across the A|B interface [W/m2] (A -> B positive)
    double interfaceFlux() const
    {
        int i = NA_-1;
        double g = 1.0/(0.5*dxA_/kA_ + 0.5*dxB_/kB_);
        return g*(T_[i] - T_[i+1]);
    }
};

} // End namespace Foam

#endif
