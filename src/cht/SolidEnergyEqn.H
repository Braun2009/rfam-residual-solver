#ifndef SolidEnergyEqn_H
#define SolidEnergyEqn_H

/*---------------------------------------------------------------------------*\
    SolidEnergyEqn - 1D solid heat conduction (Stage 4 CHT, standalone)

        rho*cp dT/dt = d/dx( k dT/dx )

    Cell-centered finite volume, implicit Backward-Euler in time (or steady),
    solved with the Thomas tridiagonal algorithm (unconditionally stable ->
    CHT-friendly: no dt < dx^2/2alpha restriction).

    Boundary conditions per face: FixedT / FixedFlux / Convection / Adiabatic.
    Convection uses a series resistance (film + half-cell conduction):
        g = 1 / (1/h + (dx/2)/k),   q_into_solid = g*(T_inf - T_boundaryCell)

    Pure std C++ (no OpenFOAM/Eigen), matches RFAM module style.
\*---------------------------------------------------------------------------*/

// ---------------------------------------------------------------------------
// Notes for use in the thrust-chamber CHT chain.
//
// 1. The implicit Backward-Euler form is what makes this usable at all inside
//    a CHT loop: the solid's thermal diffusivity is large compared to the
//    fluid's effective one, so an explicit solid solve would set the global
//    timestep. Do not "simplify" this to explicit Euler to share the fluid
//    step.
//
// 2. rho, cp and k are single scalars for the whole domain, and the capacity
//    term is rho*cp*dx/dt with dx uniform. For a multi-material wall (e.g. a
//    copper liner against a steel jacket) that is wrong at the material
//    boundary. Build one SolidEnergyEqn per material and couple them through
//    CHTCoupler rather than fitting an averaged rho*cp.
//
// 3. FixedFlux enters the solid with a POSITIVE sign: q is defined as heat
//    flowing INTO the solid. Reversing it is the classic CHT sign error and
//    produces a wall that heats the fluid when the fluid is heating it. The
//    convention is stated here and every caller that derives q from a
//    fluid-side closure must preserve it.
//
// 4. T_ is committed state. Callers must not write T_ directly (bypassing the
//    solve): that breaks the rollback contract CHTCoupler and the fluid step
//    rely on. Use the accessor and let the solve own the update.
// ---------------------------------------------------------------------------

#include <vector>
#include <algorithm>

namespace Foam
{

enum class SolidBCType { FixedT, FixedFlux, Convection, Adiabatic };

struct SolidBC
{
    SolidBCType type = SolidBCType::Adiabatic;
    double a = 0.0;   // FixedT: T_wall | FixedFlux: q'' (into solid +) | Convection: h
    double b = 0.0;   // Convection: T_inf

    static SolidBC fixedT(double T)              { return {SolidBCType::FixedT, T, 0.0}; }
    static SolidBC fixedFlux(double q)           { return {SolidBCType::FixedFlux, q, 0.0}; }
    static SolidBC convection(double h, double T){ return {SolidBCType::Convection, h, T}; }
    static SolidBC adiabatic()                   { return {SolidBCType::Adiabatic, 0.0, 0.0}; }
};

class SolidEnergyEqn
{
private:
    int    N_;
    double L_, dx_;
    double rho_, cp_, k_;
    std::vector<double> T_;

    // Thomas solver for tridiagonal system (a:sub, b:diag, c:sup, d:rhs) -> x
    static void thomas(std::vector<double> a, std::vector<double> b,
                       std::vector<double> c, std::vector<double> d,
                       std::vector<double>& x)
    {
        int n = (int)b.size();
        x.assign(n, 0.0);
        for (int i = 1; i < n; ++i)
        {
            double m = a[i] / b[i-1];
            b[i] -= m * c[i-1];
            d[i] -= m * d[i-1];
        }
        x[n-1] = d[n-1] / b[n-1];
        for (int i = n-2; i >= 0; --i) x[i] = (d[i] - c[i]*x[i+1]) / b[i];
    }

    // Boundary face conductance g [W/m2/K] and driving temperature contribution.
    // Returns g; sets Tdrive such that flux_into_solid = g*(Tdrive - T_cell).
    // For FixedFlux, g=0 and the flux is added directly via *fluxOut.
    double faceCond(const SolidBC& bc, double& Tdrive, double& fluxOut) const
    {
        Tdrive = 0.0; fluxOut = 0.0;
        switch (bc.type)
        {
            case SolidBCType::FixedT:     Tdrive = bc.a;  return k_/(0.5*dx_);
            case SolidBCType::Convection: Tdrive = bc.b;  return 1.0/(1.0/bc.a + 0.5*dx_/k_);
            case SolidBCType::FixedFlux:  fluxOut = bc.a; return 0.0;   // q'' into solid
            case SolidBCType::Adiabatic:  default:        return 0.0;
        }
    }

    // Assemble & solve. dt>0 => Backward-Euler transient; dt<=0 => steady.
    void assembleSolve(const SolidBC& left, const SolidBC& right, double dt)
    {
        std::vector<double> a(N_,0.0), b(N_,0.0), c(N_,0.0), d(N_,0.0), x;
        const double kA  = k_ / dx_;                                   // internal face conductance
        const double cap = (dt > 0.0) ? rho_*cp_*dx_/dt : 0.0;         // capacitance / dt

        double gL, TdL, qL, gR, TdR, qR;
        gL = faceCond(left,  TdL, qL);
        gR = faceCond(right, TdR, qR);

        for (int i = 0; i < N_; ++i)
        {
            double diag = cap;
            double rhs  = (dt > 0.0) ? cap*T_[i] : 0.0;

            // west face
            if (i > 0) { a[i] = -kA; diag += kA; }
            else       { diag += gL; rhs += gL*TdL + qL; }   // left boundary

            // east face
            if (i < N_-1) { c[i] = -kA; diag += kA; }
            else          { diag += gR; rhs += gR*TdR + qR; } // right boundary

            b[i] = diag; d[i] = rhs;
        }
        thomas(a, b, c, d, x);
        T_ = x;
    }

public:
    SolidEnergyEqn(int N, double L, double rho, double cp, double k, double T0 = 300.0)
        : N_(N), L_(L), dx_(L/N), rho_(rho), cp_(cp), k_(k), T_(N, T0) {}

    void setUniform(double T) { std::fill(T_.begin(), T_.end(), T); }
    void setField(const std::vector<double>& T) { if ((int)T.size()==N_) T_=T; }

    void step(double dt, const SolidBC& left, const SolidBC& right) { assembleSolve(left, right, dt); }
    void solveSteady(const SolidBC& left, const SolidBC& right)     { assembleSolve(left, right, -1.0); }

    int    N()        const { return N_; }
    double x(int i)   const { return (i + 0.5) * dx_; }
    double T(int i)   const { return T_[i]; }
    double dx()       const { return dx_; }
    double k()        const { return k_; }
    const std::vector<double>& field() const { return T_; }

    // Heat flux into the solid at the boundary face [W/m2]
    double leftFlux(const SolidBC& left) const
    {
        double Td, q; double g = faceCond(left, Td, q);
        return (left.type == SolidBCType::FixedFlux) ? q : g*(Td - T_[0]);
    }
    double rightFlux(const SolidBC& right) const
    {
        double Td, q; double g = faceCond(right, Td, q);
        return (right.type == SolidBCType::FixedFlux) ? q : g*(Td - T_[N_-1]);
    }
};

} // End namespace Foam

#endif
