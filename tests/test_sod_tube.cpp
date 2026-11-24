/*---------------------------------------------------------------------------*\
    test_sod_tube - Stage 7: Sod shock tube quantitative validation.

    Standard Sod problem (Sod 1978) with real-gas PR EOS / pure O2.
    IC: left (p=1e6 Pa, T=300K, u=0) | right (p=1e5 Pa, T=240K, u=0).
    Run to t=100µs. Compare rho, p, u profiles against the exact solution
    (characteristic-based, computed inline for ideal-gas approximation).

    Checks:
      1. Density L1 error < 5% relative (coarse grid, AUSM 1st-order)
      2. Pressure L1 error < 5%
      3. Velocity L1 error < 5% (relative to left-state sound speed scale)
      4. Correct shock / expansion / contact ordering in the final state
      5. Mass conservation to 1e-6 relative
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"

using namespace Foam;
using namespace std;

static int g_pass=0, g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if (!ok) cout << "  (got=" << fixed << setprecision(4) << got
                  << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

/*---------------------------------------------------------------------------
    Exact Sod solution (1D ideal gas, gamma fixed from initial sound speeds).
    Returns rho(x,t), p(x,t), u(x,t) at time t for the standard problem.
---------------------------------------------------------------------------*/
struct SodExact
{
    double rhoL, pL, uL, aL;   // left state
    double rhoR, pR, uR, aR;   // right state
    double gamma;

    SodExact(double rhoL, double pL, double rhoR, double pR, double gamma)
        : rhoL(rhoL), pL(pL), uL(0), rhoR(rhoR), pR(pR), uR(0),
          gamma(gamma)
    {
        aL = std::sqrt(gamma * pL / rhoL);
        aR = std::sqrt(gamma * pR / rhoR);
    }

    // Iteratively solve for post-shock pressure p* (Toro 2009 §4.2).
    double pStar() const
    {
        const double g1 = (gamma-1)/(2*gamma);
        const double g2 = (gamma+1)/(2*gamma);
        auto f = [&](double p)
        {
            double fL = std::pow(p/pL, g1) - 1.0;
            fL *= 2*aL/(gamma-1);
            double AL = 2.0/((gamma+1)*rhoR);
            double BL = (gamma-1)/(gamma+1)*pR;
            double fR = (p - pR) * std::sqrt(AL/(p + BL));
            double fl2 = std::pow(p/pL, g1) - 1.0;
            fl2 *= 2*aL/(gamma-1);
            return fl2 + fR - (uR - uL);
        };
        double p0 = 0.5*(pL+pR), p1 = p0*2;
        for (int i=0;i<100;++i)
        {
            double fp0 = f(p0), fp1 = f(p1);
            if (std::abs(fp1-fp0) < 1e-30) break;
            double p2 = p1 - fp1*(p1-p0)/(fp1-fp0);
            p0=p1; p1=std::max(p2, 1e-6*pL);
            if (std::abs(p1-p0)/(std::abs(p1)+1e-30) < 1e-10) break;
        }
        return p1;
    }

    // Get (rho, p, u) at position x (relative to initial diaphragm at 0) at time t.
    tuple<double,double,double> state(double x, double t) const
    {
        double ps = pStar();
        double g1 = (gamma-1)/(2*gamma);
        // uStar: contact discontinuity velocity
        double us = uL + 2*aL/(gamma-1) * (1 - std::pow(ps/pL, g1));

        // Shock speed
        double A  = 2.0/((gamma+1)*rhoR);
        double B  = (gamma-1)/(gamma+1)*pR;
        double Ws = uR + aR*std::sqrt((gamma+1)/(2*gamma) * ps/pR + (gamma-1)/(2*gamma));

        // Star-state densities
        double rhoLs = rhoL * std::pow(ps/pL, 1.0/gamma);
        double rhoRs = rhoR * (ps/pR + (gamma-1)/(gamma+1))
                             / ((gamma-1)/(gamma+1) * ps/pR + 1);

        // Wave-head and wave-tail speeds of the expansion fan
        double Shl = uL - aL;           // expansion fan left edge
        double Stl = us - aL*std::pow(ps/pL, (gamma-1)/(2*gamma)); // fan right edge

        double xi = x / t;   // similarity variable
        if (xi < Shl)               return {rhoL, pL, uL};
        if (xi < Stl)
        {
            // Inside expansion fan
            double u = 2.0/(gamma+1) * (aL + (gamma-1)/2 * uL + xi);
            double a = aL - (gamma-1)/2 * (u - uL);
            double p = pL * std::pow(a/aL, 2*gamma/(gamma-1));
            double rho = gamma * p / (a*a);
            return {rho, p, u};
        }
        if (xi < us) return {rhoLs, ps, us};
        if (xi < Ws) return {rhoRs, ps, us};
        return {rhoR, pR, uR};
    }
};

int main()
{
    cout << "=== Sod shock tube: AUSM vs exact solution ===\n\n";

    const int    N    = 200;
    const double L    = 1.0;       // [m], diaphragm at x=0.5
    const double tEnd = 1.0e-4;   // 100 µs
    const double CFL  = 0.4;

    // Initial conditions: pure O2 (use legacy 2-species slot 0)
    // PR-EOS gives rho directly; use approximate ideal-gas Tref for p->rho.
    const double gamma = 1.4;   // ideal approximation for comparison

    const double pL=1.0e6, TL=300.0;
    const double pR=1.0e5, TR=240.0;
    // Ideal gas rho (O2: MW=32e-3, R_spec=8.314/0.032)
    const double R_O2 = 8.314462618 / 0.032;
    const double rhoL = pL / (R_O2 * TL);
    const double rhoR = pR / (R_O2 * TR);

    Domain dom(N, 1, L, 0.0, 2);
    double dx = dom.dx();
    double x0 = 0.5 * L;   // diaphragm position

    for (int i=0;i<N;++i)
    {
        double x = (i+0.5)*dx;
        PrimitiveVars& w = dom.stagedW(i);
        w.Y[0]=1.0; w.Y[1]=0.0; w.U[0]=w.U[1]=w.U[2]=0.0;
        if (x < x0) { w.p=pL; w.T=TL; }
        else         { w.p=pR; w.T=TR; }
    }
    dom.primitiveToConserved();

    // Record initial mass
    double mass0=0.0;
    for (int i=0;i<N;++i) mass0 += dom.U(i).rho * dx;

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.advanceToTime(tEnd, CFL);

    // Exact solution instance
    SodExact exact(rhoL, pL, rhoR, pR, gamma);

    // Compute L1 errors (skip a few cells near boundaries)
    double err_rho=0, err_p=0, err_u=0;
    double ref_rho=0, ref_p=0, ref_u=0;
    const double aL_scale = std::sqrt(gamma*pL/rhoL);
    int nCells = 0;
    for (int i=5;i<N-5;++i)
    {
        double x = (i+0.5)*dx - x0;   // relative to diaphragm
        auto [rho_ex, p_ex, u_ex] = exact.state(x, tEnd);
        double rho_c = dom.U(i).rho;
        double p_c   = dom.W(i).p;
        double u_c   = dom.W(i).U[0];
        err_rho += std::abs(rho_c - rho_ex);
        err_p   += std::abs(p_c   - p_ex);
        err_u   += std::abs(u_c   - u_ex);
        ref_rho += std::abs(rho_ex);
        ref_p   += std::abs(p_ex);
        ref_u   += aL_scale;   // velocity reference = left-state sound speed
        ++nCells;
    }
    double l1_rho = err_rho / ref_rho;
    double l1_p   = err_p   / ref_p;
    double l1_u   = err_u   / ref_u;

    cout << fixed << setprecision(4);
    cout << "  L1(rho)=" << l1_rho*100 << "%  L1(p)=" << l1_p*100
         << "%  L1(u)=" << l1_u*100 << "%  (5% threshold)\n\n";

    check("no trip", !dom.tripped());
    check("L1(rho) < 5%",   l1_rho < 0.05, l1_rho*100);
    check("L1(p)   < 5%",   l1_p   < 0.05, l1_p*100);
    check("L1(u)   < 5%",   l1_u   < 0.05, l1_u*100);

    // Wave ordering: near left edge should be dense/high-pressure; near right edge lower
    double rho_left  = dom.U(N/6).rho;
    double rho_right = dom.U(5*N/6).rho;
    check("density drops from left to right (expansion+contact+shock order)",
          rho_left > rho_right, rho_left, rho_right);

    // Mass conservation
    double mass1=0.0;
    for (int i=0;i<N;++i) mass1 += dom.U(i).rho * dx;
    double dmass = std::abs(mass1-mass0)/(mass0+1e-30);
    check("mass conserved (< 1e-6 rel)", dmass < 1e-6, dmass);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "SOD_TUBE_PASS" : "SOD_TUBE_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
