#ifndef CHTOperator2D_H
#define CHTOperator2D_H

/*---------------------------------------------------------------------------*\
    CHTOperator2D - 2D conjugate heat transfer as a framework Operator.

    Couples a 2D fluid Domain (bottom face = fluid-solid interface) to a 2D
    solid wall, solved column-by-column with the existing SolidEnergyEqn
    (1D implicit Thomas). This "strip-based" approach is exact when lateral
    conduction in the solid is negligible compared to normal conduction, which
    is a good approximation for thin rocket-chamber liners.

    Per step:
      1. Read fluid T and transport props at the first interior row (j=1).
      2. Compute convective heat flux q = h*(T_fluid - T_wall_top).
      3. Solve each column of the solid with SolidEnergyEqn (top: Convection BC
         from fluid, bottom: user-specified outer BC — typically FixedT or Adiabatic).
      4. Remove the same interface heat from the first fluid interior cell.

    The fluid domain's bottom halo (j=0) remains boundary-owned and is not the
    solid interface temperature. Solid state participates in the same accepted
    step transaction as the fluid energy update.

    Parameters:
      Nsolid   — number of solid cells in y (through-thickness)
      Lsolid   — solid thickness [m]
      k_solid  — thermal conductivity [W/(m K)]
      rc_solid — rho*cp of solid [J/(m3 K)]
      h_conv   — convective heat transfer coefficient [W/(m2 K)] at fluid-solid face
      bcOuter  — SolidBC on the external (cold) side of the solid
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "Domain.H"
#include "TransportProps.H"
#include "SpeciesRegistry.H"
#include "SolidEnergyEqn.H"
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace Foam
{

class CHTOperator2D : public Operator
{
private:
    int    Nsolid_;
    double Lsolid_, k_solid_, rc_solid_;
    double rho_solid_, cp_solid_, T0_solid_;
    double h_conv_;    // user-supplied convective HTC [W/(m2 K)]
    SolidBC bcOuter_;
    const SpeciesRegistry& reg_;

    mutable std::vector<std::vector<double>> solid_T_;
    mutable bool initialized_;

    // Conservation bookkeeping [J/m2 of interface, per unit depth]:
    //   qToSolid_   = cumulative heat that left the fluid into the solid
    //   qOutOuter_  = cumulative heat that left the solid at its outer BC
    // Interface energy balance check: Δ(solid stored) == qToSolid_ - qOutOuter_.
    mutable double qToSolid_;
    mutable double qOutOuter_;

    // Whole-step transaction state. The solid field and conservation ledger
    // advance together with Domain::U/W and must be restored on a rejected
    // step (including an RK2 stage-2 failure).
    mutable std::vector<std::vector<double>> solid_TStepStart_;
    mutable bool initializedStepStart_ = false;
    mutable double qToSolidStepStart_ = 0.0;
    mutable double qOutOuterStepStart_ = 0.0;
    mutable bool stepTransactionActive_ = false;

public:
    CHTOperator2D(int Nsolid, double Lsolid,
                  double rho_solid, double cp_solid, double k_solid,
                  double h_conv,
                  const SolidBC& bcOuter,
                  const SpeciesRegistry& reg,
                  double T0_solid = 300.0)
        : Nsolid_(Nsolid), Lsolid_(Lsolid),
          k_solid_(k_solid), rc_solid_(rho_solid*cp_solid),
          h_conv_(h_conv), bcOuter_(bcOuter),
          rho_solid_(rho_solid), cp_solid_(cp_solid),
          reg_(reg), initialized_(false), T0_solid_(T0_solid),
          qToSolid_(0.0), qOutOuter_(0.0)
    {}

    // ---- Operator interface ----
    bool ownsPersistentState() const override { return true; }

    void beginStepTransaction(const Domain&) const override
    {
        if (stepTransactionActive_)
            throw std::logic_error("nested CHT step transaction");
        solid_TStepStart_ = solid_T_;
        initializedStepStart_ = initialized_;
        qToSolidStepStart_ = qToSolid_;
        qOutOuterStepStart_ = qOutOuter_;
        stepTransactionActive_ = true;
    }

    void acceptStepTransaction() const noexcept override
    {
        solid_TStepStart_.clear();
        initializedStepStart_ = false;
        qToSolidStepStart_ = 0.0;
        qOutOuterStepStart_ = 0.0;
        stepTransactionActive_ = false;
    }

    void rollbackStepTransaction() const noexcept override
    {
        if (!stepTransactionActive_) return;
        solid_T_.swap(solid_TStepStart_);
        initialized_ = initializedStepStart_;
        qToSolid_ = qToSolidStepStart_;
        qOutOuter_ = qOutOuterStepStart_;
        solid_TStepStart_.clear();
        stepTransactionActive_ = false;
    }

    void averageWithStepStartState() const override
    {
        if (!stepTransactionActive_)
            throw std::logic_error("CHT RK2 state transaction is inactive");

        bool compatibleStart = initializedStepStart_
            && solid_TStepStart_.size() == solid_T_.size();
        if (compatibleStart)
        {
            for (std::size_t i = 0; i < solid_T_.size(); ++i)
                if (solid_TStepStart_[i].size() != solid_T_[i].size())
                {
                    compatibleStart = false;
                    break;
                }
        }

        // Lazy initialization and a changed column count both reset the solid
        // to T0_solid_ in apply(); use that same physical start state for RK2.
        for (std::size_t i = 0; i < solid_T_.size(); ++i)
            for (std::size_t k = 0; k < solid_T_[i].size(); ++k)
            {
                const double startT = compatibleStart
                    ? solid_TStepStart_[i][k] : T0_solid_;
                solid_T_[i][k] = 0.5*(startT + solid_T_[i][k]);
            }

        qToSolid_ = 0.5*(qToSolidStepStart_ + qToSolid_);
        qOutOuter_ = 0.5*(qOutOuterStepStart_ + qOutOuter_);
    }

    void apply(Domain& d, double dt) const override
    {
        if (d.Ny() < 2) return;   // 1D domain: CHT2D is a no-op

        const int Nx = d.Nx();

        // Lazy-init: size solid arrays to Nx columns
        if (!initialized_ || (int)solid_T_.size() != Nx)
        {
            solid_T_.assign(Nx, std::vector<double>(Nsolid_, 300.0));
            // Initialise solid T from the fluid domain's bottom halo T
            // Initialise the solid at its TRUE initial temperature T0_solid_,
            // NOT the fluid halo (which starts at the hot fluid T). Using the
            // halo would inject spurious stored energy into the conservation
            // ledger (solidStoredDelta is measured relative to T0_solid_).
            for (int i = 0; i < Nx; ++i)
                solid_T_[i].assign(Nsolid_, T0_solid_);
            initialized_ = true;
        }

        const double dy = d.dy();

        // For each column: solve 1D solid with fluid convection BC on top,
        // then REMOVE the exact same heat from the fluid interface cell so the
        // fluid+solid system conserves energy (flux-conservative coupling).
        //
        // ORDER CONTRACT (see Operator.H): CHT ACCUMULATES into Unew.rhoE, so it
        // MUST be added to the solver AFTER ConvectionOperator (which assigns).
        for (int i = 0; i < Nx; ++i)
        {
            // Fluid-side state at the interface (j=1, first interior row)
            const PrimitiveVars& wf = d.W(i, 1);
            double T_fluid = wf.T;

            // Top solid BC: convection from fluid with user HTC. The convection
            // conductance is a SERIES resistance (film + half-cell conduction),
            // so the flux the solid actually receives is consistent with bcTop.
            SolidBC bcTop = SolidBC::convection(h_conv_, T_fluid);

            // Step the 1D solid column (SolidEnergyEqn API: N,L,rho,cp,k,T0)
            SolidEnergyEqn solid(Nsolid_, Lsolid_, rho_solid_, cp_solid_, k_solid_, T0_solid_);
            solid.setField(solid_T_[i]);   // load current temperature profile
            solid.step(dt, bcTop, bcOuter_);

            // Interface heat flux INTO the solid this step [W/m2] (fluid loses it)
            double q_int = solid.leftFlux(bcTop);
            // Heat leaving the solid at its outer BC [W/m2] (positive = out)
            double q_out = -solid.rightFlux(bcOuter_);

            for (int k = 0; k < Nsolid_; ++k) solid_T_[i][k] = solid.T(k);

            // CONSERVATION: remove q_int from the fluid interface cell (j=1).
            // Volumetric energy sink [J/m3] over this cell of height dy:
            //   ΔrhoE = -q_int * dt / dy
            // This is the mirror of ChemistrySourceOperator's explicit rhoE add.
            d.Unew(i, 1).rhoE -= q_int * dt / dy;

            // Bookkeeping: TRUE energy per unit depth [J/m], summed over columns.
            // Face length per column = dx; flux [W/m2] * dx [m] * dt [s] = J/m.
            const double dxCol = d.dx();
            qToSolid_  += q_int * dxCol * dt;
            qOutOuter_ += q_out * dxCol * dt;
        }
    }

    // Stability limit for the EXPLICIT interface heat extraction (the solid
    // solve itself is implicit/unconditionally stable, but the explicit rhoE
    // sink on the thin fluid cell j=1 is NOT). Cap the per-step energy removed
    // from the interface cell to a fraction of its sensible internal energy:
    //   dt_cht = safety * (rho * e_sens * dy) / q_int      per column, take min.
    // Mirrors ChemistrySourceOperator's tau_chem limiter. Without it, a thin
    // cell + high HTC drives rhoE negative -> T inversion fails -> runaway.
    double maxStableDt(const Domain& d, double CFL) const override
    {
        (void)CFL;
        if (d.Ny() < 2) return 1e30;
        const int Nx = d.Nx();
        const double dy = d.dy();
        const double safety = 0.2;
        double dtmin = 1e30;

        for (int i = 0; i < Nx; ++i)
        {
            const PrimitiveVars& wf = d.W(i, 1);
            double T_fluid = wf.T;
            // Wall-top T: use the stored solid profile if initialised, else the
            // solid's INITIAL temperature (first step, before apply has run).
            // NOT the halo T — the halo is seeded to the fluid T, which would
            // give a spurious zero flux and skip the dt limit on step 1.
            double T_wall = (initialized_ && (int)solid_T_.size() == Nx)
                          ? solid_T_[i][0] : T0_solid_;

            // Convective flux magnitude into the solid [W/m2] with the series
            // resistance (film + half-cell conduction), matching SolidEnergyEqn.
            double cellThk = Lsolid_ / Nsolid_;
            double g = 1.0 / (1.0/h_conv_ + 0.5*cellThk/k_solid_);
            double q_int = std::abs(g * (T_fluid - T_wall));
            if (q_int < 1.0) continue;   // negligible flux -> no dt limit

            // Sensible internal energy density of the fluid cell [J/m3].
            const thermo::EquilibriumState& q = d.Q(i,1);
            double rhoE_sens = q.rho * std::abs(q.e);
            if (rhoE_sens < 1.0) continue;

            // energy removed this step [J/m3] = q_int*dt/dy ; cap to safety*rhoE_sens
            double dt_i = safety * rhoE_sens * dy / q_int;
            dtmin = std::min(dtmin, dt_i);
        }
        return dtmin;
    }

    std::string name() const override { return "CHT2D (fluid-solid coupling)"; }

    // Solid temperature accessor for post-processing
    const std::vector<double>& solidT(int col) const { return solid_T_[col]; }

    // Mean interface (wall) temperature
    double meanWallT() const
    {
        double sum = 0.0;
        for (const auto& col : solid_T_) sum += col[0];
        return solid_T_.empty() ? 0.0 : sum / solid_T_.size();
    }

    // ---- conservation diagnostics [J/m2 of interface, per unit depth] ----
    double heatToSolid()  const { return qToSolid_;  }  // cumulative fluid -> solid
    double heatOutOuter() const { return qOutOuter_; }  // cumulative solid -> outer BC

    // Solid stored-energy CHANGE relative to the uniform initial state T0_solid_,
    // summed over all columns [J/m2 per unit depth]. dxCol = interface length per
    // column (fluid dx); each solid cell has thickness Lsolid_/Nsolid_.
    double solidStoredDelta(double dxCol) const
    {
        const double cellThk = Lsolid_ / Nsolid_;
        double dE = 0.0;
        for (const auto& col : solid_T_)
            for (int k = 0; k < Nsolid_; ++k)
                // TRUE energy per unit depth [J/m]: rc [J/m3K] * dT [K]
                //   * cellThk [m] (through-thickness) * dxCol [m] (column width).
                dE += rc_solid_ * (col[k] - T0_solid_) * cellThk * dxCol;
        return dE;
    }

    // Interface energy-balance residual [J/m2]: (heat in) - (heat out) - (stored).
    // Should be ~0 to the implicit solver's accuracy (this is the Stage-5 check).
    double conservationResidual(double dxCol) const
    {
        return qToSolid_ - qOutOuter_ - solidStoredDelta(dxCol);
    }

    // ---- Stage 5: q_wall(x) profile for post-processing / PinFin input ----

    // Local instantaneous wall heat flux [W/m²] at column i.
    // Uses the same series-resistance formula as maxStableDt to stay consistent.
    // q > 0 means heat flows FROM fluid INTO solid (normal hot-gas-cooled wall).
    double qWall(int col, const Domain& d) const
    {
        if (!initialized_ || col >= (int)solid_T_.size()) return 0.0;
        const PrimitiveVars& wf = d.W(col, 1);
        double T_fluid = wf.T;
        double T_wall  = solid_T_[col][0];
        // series conductance: film (h_conv) + half-cell solid conduction
        double cellThk = Lsolid_ / Nsolid_;
        double g = 1.0 / (1.0/h_conv_ + 0.5*cellThk/k_solid_);
        return g * (T_fluid - T_wall);
    }

    // q_wall profile over all columns [W/m²] as a vector indexed by column i.
    std::vector<double> qWallProfile(const Domain& d) const
    {
        const int Nx = d.Nx();
        std::vector<double> qw(Nx, 0.0);
        for (int i = 0; i < Nx; ++i) qw[i] = qWall(i, d);
        return qw;
    }

    // T_wall profile over all columns [K].
    std::vector<double> wallTProfile() const
    {
        std::vector<double> tw;
        tw.reserve(solid_T_.size());
        for (const auto& col : solid_T_) tw.push_back(col.empty() ? T0_solid_ : col[0]);
        return tw;
    }

    // Write q_wall(x) and T_wall(x) to a CSV for inspection.
    // Format: x[m], q_wall[W/m2], T_wall[K]
    void writeQwallCSV(const std::string& filename, const Domain& d) const
    {
        std::ofstream f(filename);
        f << "# x[m],  q_wall[W/m2],  T_wall[K]  (generated by CHTOperator2D)\n";
        f << "# T_fluid_inlet=" << (d.W(0,1).T) << "K\n";
        const int Nx = d.Nx();
        for (int i = 0; i < Nx; ++i)
        {
            double x = d.x(i, 1);
            double q = qWall(i, d);
            double Tw = solid_T_.empty() ? T0_solid_ : (i<(int)solid_T_.size() ? solid_T_[i][0] : T0_solid_);
            f << std::setprecision(6) << x << "," << q << "," << Tw << "\n";
        }
    }

    // Write a PinFin-Designer-compatible JSON (Stage 5 handoff).
    // Schema matches pinfin_designer/database/schema.py:
    //   x_positions   [mm]   (RFAM stores x in m — converted here)
    //   q_wall_target [W/m2]
    //   T_wall        [K]    (extra field; PinFin recomputes its own)
    // PinFin's Designer.load_heat_flux(x_mm, q_wall) consumes x_positions
    // and q_wall_target directly.
    void writeQwallJSON(const std::string& filename, const Domain& d) const
    {
        const int Nx = d.Nx();
        std::ofstream f(filename);
        f << std::setprecision(8);
        f << "{\n";
        f << "  \"source\": \"realFluidAllMach CHTOperator2D\",\n";
        f << "  \"units\": {\"x_positions\": \"mm\", \"q_wall_target\": \"W/m2\", \"T_wall\": \"K\"},\n";
        f << "  \"T_fluid_inlet_K\": " << d.W(0,1).T << ",\n";

        // x_positions [mm]
        f << "  \"x_positions\": [";
        for (int i = 0; i < Nx; ++i)
            f << (i?", ":"") << d.x(i,1)*1000.0;   // m -> mm
        f << "],\n";

        // q_wall_target [W/m2]
        f << "  \"q_wall_target\": [";
        for (int i = 0; i < Nx; ++i)
            f << (i?", ":"") << qWall(i, d);
        f << "],\n";

        // T_wall [K]
        f << "  \"T_wall\": [";
        for (int i = 0; i < Nx; ++i)
        {
            double Tw = (i < (int)solid_T_.size() && !solid_T_[i].empty())
                        ? solid_T_[i][0] : T0_solid_;
            f << (i?", ":"") << Tw;
        }
        f << "]\n";
        f << "}\n";
    }
};

} // End namespace Foam

#endif
