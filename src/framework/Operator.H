#ifndef Operator_H
#define Operator_H

/*---------------------------------------------------------------------------*\
    Operator - pluggable physics applied each time-step.

    Each operator READS the t^n state (domain.U / domain.W) and WRITES or
    ACCUMULATES increments into domain.Unew. The integrator's additive update
    (beginStep -> apply ops -> commit) therefore reproduces solver1d_fpv's
    single synchronous update exactly.

    maxStableDt() lets each operator contribute its own stability limit;
    SolverBase takes the min over all operators.

    Operators:
      - ConvectionOperator : AUSM finite volume, interior i=1..N-2.
                             Optional experimental acoustic scaling for the
                             CFL estimate only. It does not precondition the
                             flux residual, time derivative, or Jacobian.
      - FPVSourceOperator  : Z/C diffusion + FPV table source on rhoC
                             (reproduces solver1d_fpv source step).
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "Residual.H"
#include "AcousticStateAdmission.H"
#include "FluxScheme.H"
#include "FPVTable.H"
#include "RealFluidState.H"
#include "ExperimentalCflAcousticScaler.H"
#include "MUSCLReconstruction.H"
#include "StructuredImplicitBlock.H"
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace Foam
{

// Positive scalar bounds for a structured approximate Jacobian. The diagonal
// stores local rate bounds [1/s]; each face entry stores the magnitude of the
// corresponding negative off-diagonal coupling. LU-SGS may use these bounds
// without differentiating a noisy thermo/chemistry residual. Operators that
// cannot provide trustworthy neighbor coupling fall back to a diagonal bound.
struct StructuredImplicitStencil
{
    enum Face { I_MINUS=0,I_PLUS,J_MINUS,J_PLUS,K_MINUS,K_PLUS,FACE_COUNT };
    std::vector<double> diagonal;
    std::vector<double> face[FACE_COUNT];

    void ensureSize(std::size_t count)
    {
        if(diagonal.size()!=count) diagonal.assign(count,0.0);
        for(auto& values:face)
            if(values.size()!=count) values.assign(count,0.0);
    }
};

class Operator
{
public:
    virtual ~Operator() = default;
    // Persistent model state outside Domain participates in the same accepted
    // or rejected step as U/W.  Diagnostic counters and reconstructible
    // workspaces do not need to opt in.
    virtual bool ownsPersistentState() const { return false; }
    virtual void beginStepTransaction(const Domain&) const {}
    virtual void acceptStepTransaction() const noexcept {}
    virtual void rollbackStepTransaction() const noexcept {}
    virtual void averageWithStepStartState() const {}
    // CONTRACT (additive-update, see TimeIntegrator):
    //   - READ the t^n state via d.U(i)/d.W(i) ONLY.
    //   - WRITE results into d.Unew(i). An operator may assign (=) or accumulate
    //     (+=) into Unew, but MUST NOT clobber another operator's contribution.
    //     ConvectionOperator ASSIGNS interior cells (it is the base transport
    //     step); source operators (FPV, Chemistry) ACCUMULATE (+=). Therefore
    //     Convection MUST be added BEFORE any source operator in the chain, or
    //     the assignment wipes the source increment. (A future residual-based
    //     integrator would remove this ordering coupling by summing residuals.)
    virtual void apply(Domain& d, double dt) const = 0;
    // RESIDUAL CONTRACT (see Residual.H and ResidualIntegrator):
    //   - READ the evaluation state via d.U/W/Q; must not write d.Unew.
    //   - ACCUMULATE (+=) the spatial balance div(F) - S into R.
    //   - Must be a pure function of the current state: called more than once
    //     per step (inner iteration, line search) it must return the same R.
    // Operators that cannot satisfy this (implicit sub-solves that depend on
    // dt) must say so through supportsResidual()==false rather than silently
    // returning a dt-averaged rate.
    virtual bool supportsResidual() const { return false; }
    virtual void residual(const Domain& d, ResidualField& R) const
    {
        (void)d; (void)R;
        throw std::logic_error(
            name() + ": residual() is not implemented; this operator can only "
            "run on the additive (apply) path");
    }
    virtual double maxStableDt(const Domain& d, double CFL) const
    { (void)d; (void)CFL; return 1e30; }
    // Optional geometry-aware per-cell stability field for steady/pseudo-time
    // preconditioning. The default preserves legacy global-step behavior.
    virtual void localStableDt(
        const Domain& d,double CFL,std::vector<double>& dt) const
    {
        dt.assign(static_cast<std::size_t>(d.N()),maxStableDt(d,CFL));
    }
    // Add a scalar structured approximate-Jacobian bound for LU-SGS. The
    // default deliberately contributes diagonal damping only: inventing face
    // coupling for a non-smooth or black-box operator would be less robust
    // than omitting it. A finite local stability step supplies that bound.
    virtual void addStructuredImplicitStencil(
        const Domain& d,StructuredImplicitStencil& stencil) const
    {
        stencil.ensureSize(static_cast<std::size_t>(d.N()));
        std::vector<double> dt;
        localStableDt(d,1.0,dt);
        if(dt.size()!=static_cast<std::size_t>(d.N()))
            throw std::runtime_error(
                "operator returned invalid local stability field: "+name());
        for(int cell=0;cell<d.N();++cell)
        {
            const double value=dt[static_cast<std::size_t>(cell)];
            if(std::isfinite(value)&&value>0.0&&value<1e29)
                stencil.diagonal[static_cast<std::size_t>(cell)]+=1.0/value;
        }
    }
    // Operators with a trustworthy analytic frozen block may add the delta
    // between that block action and their already-advertised scalar neighbor
    // bound during LU-SGS sweeps. Returning false requests a per-face fallback
    // to the scalar bound, so no trusted coupling is ever silently removed.
    virtual bool supportsStructuredBlockCoupling() const { return false; }
    virtual bool supportsStructuredDiagonalBlock() const { return false; }
    // Optional local hydro diagonal block. The returned contribution is the
    // operator's physical diagonal block, while scalarContribution is the
    // scalar diagonal already advertised by addStructuredImplicitStencil.
    // The integrator replaces that scalar contribution with this block and
    // keeps all unsupported operators in the scalar diagonal.
    virtual bool addStructuredImplicitDiagonalBlock(
        const Domain& d,int cell,
        StructuredImplicitBlock& contribution,
        double& scalarContribution) const
    {
        (void)d;(void)cell;(void)contribution;(void)scalarContribution;
        return false;
    }
    virtual bool addStructuredImplicitNeighborContribution(
        const Domain& d,int cell,int neighbor,int face,
        const ConservedVars& neighborCorrection,
        ConservedVars& contribution) const
    {
        (void)d;(void)cell;(void)neighbor;(void)face;
        (void)neighborCorrection;(void)contribution;
        return false;
    }
    virtual std::string name() const = 0;
};

// --- Convection: AUSM finite-volume on interior cells (== solver1d_fpv/test) ---
class ConvectionOperator : public Operator
{
private:
    std::unique_ptr<FluxScheme> flux_;
    bool useExperimentalCflScaler_;
    bool useMUSCL_;    // Stage 7+: second-order spatial reconstruction (opt-in)
public:
    explicit ConvectionOperator(std::unique_ptr<FluxScheme> f,
                                bool useExperimentalCflScaler = false,
                                bool useMUSCL = false)
        : flux_(std::move(f)),
          useExperimentalCflScaler_(useExperimentalCflScaler),
          useMUSCL_(useMUSCL)
    {}

    void apply(Domain& d, double dt) const override
    {
        requireAdmittedAcoustics(d,"apply");
        flux_->setSpeciesNames(d.speciesNames());
        if (useMUSCL_ && (d.Ny() > 1 || d.Nz() > 1))
            throw std::runtime_error(
                "ConvectionOperator: MUSCL is implemented only for 1D; "
                "multidimensional reconstruction must not degrade to first order");
        if (useMUSCL_
            && !d.thermoModel().supportsPrimitivePTYInitialization())
            throw std::runtime_error(
                "ConvectionOperator: phase-aware MUSCL requires a complete "
                "face-state rho-e-Y recovery");
        if (d.Nz() > 1)
            apply3D(d, dt);
        else if (d.Ny() > 1)
            apply2D(d, dt);
        else
            apply1D(d, dt);
    }

    // ---- residual path (ResidualIntegrator) ----
    // Avoid inconsistent reconstruction on shared faces at the legacy
    // first-order/MUSCL transition until that path is made face-based.
    bool supportsResidual() const override { return !useMUSCL_; }

    void residual(const Domain& d, ResidualField& R) const override
    {
        if (!supportsResidual())
            throw std::logic_error("ConvectionOperator: residual requires first-order face states");
        requireAdmittedAcoustics(d,"residual");
        flux_->setSpeciesNames(d.speciesNames());
        if (R.size() != static_cast<std::size_t>(d.N()))
            throw std::runtime_error(
                "ConvectionOperator: residual field size does not match Domain");
        if (d.Nz() > 1)
            residual3D(d, R);
        else if (d.Ny() > 1)
            residual2D(d, R);
        else
            residual1D(d, R);
    }

private:
    static void requireAdmittedAcoustics(
        const Domain& d,
        const char* operation)
    {
        for (int cell = 0; cell < d.N(); ++cell)
            if (!solverConsumerAdmitsAcousticState(d.Q(cell)))
                throw std::runtime_error(
                    std::string("ConvectionOperator: ") + operation
                    + " rejected acoustic state at cell "
                    + std::to_string(cell));
    }

    thermo::EquilibriumState initializeFaceThermo(
        const PrimitiveVars& face,
        const Domain& d) const
    {
        thermo::EquilibriumResult initialized =
            d.thermoModel().initializePTY(face.p,face.T,face.Y);
        if (!initialized.success())
            throw std::runtime_error(
                "ConvectionOperator: reconstructed face thermodynamics failed");
        return std::move(*initialized.state);
    }

    // ---- 1D face fluxes: the single source of truth for both paths ----
    // useMUSCL_=false: first-order (bit-identical MVP path).
    // useMUSCL_=true:  MUSCL minmod 2nd-order reconstruction on primitives.
    void faceFluxes1D(int i, const Domain& d,
                      ConservedVars& Fl, ConservedVars& Fr) const
    {
        static const std::vector<double> normal{1.0, 0.0, 0.0};
        MUSCLReconstruction muscl(Limiter::MINMOD);
        {
            if (useMUSCL_ && i >= 2 && i < d.N()-2)
            {
                // MUSCL: reconstruct face states from 4-point stencil
                // Left face (i-1/2): L+ from cell i-1 extrapolated right
                PrimitiveVars wL_left = d.W(i-1);
                wL_left.p    = muscl.leftFace(d.W(i-2).p,    d.W(i-1).p,    d.W(i).p);
                wL_left.T    = muscl.leftFace(d.W(i-2).T,    d.W(i-1).T,    d.W(i).T);
                wL_left.U[0] = muscl.leftFace(d.W(i-2).U[0], d.W(i-1).U[0], d.W(i).U[0]);
                // R- from cell i extrapolated left
                PrimitiveVars wR_left = d.W(i);
                wR_left.p    = muscl.rightFace(d.W(i-1).p,    d.W(i).p,    d.W(i+1).p);
                wR_left.T    = muscl.rightFace(d.W(i-1).T,    d.W(i).T,    d.W(i+1).T);
                wR_left.U[0] = muscl.rightFace(d.W(i-1).U[0], d.W(i).U[0], d.W(i+1).U[0]);

                // Right face (i+1/2): L+ from cell i extrapolated right
                PrimitiveVars wL_right = d.W(i);
                wL_right.p    = muscl.leftFace(d.W(i-1).p,    d.W(i).p,    d.W(i+1).p);
                wL_right.T    = muscl.leftFace(d.W(i-1).T,    d.W(i).T,    d.W(i+1).T);
                wL_right.U[0] = muscl.leftFace(d.W(i-1).U[0], d.W(i).U[0], d.W(i+1).U[0]);
                // R- from cell i+1 extrapolated left
                PrimitiveVars wR_right = d.W(i+1);
                wR_right.p    = muscl.rightFace(d.W(i).p,    d.W(i+1).p,    d.W(i+2).p);
                wR_right.T    = muscl.rightFace(d.W(i).T,    d.W(i+1).T,    d.W(i+2).T);
                wR_right.U[0] = muscl.rightFace(d.W(i).U[0], d.W(i+1).U[0], d.W(i+2).U[0]);

                const thermo::EquilibriumState qLLeft =
                    initializeFaceThermo(wL_left,d);
                const thermo::EquilibriumState qRLeft =
                    initializeFaceThermo(wR_left,d);
                const thermo::EquilibriumState qLRight =
                    initializeFaceThermo(wL_right,d);
                const thermo::EquilibriumState qRRight =
                    initializeFaceThermo(wR_right,d);
                flux_->computeFluxFromEquilibrium(
                    wL_left,qLLeft,wR_left,qRLeft,normal,Fl);
                flux_->computeFluxFromEquilibrium(
                    wL_right,qLRight,wR_right,qRRight,normal,Fr);
            }
            else
            {
                // First-order (legacy, bit-identical MVP)
                flux_->computeFluxFromEquilibrium(
                    d.W(i-1),d.Q(i-1),d.W(i),d.Q(i),normal,Fl);
                flux_->computeFluxFromEquilibrium(
                    d.W(i),d.Q(i),d.W(i+1),d.Q(i+1),normal,Fr);
            }
        }
    }

    // ---- 1D interior update (additive path, arithmetic unchanged) ----
    void apply1D(Domain& d, double dt) const
    {
        const int N = d.N();
        const double dx = d.dx();

        for (int i = 1; i < N-1; ++i)
        {
            ConservedVars Fl(d.nSpecies()), Fr(d.nSpecies());
            faceFluxes1D(i, d, Fl, Fr);

            const ConservedVars& u = d.U(i);
            ConservedVars& un = d.Unew(i);
            un.rho  = u.rho  - dt/dx * (Fr.rho  - Fl.rho);
            un.rhoE = u.rhoE - dt/dx * (Fr.rhoE - Fl.rhoE);
            for (int j = 0; j < 3; ++j)
                un.rhoU[j] = u.rhoU[j] - dt/dx * (Fr.rhoU[j] - Fl.rhoU[j]);
            for (size_t j = 0; j < u.rhoY.size(); ++j)
                un.rhoY[j] = u.rhoY[j] - dt/dx * (Fr.rhoY[j] - Fl.rhoY[j]);
            un.rhoZ = u.rhoZ - dt/dx * (Fr.rhoZ - Fl.rhoZ);
            un.rhoC = u.rhoC - dt/dx * (Fr.rhoC - Fl.rhoC);
        }
    }

    // ---- 1D spatial residual ----
    void residual1D(const Domain& d, ResidualField& R) const
    {
        const int N = d.N();
        const double invDx = 1.0/d.dx();
        for (int i = 1; i < N-1; ++i)
        {
            ConservedVars Fl(d.nSpecies()), Fr(d.nSpecies());
            faceFluxes1D(i, d, Fl, Fr);
            addFluxDivergence(R[i], Fr, Fl, invDx);
        }
    }

    // ---- 2D face fluxes: shared by apply2D and residual2D ----
    void faceFluxes2D(int i, int j, const Domain& d,
                      ConservedVars& FxL, ConservedVars& FxR,
                      ConservedVars& FyL, ConservedVars& FyR) const
    {
        static const std::vector<double> nx{1.0,0.0,0.0}, ny{0.0,1.0,0.0};
        flux_->computeFluxFromEquilibrium(
            d.W(i-1,j),d.Q(i-1,j),d.W(i,j),d.Q(i,j),nx,FxL);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j),d.Q(i,j),d.W(i+1,j),d.Q(i+1,j),nx,FxR);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j-1),d.Q(i,j-1),d.W(i,j),d.Q(i,j),ny,FyL);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j),d.Q(i,j),d.W(i,j+1),d.Q(i,j+1),ny,FyR);
    }

    // ---- 2D interior update: Cartesian FV (dimensionally split in-place) ----
    // Each interior cell (i=1..Nx-2, j=1..Ny-2) accumulates x-fluxes and
    // y-fluxes via the same AUSM scheme with different face normals.
    // Boundary cells (halo layer i=0/Nx-1, j=0/Ny-1) are left at Unew = U
    // (beginStep copies U); the calling solver applies BCs separately.
    void apply2D(Domain& d, double dt) const
    {
        const int  Nx = d.Nx(), Ny = d.Ny();
        const int  nS = d.nSpecies();
        const double dx = d.dx(), dy = d.dy();

        for (int i = 1; i < Nx-1; ++i)
        for (int j = 1; j < Ny-1; ++j)
        {
            ConservedVars FxL(nS), FxR(nS), FyL(nS), FyR(nS);
            faceFluxes2D(i,j,d,FxL,FxR,FyL,FyR);

            const ConservedVars& u = d.U(i,j);
            ConservedVars& un = d.Unew(i,j);

            // FV update: Unew = U - dt*(dFx/dx + dFy/dy)
            un.rho  = u.rho  - dt*(FxR.rho -FxL.rho )/dx - dt*(FyR.rho -FyL.rho )/dy;
            un.rhoE = u.rhoE - dt*(FxR.rhoE-FxL.rhoE)/dx - dt*(FyR.rhoE-FyL.rhoE)/dy;
            for (int k = 0; k < 3; ++k)
                un.rhoU[k] = u.rhoU[k] - dt*(FxR.rhoU[k]-FxL.rhoU[k])/dx
                                        - dt*(FyR.rhoU[k]-FyL.rhoU[k])/dy;
            for (int k = 0; k < nS; ++k)
                un.rhoY[k] = u.rhoY[k] - dt*(FxR.rhoY[k]-FxL.rhoY[k])/dx
                                        - dt*(FyR.rhoY[k]-FyL.rhoY[k])/dy;
            un.rhoZ = u.rhoZ - dt*(FxR.rhoZ-FxL.rhoZ)/dx - dt*(FyR.rhoZ-FyL.rhoZ)/dy;
            un.rhoC = u.rhoC - dt*(FxR.rhoC-FxL.rhoC)/dx - dt*(FyR.rhoC-FyL.rhoC)/dy;
        }
    }

    // ---- 2D spatial residual ----
    void residual2D(const Domain& d, ResidualField& R) const
    {
        const int  Nx = d.Nx(), Ny = d.Ny(), nS = d.nSpecies();
        const double invDx = 1.0/d.dx(), invDy = 1.0/d.dy();
        for (int i = 1; i < Nx-1; ++i)
        for (int j = 1; j < Ny-1; ++j)
        {
            ConservedVars FxL(nS), FxR(nS), FyL(nS), FyR(nS);
            faceFluxes2D(i,j,d,FxL,FxR,FyL,FyR);
            ConservedVars& r = R[d.index(i,j)];
            addFluxDivergence(r, FxR, FxL, invDx);
            addFluxDivergence(r, FyR, FyL, invDy);
        }
    }

    // ---- 3D face fluxes: shared by apply3D and residual3D ----
    void faceFluxes3D(int i,int j,int k,const Domain& d,
                      ConservedVars& FxL,ConservedVars& FxR,
                      ConservedVars& FyL,ConservedVars& FyR,
                      ConservedVars& FzL,ConservedVars& FzR) const
    {
        static const std::vector<double> nx{1.0,0.0,0.0}, ny{0.0,1.0,0.0}, nz{0.0,0.0,1.0};
        flux_->computeFluxFromEquilibrium(
            d.W(i-1,j,k),d.Q(i-1,j,k),d.W(i,j,k),d.Q(i,j,k),nx,FxL);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j,k),d.Q(i,j,k),d.W(i+1,j,k),d.Q(i+1,j,k),nx,FxR);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j-1,k),d.Q(i,j-1,k),d.W(i,j,k),d.Q(i,j,k),ny,FyL);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j,k),d.Q(i,j,k),d.W(i,j+1,k),d.Q(i,j+1,k),ny,FyR);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j,k-1),d.Q(i,j,k-1),d.W(i,j,k),d.Q(i,j,k),nz,FzL);
        flux_->computeFluxFromEquilibrium(
            d.W(i,j,k),d.Q(i,j,k),d.W(i,j,k+1),d.Q(i,j,k+1),nz,FzR);
    }

    void apply3D(Domain& d, double dt) const
    {
        const int Nx=d.Nx(),Ny=d.Ny(),Nz=d.Nz(),nS=d.nSpecies();
        const double dx=d.dx(),dy=d.dy(),dz=d.dz();

        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<Nx-1;++i)
        for(int j=1;j<Ny-1;++j)
        for(int k=1;k<Nz-1;++k)
        {
            ConservedVars FxL(nS),FxR(nS),FyL(nS),FyR(nS),FzL(nS),FzR(nS);
            faceFluxes3D(i,j,k,d,FxL,FxR,FyL,FyR,FzL,FzR);
            const ConservedVars& u=d.U(i,j,k);
            ConservedVars& un=d.Unew(i,j,k);
            un.rho =u.rho -dt*(FxR.rho -FxL.rho )/dx-dt*(FyR.rho -FyL.rho )/dy-dt*(FzR.rho -FzL.rho )/dz;
            un.rhoE=u.rhoE-dt*(FxR.rhoE-FxL.rhoE)/dx-dt*(FyR.rhoE-FyL.rhoE)/dy-dt*(FzR.rhoE-FzL.rhoE)/dz;
            for(int m=0;m<3;++m)
                un.rhoU[m]=u.rhoU[m]-dt*(FxR.rhoU[m]-FxL.rhoU[m])/dx
                                     -dt*(FyR.rhoU[m]-FyL.rhoU[m])/dy
                                     -dt*(FzR.rhoU[m]-FzL.rhoU[m])/dz;
            for(int m=0;m<nS;++m)
                un.rhoY[m]=u.rhoY[m]-dt*(FxR.rhoY[m]-FxL.rhoY[m])/dx
                                     -dt*(FyR.rhoY[m]-FyL.rhoY[m])/dy
                                     -dt*(FzR.rhoY[m]-FzL.rhoY[m])/dz;
            un.rhoZ=u.rhoZ-dt*(FxR.rhoZ-FxL.rhoZ)/dx-dt*(FyR.rhoZ-FyL.rhoZ)/dy-dt*(FzR.rhoZ-FzL.rhoZ)/dz;
            un.rhoC=u.rhoC-dt*(FxR.rhoC-FxL.rhoC)/dx-dt*(FyR.rhoC-FyL.rhoC)/dy-dt*(FzR.rhoC-FzL.rhoC)/dz;
        }
    }

    // ---- 3D spatial residual ----
    void residual3D(const Domain& d, ResidualField& R) const
    {
        const int Nx=d.Nx(),Ny=d.Ny(),Nz=d.Nz(),nS=d.nSpecies();
        const double invDx=1.0/d.dx(),invDy=1.0/d.dy(),invDz=1.0/d.dz();

        // Keep exceptions inside the transaction's evaluation thread.
        for(int i=1;i<Nx-1;++i)
        for(int j=1;j<Ny-1;++j)
        for(int k=1;k<Nz-1;++k)
        {
            ConservedVars FxL(nS),FxR(nS),FyL(nS),FyR(nS),FzL(nS),FzR(nS);
            faceFluxes3D(i,j,k,d,FxL,FxR,FyL,FyR,FzL,FzR);
            ConservedVars& r = R[d.index(i,j,k)];
            addFluxDivergence(r, FxR, FxL, invDx);
            addFluxDivergence(r, FyR, FyL, invDy);
            addFluxDivergence(r, FzR, FzL, invDz);
        }
    }

public:
    double maxStableDt(const Domain& d, double CFL) const override
    {
        requireAdmittedAcoustics(d,"CFL scan");
        flux_->setSpeciesNames(d.speciesNames());
        double dtmin = 1e30;
        const int total = d.N();  // works for both 1D and 2D (Nx*Ny)
        const double dx = d.dx(), dy = (d.Ny()>1) ? d.dy() : d.dx();
        const double dz = (d.Nz()>1) ? d.dz() : d.dx();
        #pragma omp parallel for reduction(min:dtmin) schedule(static)
        for (int k = 0; k < total; ++k)
        {
            // Access via 1D index (linear storage, same for 1D and 2D)
            const PrimitiveVars& w = d.W(k);
            const thermo::EquilibriumState& q = d.Q(k);
            double ux = std::abs(w.U[0]);
            double uy = std::abs(w.U[1]);
            double uz = std::abs(w.U[2]);
            double a = q.acoustic.soundSpeed();
            if (useExperimentalCflScaler_)
            {
                ExperimentalCflAcousticScaler scaler;
                scaler.setEnabled(true);
                const double velocityMagnitude =
                    std::sqrt(w.U[0]*w.U[0] + w.U[1]*w.U[1]
                        + w.U[2]*w.U[2]);
                scaler.updateForMach(velocityMagnitude/a);
                a = scaler.cflAcousticSpeed(a, w.U);
            }
            double dt_x = CFL * dx / (ux + a + 1e-30);
            double dt_y = CFL * dy / (uy + a + 1e-30);
            double dt_z = CFL * dz / (uz + a + 1e-30);
            double localDt = dt_x;
            if (d.Ny() > 1) localDt = std::min(localDt,dt_y);
            if (d.Nz() > 1) localDt = std::min(localDt,dt_z);
            dtmin = std::min(dtmin,localDt);
        }
        return dtmin;
    }
    std::string name() const override
    {
        return useExperimentalCflScaler_
            ? "Convection (AUSM; experimental CFL scaler only)"
            : "Convection (AUSM)";
    }
};

// --- FPV source + Z/C diffusion (LEGACY REGRESSION BASELINE ONLY) ---
// DO NOT use in Vortex production solver chains.
// Production combustion → ChemistrySourceOperator (TwoStepCO + CVODE).
// This operator exists solely to preserve solver1d_fpv iter=293 baseline.
// Arch note: FPV uses (Z, C) variables incompatible with rhoY_i finite-rate path.
class FPVSourceOperator : public Operator
{
private:
    const FPVTable* table_;
    double D_Z_, D_C_;
public:
    FPVSourceOperator(const FPVTable* table, double D_Z = 1e-5, double D_C = 1e-5)
        : table_(table), D_Z_(D_Z), D_C_(D_C) {}

    void apply(Domain& d, double dt) const override
    {
        const int N = d.N();
        const double dx = d.dx();
        for (int i = 1; i < N-1; ++i)
        {
            // Z diffusion (passive scalar, no source)
            double Zl = d.U(i-1).rhoZ / d.U(i-1).rho;
            double Zc = d.U(i).rhoZ   / d.U(i).rho;
            double Zr = d.U(i+1).rhoZ / d.U(i+1).rho;
            double d2Z = (Zr - 2.0*Zc + Zl) / (dx*dx);
            d.Unew(i).rhoZ += dt * d.U(i).rho * D_Z_ * d2Z;

            // C diffusion + FPV source
            double Cl = d.U(i-1).rhoC / d.U(i-1).rho;
            double Cc = d.U(i).rhoC   / d.U(i).rho;
            double Cr = d.U(i+1).rhoC / d.U(i+1).rho;
            double d2C = (Cr - 2.0*Cc + Cl) / (dx*dx);

            RealFluidState s;
            s.updatePTY(d.W(i).p, d.W(i).T, d.W(i).Y);
            double h = s.h();
            double omega_C = table_->getOmegaDot(d.W(i).Z, d.W(i).C, h);

            d.Unew(i).rhoC += dt * (d.U(i).rho * D_C_ * d2C + omega_C);
        }
    }

    double maxStableDt(const Domain& d, double CFL) const override
    {
        (void)CFL;
        return 0.5 * d.dx() * d.dx() / std::max(D_Z_, D_C_);
    }
    std::string name() const override { return "FPV source + Z/C diffusion"; }
};

} // End namespace Foam

#endif
