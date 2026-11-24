#ifndef CVODEWrapper_H
#define CVODEWrapper_H

/*---------------------------------------------------------------------------*\
    CVODEWrapper - SUNDIALS CVODE interface for stiff chemistry ODEs

    Integrates the 0-D reactor state vector [Y_0..Y_{nS-1}, T] using CVODE:
      - BDF (Backward Differentiation Formula) for stiff systems
      - Newton iteration with a dense, difference-quotient Jacobian
      - Adaptive internal time stepping with built-in error control

    Drop-in stiff alternative to RosenbrockSolver. Constant-(p, rho) reactor:
        dY_i/dt = omega_i / rho          (species mass fractions)
        dT/dt   = qdot   / (rho * Cv)    (energy, ideal-Cv closure)

    NOTE: the ideal-gas Cv closure mirrors RosenbrockSolver for apples-to-apples
    0-D validation. Real-fluid T inversion (RealFluidState::updateRhoEY) belongs
    to the Stage-4 CFD coupling, NOT to the chemistry layer. Chemistry only ever
    supplies source terms (omega_i, qdot) here.

    Requires SUNDIALS >= 7.0 (SUNContext, SUN_COMM_NULL, sunrealtype).
    Link: -lsundials_cvode -lsundials_core -lsundials_nvecserial
\*---------------------------------------------------------------------------*/

#include "ChemistryModel.H"
#include <vector>
#include <stdexcept>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

namespace Foam
{

class CVODEWrapper
{
private:
    ChemistryModel* chem_;
    int nS_;             // number of species
    int nEq_;            // nS + 1 (state = [Y..., T])

    double p_, rho_, Cv_;
    double rtol_, atol_;

    SUNContext      sunctx_;
    void*           cvode_mem_;
    N_Vector        y_;
    SUNMatrix       A_;
    SUNLinearSolver LS_;
    bool            initialized_;
    std::vector<double> rhsY_,rhsOmega_;

    // RHS: dy/dt = f(t, y) ; user_data == this
    static int rhsFn(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data)
    {
        CVODEWrapper* self = static_cast<CVODEWrapper*>(user_data);
        const int nS = self->nS_;

        const sunrealtype* yd = N_VGetArrayPointer(y);
        sunrealtype*       fd = N_VGetArrayPointer(ydot);

        // Unpack with the same safeguards RosenbrockSolver uses
        for (int i = 0; i < nS; ++i) self->rhsY_[i] = (yd[i] > 0.0 ? yd[i] : 0.0);
        double T = (yd[nS] > 100.0 ? yd[nS] : 100.0);

        double qdot;
        self->chem_->computeSourceTerms(T,self->p_,self->rho_,self->rhsY_,self->rhsOmega_,qdot);

        for (int i = 0; i < nS; ++i) fd[i] = self->rhsOmega_[i] / self->rho_;
        fd[nS] = qdot / (self->rho_ * self->Cv_);
        return 0;
    }

    void freeAll()
    {
        if (cvode_mem_) { CVodeFree(&cvode_mem_);    cvode_mem_ = nullptr; }
        if (LS_)        { SUNLinSolFree(LS_);        LS_        = nullptr; }
        if (A_)         { SUNMatDestroy(A_);         A_         = nullptr; }
        if (y_)         { N_VDestroy(y_);            y_         = nullptr; }
        if (sunctx_)    { SUNContext_Free(&sunctx_); sunctx_    = nullptr; }
        initialized_ = false;
    }

public:
    CVODEWrapper(ChemistryModel* chem, double Cv = 1500.0)
    :
        chem_(chem), nS_(chem->nSpecies()), nEq_(chem->nSpecies() + 1),
        p_(0.0), rho_(0.0), Cv_(Cv), rtol_(1e-6), atol_(1e-10),
        sunctx_(nullptr), cvode_mem_(nullptr),
        y_(nullptr), A_(nullptr), LS_(nullptr), initialized_(false),
        rhsY_(nS_,0.0),rhsOmega_(nS_,0.0)
    {}

    ~CVODEWrapper() { freeAll(); }

    void setTolerances(double rtol, double atol) { rtol_ = rtol; atol_ = atol; }

    // Update the ideal-gas Cv closure used in dT/dt = qdot/(rho*Cv).
    // Call BEFORE init() so the RHS function sees the new value.
    // Stage 4: FiniteRateReactor may call this with NASA-7 Cv_mix(T,Y)
    // to replace the placeholder 1500 J/(kg K).
    void setCv(double cv) { Cv_ = cv; }

    // Initialise the integrator at (t0, state). Constant-(p, rho) reactor.
    void init(double t0, double p, double rho,
              const std::vector<double>& Y0, double T0)
    {
        p_ = p; rho_ = rho;

        if (initialized_)
        {
            sunrealtype* yd = N_VGetArrayPointer(y_);
            for (int i = 0; i < nS_; ++i) yd[i] = Y0[i];
            yd[nS_] = T0;
            if (CVodeReInit(cvode_mem_, t0, y_))
                throw std::runtime_error("CVODEWrapper: CVodeReInit failed");
            CVodeSStolerances(cvode_mem_, rtol_, atol_);
            CVodeSetUserData(cvode_mem_, this);
            CVodeSetMaxNumSteps(cvode_mem_, 1000000);
            return;
        }

        if (SUNContext_Create(SUN_COMM_NULL, &sunctx_))
            throw std::runtime_error("CVODEWrapper: SUNContext_Create failed");

        y_ = N_VNew_Serial(nEq_, sunctx_);
        if (!y_) throw std::runtime_error("CVODEWrapper: N_Vector allocation failed");
        sunrealtype* yd = N_VGetArrayPointer(y_);
        for (int i = 0; i < nS_; ++i) yd[i] = Y0[i];
        yd[nS_] = T0;

        cvode_mem_ = CVodeCreate(CV_BDF, sunctx_);
        if (!cvode_mem_)
            throw std::runtime_error("CVODEWrapper: CVodeCreate failed");

        if (CVodeInit(cvode_mem_, &CVODEWrapper::rhsFn, t0, y_))
            throw std::runtime_error("CVODEWrapper: CVodeInit failed");

        CVodeSStolerances(cvode_mem_, rtol_, atol_);
        CVodeSetUserData(cvode_mem_, this);

        A_  = SUNDenseMatrix(nEq_, nEq_, sunctx_);
        LS_ = SUNLinSol_Dense(y_, A_, sunctx_);
        if (CVodeSetLinearSolver(cvode_mem_, LS_, A_))
            throw std::runtime_error("CVODEWrapper: CVodeSetLinearSolver failed");

        CVodeSetMaxNumSteps(cvode_mem_, 1000000);
        initialized_ = true;
    }

    // Advance the (persisted) integrator to tout; unpack into Y, T.
    // Returns false on a CVODE failure flag.
    bool advanceTo(double tout, std::vector<double>& Y, double& T)
    {
        if (!initialized_)
            throw std::runtime_error("CVODEWrapper: advanceTo before init");

        sunrealtype tret;
        int flag = CVode(cvode_mem_, tout, y_, &tret, CV_NORMAL);
        if (flag < 0) return false;

        const sunrealtype* yd = N_VGetArrayPointer(y_);
        Y.resize(nS_);
        for (int i = 0; i < nS_; ++i) Y[i] = yd[i];
        T = yd[nS_];
        return true;
    }

    // Diagnostics (adaptive-stepping cost)
    long numSteps() const
    { long n = 0; if (cvode_mem_) CVodeGetNumSteps(cvode_mem_, &n); return n; }

    long numRhsEvals() const
    { long n = 0; if (cvode_mem_) CVodeGetNumRhsEvals(cvode_mem_, &n); return n; }

    long numJacEvals() const
    { long n = 0; if (cvode_mem_) CVodeGetNumJacEvals(cvode_mem_, &n); return n; }
};

} // End namespace Foam

#endif
