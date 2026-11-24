#ifndef TetraFiniteVolume_H
#define TetraFiniteVolume_H

#include "TetraMesh.H"
#include "TetraViscousFlux.H"
#include "FluxScheme.H"
#include "EquilibriumThermoModel.H"
#include "TetraFaceInvariant.H"
#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Foam
{

// Boundary ghost states are supplied per boundary face by the driver. The
// kernel never invents injector physics; it only requires a complete primitive
// and equilibrium pair for the ghost cell.
using TetraBoundaryProvider = std::function<void(
    std::size_t cell,std::size_t face,
    PrimitiveVars& ghostW,thermo::EquilibriumState& ghostQ)>;

// Optional physical boundary flux, per unit face area. Return false to use
// the numerical ghost-state flux. Useful for prescribed conservative inlets.
using TetraBoundaryFluxProvider = std::function<bool(
    std::size_t cell,std::size_t face,ConservedVars& flux)>;

// First failure inside the parallel residual loop, reported instead of
// letting an exception escape an OpenMP region (which would terminate).
struct TetraKernelFailure
{
    bool failed=false;
    std::size_t cell=0;
    std::string detail;
};

// Unstructured cell-centred finite-volume kernel for linear tetrahedra.
//
// Evaluate each shared face exactly once, then gather it with opposite signs
// into its owner and neighbour. No state-admission fallback is permitted.
//
// WAVE SPEED OWNERSHIP. Every spectral-radius evaluation in this class -- and in
// TetraLUSGS, which consumes the same quantity -- goes through
// tetraFaceInvariant() in TetraFaceInvariant.H. Do not write `|u_n| + a` inline
// anywhere in the tetra solver; that duplication is what let the LU-SGS bound
// and the live flux branch drift apart. See that header for the full rationale.
class TetraFiniteVolume
{
    const TetraMesh& mesh_;
    std::vector<PrimitiveVars> W_;
    std::vector<ConservedVars> U_;
    std::vector<thermo::EquilibriumState> Q_;
    std::unique_ptr<FluxScheme> flux_;
    mutable std::function<void(std::size_t, std::size_t)> progressCallback_;
public:
    // Optional progress sink for the face loop. The kernel never includes the
    // telemetry header; the driver passes a lambda. Called roughly every 1024
    // faces, i.e. ~3.7k times for the 3.8M-face A grid: measurable nowhere,
    // but enough to tell "slow" from "stuck".
    void setProgressCallback(std::function<void(std::size_t, std::size_t)> callback);

    // Forwarded to the AUSM flux when it supports the interface-consistency
    // hook (see AUSMPlusUpFlux::setConsistentDensityHook).
    void setConsistentDensityHook(std::function<bool(double,double,double,double&)> hook);
    std::size_t consistentDensityFallbacks() const;
private:
    std::vector<std::string> species_;
    std::vector<double> areas_;
    std::vector<std::vector<double>> normals_;
    mutable std::vector<ConservedVars> faceFlux_;
    TetraBoundaryFluxProvider boundaryFlux_;
    std::function<std::pair<double,double>(const PrimitiveVars&,
                                          const thermo::EquilibriumState&)> molecularTransport_;
    // Optional turbulent addition (mu_t, lambda_t) per cell. Kept separate from
    // molecularTransport_ because the provider needs the cell index and the
    // eddy viscosity is a transported field, not a pure state function.
    std::function<std::pair<double,double>(std::size_t,const PrimitiveVars&,
                                          const thermo::EquilibriumState&)> turbulentTransport_;
    mutable std::vector<double> molecularSpectral_;
    // Last molecular (mu, lambda) arrays actually used by the residual. The
    // turbulence operator needs mu to build the omega wall value and its own
    // stability limit, and recomputing CoolProp there would double the cost.
    mutable std::vector<double> lastMu_, lastLambda_;

public:
    TetraFiniteVolume(const TetraMesh& mesh,std::vector<std::string> species,
                      std::unique_ptr<FluxScheme> flux);

    void setState(std::size_t cell,const PrimitiveVars& w,const ConservedVars& u,
                  const thermo::EquilibriumState& q);
    void setPrimitive(std::size_t cell,const PrimitiveVars& w,const thermo::EquilibriumState& q);
    // Transaction exchange: caller retains the old state without EOS recovery.
    void exchangeState(std::size_t cell,PrimitiveVars& w,ConservedVars& u,
                       thermo::EquilibriumState& q);

    const PrimitiveVars& W(std::size_t i) const { return W_.at(i); }
    const ConservedVars& U(std::size_t i) const { return U_.at(i); }
    const thermo::EquilibriumState& Q(std::size_t i) const { return Q_.at(i); }
    const TetraMesh& mesh() const { return mesh_; }
    const std::vector<std::string>& species() const { return species_; }

    std::size_t cellCount() const { return mesh_.cellCount(); }
    // Retained for old output compatibility; strict kernel never falls back.
    mutable std::atomic<std::size_t> boundaryFluxFallbacks{0};
    double faceArea(std::size_t f) const { return areas_.at(f); }
    const std::vector<double>& faceNormal(std::size_t f) const { return normals_.at(f); }
    const ConservedVars& lastFaceFlux(std::size_t f) const { return faceFlux_.at(f); }
    void setBoundaryFluxProvider(TetraBoundaryFluxProvider provider)
    { boundaryFlux_=std::move(provider); }
    // Opt-in molecular Navier-Stokes residual; the provider must fail closed
    // for any state whose transport closure has not been validated.
    void setMolecularTransportProvider(std::function<std::pair<double,double>(
        const PrimitiveVars&,const thermo::EquilibriumState&)> provider)
    { molecularTransport_=std::move(provider); molecularSpectral_.clear(); }
    double molecularSpectral(std::size_t cell) const
    { return molecularSpectral_.empty()?0.:molecularSpectral_.at(cell); }
    // Adds (mu_t, lambda_t) to the molecular transport of every cell so the
    // eddy viscosity enters the viscous flux AND the LU-SGS diagonal through
    // the same mu/lambda arrays.
    void setTurbulentTransportProvider(std::function<std::pair<double,double>(
        std::size_t,const PrimitiveVars&,const thermo::EquilibriumState&)> provider)
    { turbulentTransport_=std::move(provider); }
    bool hasTurbulentTransport() const { return static_cast<bool>(turbulentTransport_); }
    const std::vector<double>& lastMu() const { return lastMu_; }
    const std::vector<double>& lastLambda() const { return lastLambda_; }

    // Residual R_c = sum_faces A_outward . F  [units: per second, already
    // multiplied by face area]. Divide by volume to get dU/dt.
    void computeResidual(const TetraBoundaryProvider& boundary,
                         std::vector<ConservedVars>& residual,
                         TetraKernelFailure* failure=nullptr) const;

    // Largest local explicit step CFL*V/sum((|u|+a)|A|) over all cells; the
    // per-cell limit field is what a local-time-stepping solver should use.
    void localTimeStep(double cfl,std::vector<double>& dt) const;
    // Boundary-aware physical stability step. Boundary traces participate in
    // the spectral bound; they may be much faster than the adjacent cell.
    void localTimeStep(double cfl,std::vector<double>& dt,
                       const TetraBoundaryProvider& boundary) const;

    // Mass flux integrated over a face, including a prescribed inlet flux.
    double faceMassFlow(std::size_t face,const PrimitiveVars& inside,
                        const PrimitiveVars& outside,
                        const thermo::EquilibriumState& insideQ,
                        const thermo::EquilibriumState& outsideQ) const;
};

} // namespace Foam
#endif
