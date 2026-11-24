#ifndef TurkelWeissPreconditioner_H
#define TurkelWeissPreconditioner_H

#include "RealFluidState.H"
#include "ConservedVars.H"
#include <vector>
#include <cmath>
#include <algorithm>

namespace Foam
{

/**
 * @brief Turkel-Weiss low-Mach preconditioning for real-fluid Euler equations
 * 
 * This class implements the local preconditioning matrix Γ that transforms
 * the pseudo-time derivative to accelerate convergence at low Mach numbers.
 * 
 * Key principle:
 *   Original system:  ∂U/∂t + ∇·F(U) = 0
 *   Preconditioned:   Γ·∂U/∂τ + ∇·F(U) = 0
 * 
 * where τ is pseudo-time and Γ scales the acoustic modes relative to 
 * convective modes, removing the O(1/M²) condition number stiffness.
 * 
 * Reference: Turkel, E. "Preconditioned methods for solving the incompressible
 * and low speed compressible equations." Journal of Computational Physics (1987)
 * 
 * @note This transforms ONLY the pseudo-time derivative. The flux Jacobian,
 * spatial discretization, and physical sound speed remain unchanged.
 */
class TurkelWeissPreconditioner
{
private:
    // Configuration
    double machCutoff_;     ///< Minimum reference Mach (prevents singularity at M→0)
    bool enabled_;          ///< Runtime enable/disable flag
    
    // Cached local state
    double localMach_;      ///< Current cell Mach number
    double beta_;           ///< Preconditioning parameter: β = max(M, M_cutoff)
    double betaSq_;         ///< β²
    double soundSpeed_;     ///< Physical sound speed (m/s)
    
    static constexpr double DEFAULT_MACH_CUTOFF = 0.05;
    
public:
    TurkelWeissPreconditioner()
        : machCutoff_(DEFAULT_MACH_CUTOFF),
          enabled_(true),
          localMach_(0.0),
          beta_(machCutoff_),
          betaSq_(machCutoff_ * machCutoff_),
          soundSpeed_(1.0)
    {}
    
    // Configuration
    void setEnabled(bool enable) { enabled_ = enable; }
    bool isEnabled() const { return enabled_; }
    
    void setMachCutoff(double cutoff) 
    { 
        machCutoff_ = std::max(cutoff, 1e-3);  // Safety floor
        if (beta_ < machCutoff_) {
            beta_ = machCutoff_;
            betaSq_ = beta_ * beta_;
        }
    }
    
    double getMachCutoff() const { return machCutoff_; }
    
    /**
     * @brief Update local preconditioning parameters for a cell
     * 
     * @param velocity Local velocity vector (m/s)
     * @param soundSpeed Physical equilibrium sound speed (m/s)
     */
    void updateLocal(const std::vector<double>& velocity, double soundSpeed)
    {
        soundSpeed_ = std::max(soundSpeed, 1.0);  // Safety floor
        
        double velMag = 0.0;
        for (const auto& v : velocity)
            velMag += v * v;
        velMag = std::sqrt(velMag);
        
        localMach_ = velMag / soundSpeed_;
        beta_ = std::max(localMach_, machCutoff_);
        betaSq_ = beta_ * beta_;
    }

    // Accessors
    double localMach() const { return localMach_; }
    double beta() const { return beta_; }
    double betaSquared() const { return betaSq_; }
    
    /**
     * @brief Apply preconditioning matrix to pseudo-time derivative: dU_precond = Γ·dU
     * 
     * The preconditioning matrix Γ for conserved variables U = (ρ, ρu, ρv, ρw, ρE, ρY_i, ρZ, ρC):
     * 
     *   Γ = diag(1, β²/a², β²/a², β²/a², 1, 1, ..., 1)
     * 
     * This scales the momentum equations by β²/a², which at low Mach numbers (β ≈ M_cutoff << 1)
     * dramatically reduces the pseudo-time stiffness between acoustic and convective modes.
     * 
     * @param dU Input: unpreconditioned pseudo-time derivative
     * @return Output: preconditioned derivative Γ·dU
     */
    ConservedVars applyToPseudoTimeDerivative(const ConservedVars& dU) const
    {
        if (!enabled_)
            return dU;
        
        ConservedVars result = dU;
        
        // Density: unchanged
        // result.rho = dU.rho;
        
        // Momentum: scale by β²/a²
        const double momentumScale = betaSq_ / (soundSpeed_ * soundSpeed_);
        for (std::size_t i = 0; i < result.rhoU.size(); ++i)
            result.rhoU[i] *= momentumScale;
        
        // Energy: unchanged
        // result.rhoE = dU.rhoE;
        
        // Species, progress variables: unchanged
        // (rhoY, rhoZ, rhoC remain as-is)
        
        return result;
    }
    
    /**
     * @brief Compute the effective pseudo-timestep scale factor
     * 
     * For stability, the local pseudo-timestep should be scaled by the 
     * preconditioning parameter:
     * 
     *   dt_eff = β · dt_acoustic
     * 
     * This prevents the preconditioned system from violating CFL constraints
     * on the (now dominant) convective waves.
     * 
     * @return Scale factor to multiply dt by (dimensionless)
     */
    double pseudoTimestepScale() const
    {
        return enabled_ ? beta_ : 1.0;
    }
    
    /**
     * @brief Inverse preconditioning: dU = Γ⁻¹·(Γ·dU)
     * 
     * Used in implicit solvers where we need to recover the physical update
     * from the preconditioned linear system solution.
     */
    ConservedVars applyInverse(const ConservedVars& gammaDU) const
    {
        if (!enabled_)
            return gammaDU;
        
        ConservedVars result = gammaDU;
        
        const double invMomentumScale = (soundSpeed_ * soundSpeed_) / betaSq_;
        for (std::size_t i = 0; i < result.rhoU.size(); ++i)
            result.rhoU[i] *= invMomentumScale;
        
        return result;
    }
    
    /**
     * @brief Capability labels for runtime introspection
     */
    static constexpr bool modifiesFluxResidual() { return false; }
    static constexpr bool modifiesPseudoTimeDerivative() { return true; }
    static constexpr bool modifiesImplicitJacobian() { return false; }  // Not yet implemented
    
    static constexpr const char* capabilityLabel()
    {
        return "turkel-weiss-low-mach-preconditioning";
    }
};

} // End namespace Foam

#endif
