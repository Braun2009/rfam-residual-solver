#include <gtest/gtest.h>
#include "TurkelWeissPreconditioner.H"
#include "ConservedVars.H"
#include <cmath>

using namespace Foam;

TEST(TurkelWeissPreconditioner, BasicConstruction) {
    TurkelWeissPreconditioner precond;
    EXPECT_TRUE(precond.isEnabled());
    EXPECT_NEAR(precond.getMachCutoff(), 0.05, 1e-10);
}

TEST(TurkelWeissPreconditioner, MachNumberUpdate) {
    TurkelWeissPreconditioner precond;
    precond.setMachCutoff(0.05);
    
    // High Mach case: M = 0.5
    std::vector<double> velocity = {270.0, 0.0, 0.0};
    double soundSpeed = 540.0;
    precond.updateLocal(velocity, soundSpeed);
    
    EXPECT_NEAR(precond.localMach(), 0.5, 1e-6);
    EXPECT_NEAR(precond.beta(), 0.5, 1e-6);
    EXPECT_NEAR(precond.pseudoTimestepScale(), 0.5, 1e-6);
}

TEST(TurkelWeissPreconditioner, LowMachCutoff) {
    TurkelWeissPreconditioner precond;
    precond.setMachCutoff(0.05);
    
    // Very low Mach: M = 0.02, should be clamped to cutoff
    std::vector<double> velocity = {10.8, 0.0, 0.0};  // ~10.8 m/s
    double soundSpeed = 540.0;
    precond.updateLocal(velocity, soundSpeed);
    
    EXPECT_NEAR(precond.localMach(), 0.02, 1e-6);
    EXPECT_NEAR(precond.beta(), 0.05, 1e-6);  // Clamped to cutoff
    EXPECT_NEAR(precond.pseudoTimestepScale(), 0.05, 1e-6);
}

TEST(TurkelWeissPreconditioner, MomentumScaling) {
    TurkelWeissPreconditioner precond;
    precond.setMachCutoff(0.05);
    
    std::vector<double> velocity = {27.0, 0.0, 0.0};  // M = 0.05
    double soundSpeed = 540.0;
    precond.updateLocal(velocity, soundSpeed);
    
    ConservedVars dU;
    dU.rho = 1.0;
    dU.rhoU = {100.0, 50.0, 25.0};
    dU.rhoE = 1000.0;
    dU.rhoY = {};
    dU.rhoZ = 0.0;
    dU.rhoC = 0.0;
    
    ConservedVars gammaDU = precond.applyToPseudoTimeDerivative(dU);
    
    // Density and energy unchanged
    EXPECT_DOUBLE_EQ(gammaDU.rho, dU.rho);
    EXPECT_DOUBLE_EQ(gammaDU.rhoE, dU.rhoE);
    
    // Momentum scaled by β²/a²
    double expectedScale = (0.05 * 0.05) / (540.0 * 540.0);
    EXPECT_NEAR(gammaDU.rhoU[0], dU.rhoU[0] * expectedScale, 1e-10);
    EXPECT_NEAR(gammaDU.rhoU[1], dU.rhoU[1] * expectedScale, 1e-10);
    EXPECT_NEAR(gammaDU.rhoU[2], dU.rhoU[2] * expectedScale, 1e-10);
}

TEST(TurkelWeissPreconditioner, InverseOperation) {
    TurkelWeissPreconditioner precond;
    precond.setMachCutoff(0.05);
    
    std::vector<double> velocity = {54.0, 0.0, 0.0};  // M = 0.1
    double soundSpeed = 540.0;
    precond.updateLocal(velocity, soundSpeed);
    
    ConservedVars dU;
    dU.rho = 1141.0;
    dU.rhoU = {100.0, 50.0, 25.0};
    dU.rhoE = 50000.0;
    dU.rhoY = {};
    dU.rhoZ = 0.0;
    dU.rhoC = 0.0;
    
    ConservedVars gammaDU = precond.applyToPseudoTimeDerivative(dU);
    ConservedVars recovered = precond.applyInverse(gammaDU);
    
    // Should recover original (within floating-point error)
    EXPECT_NEAR(recovered.rho, dU.rho, 1e-6);
    EXPECT_NEAR(recovered.rhoU[0], dU.rhoU[0], 1e-6);
    EXPECT_NEAR(recovered.rhoU[1], dU.rhoU[1], 1e-6);
    EXPECT_NEAR(recovered.rhoU[2], dU.rhoU[2], 1e-6);
    EXPECT_NEAR(recovered.rhoE, dU.rhoE, 1e-6);
}

TEST(TurkelWeissPreconditioner, DisabledPassthrough) {
    TurkelWeissPreconditioner precond;
    precond.setEnabled(false);
    
    std::vector<double> velocity = {27.0, 0.0, 0.0};
    double soundSpeed = 540.0;
    precond.updateLocal(velocity, soundSpeed);
    
    ConservedVars dU;
    dU.rho = 1.0;
    dU.rhoU = {100.0, 50.0, 25.0};
    dU.rhoE = 1000.0;
    dU.rhoY = {};
    dU.rhoZ = 0.0;
    dU.rhoC = 0.0;
    
    ConservedVars gammaDU = precond.applyToPseudoTimeDerivative(dU);
    
    // Should be unchanged
    EXPECT_DOUBLE_EQ(gammaDU.rho, dU.rho);
    EXPECT_DOUBLE_EQ(gammaDU.rhoU[0], dU.rhoU[0]);
    EXPECT_DOUBLE_EQ(gammaDU.rhoU[1], dU.rhoU[1]);
    EXPECT_DOUBLE_EQ(gammaDU.rhoU[2], dU.rhoU[2]);
    EXPECT_DOUBLE_EQ(gammaDU.rhoE, dU.rhoE);
}
