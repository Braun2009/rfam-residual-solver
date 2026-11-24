// Quick compilation test for TurkelWeissPreconditioner
#include "src/numerics/preconditioning/TurkelWeissPreconditioner.H"
#include "src/variables/ConservedVars.H"
#include <iostream>

int main() {
    using namespace Foam;
    
    TurkelWeissPreconditioner precond;
    precond.setMachCutoff(0.05);
    precond.setEnabled(true);
    
    std::vector<double> velocity = {10.0, 5.0, 0.0};  // 10 m/s in x, 5 m/s in y
    double soundSpeed = 540.0;  // LOX
    
    precond.updateLocal(velocity, soundSpeed);
    
    std::cout << "Local Mach: " << precond.localMach() << std::endl;
    std::cout << "Beta: " << precond.beta() << std::endl;
    std::cout << "Pseudo-timestep scale: " << precond.pseudoTimestepScale() << std::endl;
    
    // Test preconditioning a dummy conserved variable
    ConservedVars dU;
    dU.rho = 1.0;
    dU.rhoU = {100.0, 50.0, 0.0};
    dU.rhoE = 1000.0;
    
    ConservedVars precondDU = precond.applyToPseudoTimeDerivative(dU);
    
    std::cout << "Original momentum[0]: " << dU.rhoU[0] << std::endl;
    std::cout << "Preconditioned momentum[0]: " << precondDU.rhoU[0] << std::endl;
    std::cout << "Scale factor: " << precondDU.rhoU[0] / dU.rhoU[0] << std::endl;
    
    return 0;
}
