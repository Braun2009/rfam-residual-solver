/* Phase 7: compile/smoke for ChemistrySourceOperator (framework wrapper).
   Verifies it builds against framework Operator + FiniteRateReactor + SUNDIALS.
   Behavioral validation of the chemistry itself = test_phase7_chem0d.
   Full 1D reacting-flow (apply over a Domain + EOS reback) = Step 2. */
#include <iostream>
#include <memory>
#include "src/thermo/SpeciesRegistry.H"
#include "src/combustion/ChemistrySourceOperator.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"

using namespace Foam;

int main()
{
    SpeciesRegistry reg;
    TwoStepCO chem;
    std::unique_ptr<Operator> op(new ChemistrySourceOperator(&chem, reg));
    std::cout << "operator: " << op->name() << "\n";
    std::cout << "PHASE7_OPERATOR_COMPILE_OK\n";
    return 0;
}
