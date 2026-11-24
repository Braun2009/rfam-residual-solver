#ifndef PengRobinsonEOS_H
#define PengRobinsonEOS_H

#include <limits>
#include <vector>
#include "SpeciesCriticalDB.H"
#include "BinaryInteractionDB.H"
#include "VdWOneFluidMixing.H"

namespace Foam
{

class PengRobinsonEOS
{
public:
    struct RootResult
    {
        double Z = std::numeric_limits<double>::quiet_NaN();
        double G = std::numeric_limits<double>::quiet_NaN();
        bool valid = false;
        int numRoots = 0;
    };

    // GAS and LIQUID select the largest and smallest physical cubic roots.
    // AUTO selects the physical root with minimum residual Gibbs energy.
    enum class PhaseHint { GAS, LIQUID, AUTO };

private:
    const SpeciesCriticalDB& criticalData_;
    const BinaryInteractionDB& kij_;
    VdWOneFluidMixing mixingRule_;

    static constexpr double Tr_lower = 0.95;
    static constexpr double Tr_upper = 1.05;
    static constexpr double Pr_lower = 0.9;
    static constexpr double Pr_upper = 1.5;

    double computeGibbs(double Z, double A, double B, double T) const;
    double selectRootGibbs(double r1, double r2, double r3, int numRoots, double A, double B, double T, PhaseHint hint) const;

public:
    PengRobinsonEOS(const SpeciesCriticalDB& criticalData, const BinaryInteractionDB& kij);

    RootResult solveZ(double p, double T, const std::vector<double>& Y, PhaseHint hint = PhaseHint::AUTO) const;
    double density(double p, double T, const std::vector<double>& Y) const;
    double pressure(double rho, double T, const std::vector<double>& Y) const;
    double enthalpyDeparture(double p, double T, double Z, const std::vector<double>& Y) const;
    double entropyDeparture(double p, double T, double Z, const std::vector<double>& Y) const;
    double cvDeparture(double p, double T, double Z, const std::vector<double>& Y) const;
    bool isNearCritical(double p, double T, const std::vector<double>& Y) const;
    double pseudoCriticalTemperature(const std::vector<double>& Y) const;
    double pseudoCriticalPressure(const std::vector<double>& Y) const;

    // Mixture molecular weight [kg/mol] via the van der Waals one-fluid rule.
    double MW(const std::vector<double>& Y) const;

    // Phase 7 Step 2: route the caller's species order (global SpeciesRegistry
    // order for reacting flow) down to the mixing rule. Throws on an unknown /
    // non-physical species (no silent fallback). Opt-in: default construction
    // stays {O2,C3H8} so legacy EOS results are unchanged.
    void setSpeciesNames(const std::vector<std::string>& names)
    {
        mixingRule_.setSpeciesNames(names);
        // Enable true quadratic VdW mixing ONLY on the external opt-in path.
        // The internal ctor call keeps multiSpeciesEnabled_=false so that the
        // default {O2,C3H8} EOS results stay byte-for-byte unchanged.
        mixingRule_.enableMultiSpeciesMixing();
    }

    const SpeciesCriticalDB& criticalData() const { return criticalData_; }
};

} // End namespace Foam

#endif
