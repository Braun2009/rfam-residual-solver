#ifndef ZeroD_Reactor_H
#define ZeroD_Reactor_H

#include "ChemistryModel.H"
#include <vector>
#include <algorithm>

namespace Foam
{

class ZeroD_Reactor
{
private:
    double p0_, T_, rho_, Cv_;
    std::vector<double> Y_;
    ChemistryModel* chem_;
    double t_;
    bool ignited_;
    static constexpr double T_IGN = 1500.0;

public:
    ZeroD_Reactor(ChemistryModel* chem, double Cv = 1500.0)
        : chem_(chem), Cv_(Cv), t_(0.0), ignited_(false) {}

    void initialize(double p0, double T0, const std::vector<double>& Y0)
    {
        p0_ = p0; T_ = T0; Y_ = Y0; t_ = 0.0; ignited_ = false;
        // Ideal gas density estimate
        double MW_inv = 0.0;
        const auto& MW = chem_->molecularWeights();
        for (int i = 0; i < chem_->nSpecies(); i++)
            MW_inv += Y0[i] / MW[i];
        rho_ = p0 * MW_inv / (8314.0 * T0);
    }

    void step(double dt)
    {
        std::vector<double> omega_i;
        double qdot;
        chem_->computeSourceTerms(T_, p0_, rho_, Y_, omega_i, qdot);

        for (int i = 0; i < chem_->nSpecies(); i++)
        {
            Y_[i] += dt * omega_i[i] / rho_;
            if (Y_[i] < 0.0) Y_[i] = 0.0;
        }
        // Normalize
        double sum = 0.0;
        for (auto& y : Y_) sum += y;
        if (sum > 0.0) for (auto& y : Y_) y /= sum;

        T_ += dt * qdot / (rho_ * Cv_);
        t_ += dt;

        if (!ignited_ && T_ > T_IGN) ignited_ = true;
    }

    double runToIgnition(double t_max, double dt = 1e-8)
    {
        while (t_ < t_max)
        {
            step(dt);
            if (ignited_) return t_;
        }
        return -1.0;
    }

    double T() const { return T_; }
    double t() const { return t_; }
    bool ignited() const { return ignited_; }
    const std::vector<double>& Y() const { return Y_; }
};

} // End namespace Foam

#endif
