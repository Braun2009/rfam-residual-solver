#ifndef TetraSteadyMonitor_H
#define TetraSteadyMonitor_H

#include "TetraFiniteVolume.H"
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Foam
{
// Fixed engineering scales: mdot, mdot*a_ref, mdot*h_ref and species mdot.
// Local integrated |R| prevents cancellations from masquerading as steady
// flow. Boundary budgets are an independent check of those local residuals.
// Numerical convergence is not a geometry/transport/production certification.
class TetraSteadyMonitor
{
public:
    struct Report
    {
        std::array<double,7> local{}; // mass, mx,my,mz,energy,O2,C3H8
        std::array<double,7> netBoundary{};
        double maxLocal=0.,maxBudget=0.,change=0.;
        int stableSamples=0;
        bool converged=false;
    };
private:
    std::array<double,7> scales_;
    int stable_=0,lastIteration_=-1;
public:
    TetraSteadyMonitor(double massFlow,double sound,double energy,
                       double oxygenFlow,double fuelFlow)
    {
        if(!std::isfinite(massFlow)||massFlow<=0.||!std::isfinite(sound)||sound<=0.
           ||!std::isfinite(energy)||energy<=0.
           ||!std::isfinite(oxygenFlow)||oxygenFlow<0.
           ||!std::isfinite(fuelFlow)||fuelFlow<0.)
            throw std::invalid_argument("invalid steady monitor scale");
        scales_={massFlow,massFlow*sound,massFlow*sound,massFlow*sound,massFlow*energy,
                 std::max(oxygenFlow,1e-6*massFlow),std::max(fuelFlow,1e-6*massFlow)};
    }
    static std::array<double,7> values(const ConservedVars& u)
    {
        // Species-general (2026-10-06): the residual budget tracks the first
        // two species densities in the caller's own order (cold flow: O2,C3H8;
        // reacting: C3H8,O2). The caller passes the matching reference flows.
        if(u.rhoU.size()!=3||u.rhoY.size()<2)
            throw std::runtime_error("steady residual requires at least the fuel/oxidizer pair");
        return {u.rho,u.rhoU[0],u.rhoU[1],u.rhoU[2],u.rhoE,u.rhoY[0],u.rhoY[1]};
    }
    Report evaluate(const TetraFiniteVolume& fv,const std::vector<ConservedVars>& residual,
                    double change,std::size_t limitedCells,int iteration)
    {
        if(residual.size()!=fv.cellCount())throw std::invalid_argument("steady residual cell count mismatch");
        if(iteration<=lastIteration_)throw std::invalid_argument("steady monitor requires increasing iterations");
        lastIteration_=iteration;
        Report report;report.change=change;
        for(const auto& u:residual)
        {
            const auto r=values(u);
            for(std::size_t j=0;j<r.size();++j)
            {
                if(!std::isfinite(r[j]))throw std::runtime_error("nonfinite steady residual");
                report.local[j]+=std::abs(r[j])/scales_[j];
            }
        }
        for(std::size_t f=0;f<fv.mesh().faceCount();++f)
        {
            if(fv.mesh().face(f).neighbour>=0)continue;
            const auto flux=values(fv.lastFaceFlux(f));
            for(std::size_t j=0;j<flux.size();++j)
            {
                if(!std::isfinite(flux[j]))throw std::runtime_error("nonfinite steady boundary budget");
                report.netBoundary[j]+=fv.faceArea(f)*flux[j]/scales_[j];
            }
        }
        for(std::size_t j=0;j<scales_.size();++j)
        {
            report.maxLocal=std::max(report.maxLocal,report.local[j]);
            report.maxBudget=std::max(report.maxBudget,std::abs(report.netBoundary[j]));
        }
        const bool passed=iteration>0&&std::isfinite(change)&&change>=0.&&change<=1e-5
            &&report.maxLocal<=1e-5&&report.maxBudget<=1e-3&&limitedCells==0;
        stable_=passed?stable_+1:0;report.stableSamples=stable_;
        report.converged=stable_>=5;return report;
    }
};
}
#endif
