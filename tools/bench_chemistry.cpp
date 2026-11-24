// Chemistry-only throughput benchmark.
//
// Answers one question before any mechanism reduction work: how much of a
// solver step is the reaction source term, for the 6-species two-step mechanism
// against the 59-species San Diego mechanism?
//
// Runs the same number of cell evaluations for both, on a representative
// chamber state, and reports microseconds per cell evaluation and the implied
// seconds per solver step for a grid of NC cells at S chemistry sub-steps.
#include "TwoStepCO.H"
#include "SandiegoChemistry.H"
#include "SandiegoReducedChemistry.H"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

namespace
{
struct Case
{
    const char* name;
    Foam::ChemistryModel* model;
};

double run(Foam::ChemistryModel& model,std::size_t ncells,int substeps,
           double& qdotOut)
{
    const int ns=model.nSpecies();
    std::vector<double> Y(static_cast<std::size_t>(ns),0.0);
    std::vector<double> omega(static_cast<std::size_t>(ns),0.0);
    // Burnt-gas-like state: hot, compressed, single phase.
    const double T=3400.0,p=4.0e6,rho=4.6;
    // Distribute over the species the mechanism actually has, weighting the
    // products so the rate evaluation sees a chemically plausible mixture.
    double sum=0.0;
    for(int s=0;s<ns;++s)
    {
        const std::string& n=model.speciesNames()[static_cast<std::size_t>(s)];
        double w=0.0;
        if(n=="CO2")w=0.55;
        else if(n=="H2O")w=0.30;
        else if(n=="O2")w=0.06;
        else if(n=="C3H8")w=0.02;
        else if(n=="N2")w=0.07;
        else w=1e-6;
        Y[static_cast<std::size_t>(s)]=w;
        sum+=w;
    }
    for(double& y:Y) y/=sum;

    double qdot=0.0;
    const auto start=std::chrono::steady_clock::now();
    for(int step=0;step<substeps;++step)
        for(std::size_t c=0;c<ncells;++c)
            model.computeSourceTerms(T,p,rho,Y,omega,qdot);
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    qdotOut=qdot;
    return seconds;
}
}

int main(int argc,char** argv)
{
    // The A grid has 1.76M cells; keep the default workload small enough that
    // the benchmark itself stays inside a few minutes per mechanism.
    std::size_t ncells=1761002;
    int substeps=2;
    for(int i=1;i<argc;++i)
    {
        const std::string a=argv[i];
        if(a=="--cells"&&i+1<argc)ncells=static_cast<std::size_t>(std::atof(argv[++i]));
        else if(a=="--substeps"&&i+1<argc)substeps=std::atoi(argv[++i]);
    }
    std::printf("cells=%zu substeps=%d\n",ncells,substeps);

    Foam::TwoStepCO twoStep;
    Foam::SandiegoChemistry sandiego;
    Foam::SandiegoReducedChemistry reduced;
    Case cases[3]={{ "two-step", &twoStep },{ "sandiego", &sandiego },
                   { "reduced", &reduced }};

    for(const Case& c:cases)
    {
        double qdot=0.0;
        const double seconds=run(*c.model,ncells,substeps,qdot);
        const double perCellUs=seconds/double(ncells)/double(substeps)*1e6;
        std::printf("%-9s species=%2d reactions=%3d  total=%8.2f s  "
                    "per_cell=%6.2f us  qdot=%.3e W/m3\n",
                    c.name,c.model->nSpecies(),c.model->nReactions(),
                    seconds,perCellUs,qdot);
    }
    return 0;
}
