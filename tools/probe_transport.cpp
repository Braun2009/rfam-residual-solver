// Probe CoolProp transport properties for the injector states.
#include "CoolProp.h"
#include "AbstractState.h"
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace CoolProp;
static void show(const char* name, const std::string& fluid, const std::vector<double>& y,
                 double p, double T)
{
    try
    {
        auto* s = AbstractState::factory("HEOS", fluid);
        if (y.size() == 2) s->set_mass_fractions(y);
        s->update(PT_INPUTS, p, T);
        std::printf("%-28s p=%.4gMPa T=%.4gK rho=%12.6g mu=%12.6g lambda=%12.6g cp=%10.6g phase=%d\n",
                    name, p/1e6, T, s->rhomass(), s->viscosity(), s->conductivity(),
                    s->cpmass(), int(s->phase()));
        delete s;
    }
    catch (const std::exception& ex)
    {
        std::printf("%-28s FAIL %s\n", name, ex.what());
    }
}
int main()
{
    show("LOX pure @130K/6MPa", "Oxygen", {}, 6e6, 130.);
    show("C3H8 pure @450K/4.5MPa", "Propane", {}, 4.5e6, 450.);
    show("mix Y=.9295 @250K/4MPa", "Oxygen&Propane", {0.92953, 0.07047}, 4e6, 250.);
    show("mix Y=.8728 @250K/4MPa", "Oxygen&Propane", {0.87283, 0.12717}, 4e6, 250.);
    show("mix Y=.7273 @235K/4.1MPa", "Oxygen&Propane", {0.72727, 0.27273}, 4.11e6, 235.);
    show("mix Y=.7273 @300K/4.1MPa", "Oxygen&Propane", {0.72727, 0.27273}, 4.11e6, 300.);
    show("mix Y=.9295 @165K/3.8MPa", "Oxygen&Propane", {0.92953, 0.07047}, 3.8e6, 165.);
    return 0;
}
