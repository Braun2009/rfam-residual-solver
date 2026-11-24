/*---------------------------------------------------------------------------*\
  realFluidAllMach - Species Critical Properties Database
\*---------------------------------------------------------------------------*/

#ifndef SpeciesCriticalDB_H
#define SpeciesCriticalDB_H

#include <map>
#include <string>

namespace Foam
{

struct CriticalProperties
{
    double Tc;      // Critical temperature [K]
    double Pc;      // Critical pressure [Pa]
    double omega;   // Acentric factor [-]
    double MW;      // Molecular weight [kg/mol]
};

class SpeciesCriticalDB
{
private:
    std::map<std::string, CriticalProperties> data_;

public:
    SpeciesCriticalDB();

    void add(const std::string& name, double Tc, double Pc, double omega, double MW);

    const CriticalProperties& get(const std::string& name) const;

    bool has(const std::string& name) const;

    // Common species
    static SpeciesCriticalDB createDefault();
};

} // End namespace Foam

#endif
