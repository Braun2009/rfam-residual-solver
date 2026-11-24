/*---------------------------------------------------------------------------*\
  realFluidAllMach - Species Critical Properties Database Implementation
\*---------------------------------------------------------------------------*/

#include "SpeciesCriticalDB.H"
#include <stdexcept>

namespace Foam
{

SpeciesCriticalDB::SpeciesCriticalDB()
{}

void SpeciesCriticalDB::add(
    const std::string& name,
    double Tc,
    double Pc,
    double omega,
    double MW
)
{
    CriticalProperties props;
    props.Tc = Tc;
    props.Pc = Pc;
    props.omega = omega;
    props.MW = MW;
    data_[name] = props;
}

const CriticalProperties& SpeciesCriticalDB::get(const std::string& name) const
{
    auto it = data_.find(name);
    if (it == data_.end())
    {
        throw std::runtime_error("Species not found in database: " + name);
    }
    return it->second;
}

bool SpeciesCriticalDB::has(const std::string& name) const
{
    return data_.find(name) != data_.end();
}

SpeciesCriticalDB SpeciesCriticalDB::createDefault()
{
    SpeciesCriticalDB db;

    // Common rocket propellants
    db.add("O2",    154.58,  5.043e6,  0.0222,   0.0319988);
    db.add("N2",    126.19,  3.396e6,  0.0372,   0.0280134);
    db.add("H2",    32.98,   1.293e6,  -0.2160,  0.00201588);
    db.add("CH4",   190.56,  4.599e6,  0.0115,   0.01604246);
    db.add("C2H6",  305.32,  4.872e6,  0.0995,   0.03006904);
    db.add("C3H8",  369.89,  4.251e6,  0.1523,   0.04409562);
    db.add("H2O",   647.10,  22.064e6, 0.3443,   0.01801528);
    db.add("CO",    132.86,  3.494e6,  0.0497,   0.0280101);
    db.add("CO2",   304.13,  7.377e6,  0.2236,   0.0440095);
    db.add("NH3",   405.40,  11.357e6, 0.2526,   0.017031);
    db.add("RP1",   658.0,   1.400e6,  0.5800,   0.170000);  // Approximate

    return db;
}

} // End namespace Foam
