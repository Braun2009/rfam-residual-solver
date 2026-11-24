#ifndef ArrheniusRate_H
#define ArrheniusRate_H

#include <cmath>

namespace Foam
{

/*---------------------------------------------------------------------------*\
    ArrheniusRate - Arrhenius reaction rate calculator
    
    k(T) = A * T^n * exp(-Ea / (R*T))
    
    where:
        A  = pre-exponential factor
        n  = temperature exponent
        Ea = activation energy [J/mol]
        R  = universal gas constant = 8.314 J/(mol·K)
\*---------------------------------------------------------------------------*/

class ArrheniusRate
{
private:
    double A_;     // Pre-exponential factor [units depend on reaction order]
    double n_;     // Temperature exponent [-]
    double Ea_;    // Activation energy [J/mol]
    
    // Keep the kinetic rate constant on the same CODATA basis as NASA7Thermo,
    // Cantera, and the equilibrium-constant helper.  The old 8.314 truncation
    // introduced a temperature-dependent implementation mismatch.
    static constexpr double R_UNIVERSAL = 8.314462618;  // J/(mol·K)
    static constexpr double SMALL_T = 1.0;        // Avoid division by zero

public:
    // Constructor
    ArrheniusRate(double A = 0.0, double n = 0.0, double Ea = 0.0)
        : A_(A), n_(n), Ea_(Ea)
    {
    }
    
    // Compute rate constant at given temperature
    double rate(double T) const
    {
        if (T < SMALL_T) return 0.0;
        
        double k = A_ * std::pow(T, n_) * std::exp(-Ea_ / (R_UNIVERSAL * T));
        
        return k;
    }
    
    // Getters
    double A() const { return A_; }
    double n() const { return n_; }
    double Ea() const { return Ea_; }
    
    // Setters (for parameter tuning)
    void setA(double A) { A_ = A; }
    void setN(double n) { n_ = n; }
    void setEa(double Ea) { Ea_ = Ea; }
};

} // End namespace Foam

#endif
