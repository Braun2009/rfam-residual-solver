#ifndef StructuredImplicitBlock_H
#define StructuredImplicitBlock_H

/*---------------------------------------------------------------------------*\
    StructuredImplicitBlock - small local conservative-variable block.

    The five hydro variables are ordered as rho, rhoU[0..2], rhoE.  Species
    and auxiliary conserved variables remain scalar in the first block
    preconditioner.  Entries use physical units; the integrator normalizes
    rows and columns with its existing per-cell residual scales before solve.
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>

namespace Foam
{

struct StructuredImplicitBlock
{
    static constexpr int HYDRO_COMPONENTS=5;
    double value[HYDRO_COMPONENTS][HYDRO_COMPONENTS]{};

    void setScalarIdentity(double scalar)
    {
        for(int row=0;row<HYDRO_COMPONENTS;++row)
        for(int column=0;column<HYDRO_COMPONENTS;++column)
            value[row][column]=row==column?scalar:0.0;
    }

    void add(int row,int column,double contribution)
    {
        if(row>=0&&row<HYDRO_COMPONENTS
            &&column>=0&&column<HYDRO_COMPONENTS
            &&std::isfinite(contribution))
            value[row][column]+=contribution;
    }

    bool finite() const
    {
        for(int row=0;row<HYDRO_COMPONENTS;++row)
        for(int column=0;column<HYDRO_COMPONENTS;++column)
            if(!std::isfinite(value[row][column])) return false;
        return true;
    }
};

} // End namespace Foam

#endif
