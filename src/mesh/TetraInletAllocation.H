#ifndef RFAM_TETRA_INLET_ALLOCATION_H
#define RFAM_TETRA_INLET_ALLOCATION_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{
// Redistribute a requested mass flow across faces without exceeding the
// local supply capacities. Among unconstrained faces the mass flux is uniform.
inline std::vector<double> allocateTetraInletFlow(
    const std::vector<double>& areas, const std::vector<double>& capacities,
    double requested)
{
    if(areas.size()!=capacities.size()||!std::isfinite(requested)||requested<0.)
        throw std::invalid_argument("invalid tetra inlet allocation");
    double totalCapacity=0., totalArea=0.;
    for(std::size_t i=0;i<areas.size();++i)
    {
        if(!std::isfinite(areas[i])||areas[i]<=0.||!std::isfinite(capacities[i])||capacities[i]<0.)
            throw std::invalid_argument("invalid tetra inlet face capacity");
        totalArea+=areas[i];totalCapacity+=capacities[i];
    }
    std::vector<double> flow(areas.size(),0.);
    if(requested==0.||totalCapacity==0.)return flow;
    if(requested>=totalCapacity)return capacities;
    double low=0.,high=requested/totalArea;
    for(std::size_t i=0;i<areas.size();++i)
        if(capacities[i]>0.)high=std::max(high,capacities[i]/areas[i]);
    for(int iter=0;iter<80;++iter)
    {
        const double mid=low+.5*(high-low);
        double delivered=0.;
        for(std::size_t i=0;i<areas.size();++i)
            delivered+=std::min(capacities[i],mid*areas[i]);
        if(delivered<requested)low=mid;else high=mid;
    }
    for(std::size_t i=0;i<areas.size();++i)
        flow[i]=std::min(capacities[i],high*areas[i]);
    return flow;
}
}
#endif
