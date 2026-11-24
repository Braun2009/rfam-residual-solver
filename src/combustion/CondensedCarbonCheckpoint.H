#ifndef CondensedCarbonCheckpoint_H
#define CondensedCarbonCheckpoint_H

#include "CondensedCarbonOperator.H"

#include <cstdint>
#include <string>

namespace Foam
{

struct CondensedCarbonCheckpointMetadata
{
    double time = 0.0;
    std::uint64_t iteration = 0;
    std::string coupledFluidCheckpointId;
};

class CondensedCarbonCheckpoint
{
public:
    static void write(
        const std::string& path,
        const Domain& domain,
        const SectorCondensedCarbonOperator& carbon,
        const CondensedCarbonCheckpointMetadata& metadata);

    static CondensedCarbonCheckpointMetadata read(
        const std::string& path,
        const Domain& domain,
        SectorCondensedCarbonOperator& carbon,
        const std::string& expectedFluidCheckpointId);
};

} // End namespace Foam

#endif
