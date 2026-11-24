#ifndef CondensedCarbonOutput_H
#define CondensedCarbonOutput_H

#include "CondensedCarbonOperator.H"

#include <string>

namespace Foam
{

class CondensedCarbonOutput
{
public:
    static void writeSectorFieldCsv(
        const std::string& path,
        const Domain& domain,
        const SectorMeshMetrics& mesh,
        const SectorCondensedCarbonOperator& carbon);

    static void writeWallInventoryCsv(
        const std::string& path,
        const Domain& domain,
        const SectorMeshMetrics& mesh,
        const SectorCondensedCarbonOperator& carbon);

    static void writeSectorVtk(
        const std::string& path,
        const Domain& domain,
        const SectorMeshMetrics& mesh,
        const SectorCondensedCarbonOperator& carbon);
};

} // End namespace Foam

#endif
