/*---------------------------------------------------------------------------*\
  realFluidAllMach - Binary Interaction Parameters Database
\*---------------------------------------------------------------------------*/

#ifndef BinaryInteractionDB_H
#define BinaryInteractionDB_H

#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

enum class BinaryInteractionEvidenceClass
{
    LEGACY_UNPROVENANCED,
    TEST_ONLY,
    PROVISIONAL,
    VALIDATED
};

const char* binaryInteractionEvidenceClassName(
    BinaryInteractionEvidenceClass evidence);

enum class BinaryInteractionZeroClass
{
    NOT_ZERO,
    UNVALIDATED_ZERO,
    VALIDATED_EXPLICIT_ZERO
};

const char* binaryInteractionZeroClassName(
    BinaryInteractionZeroClass zeroClass);

struct BinaryInteractionApplicability
{
    double minimumTemperature = std::numeric_limits<double>::quiet_NaN();
    double maximumTemperature = std::numeric_limits<double>::quiet_NaN();
    double minimumPressure = std::numeric_limits<double>::quiet_NaN();
    double maximumPressure = std::numeric_limits<double>::quiet_NaN();

    bool valid() const;
    bool covers(const BinaryInteractionApplicability& required) const;
};

struct BinaryInteractionMetadata
{
    BinaryInteractionEvidenceClass evidence =
        BinaryInteractionEvidenceClass::LEGACY_UNPROVENANCED;
    BinaryInteractionZeroClass zeroClass =
        BinaryInteractionZeroClass::NOT_ZERO;
    std::string sourceId;
    std::string datasetId;
    BinaryInteractionApplicability applicability;
    double absoluteUncertainty = std::numeric_limits<double>::quiet_NaN();
};

struct BinaryInteractionEntry
{
    double value = std::numeric_limits<double>::quiet_NaN();
    BinaryInteractionMetadata metadata;
};

enum class BinaryInteractionReadinessStatus
{
    READY,
    MISSING,
    LEGACY_UNPROVENANCED,
    TEST_ONLY,
    PROVISIONAL,
    INVALID_METADATA,
    UNVALIDATED_ZERO,
    OUTSIDE_APPLICABILITY
};

const char* binaryInteractionReadinessStatusName(
    BinaryInteractionReadinessStatus status);

struct BinaryInteractionPairReadiness
{
    std::string first;
    std::string second;
    BinaryInteractionReadinessStatus status =
        BinaryInteractionReadinessStatus::MISSING;
};

struct BinaryInteractionReadinessAudit
{
    bool requestValid = false;
    bool requiredEnvelopeSpecified = false;
    BinaryInteractionApplicability requiredEnvelope;
    std::string requestError;
    std::vector<BinaryInteractionPairReadiness> pairs;

    bool productionReady() const;
    std::size_t count(BinaryInteractionReadinessStatus status) const;
    std::vector<std::pair<std::string, std::string>> pairsWithStatus(
        BinaryInteractionReadinessStatus status) const;
    std::string report() const;
};

class BinaryInteractionDB
{
private:
    std::map<
        std::pair<std::string, std::string>,
        BinaryInteractionEntry> entries_;

    std::pair<std::string, std::string> makeKey(
        const std::string& name1,
        const std::string& name2
    ) const;

public:
    BinaryInteractionDB();

    void add(const std::string& name1, const std::string& name2, double kij);

    void add(
        const std::string& name1,
        const std::string& name2,
        double kij,
        const BinaryInteractionMetadata& metadata);

    double get(const std::string& name1, const std::string& name2) const;

    bool has(const std::string& name1, const std::string& name2) const;

    const BinaryInteractionEntry* find(
        const std::string& name1,
        const std::string& name2) const;

    BinaryInteractionReadinessAudit auditProductionReadiness(
        const std::vector<std::string>& potentiallyActiveSpecies) const;

    BinaryInteractionReadinessAudit auditProductionReadiness(
        const std::vector<std::string>& potentiallyActiveSpecies,
        const BinaryInteractionApplicability& requiredEnvelope) const;

    // Common binary pairs
    static BinaryInteractionDB createDefault();
};

} // End namespace Foam

#endif
