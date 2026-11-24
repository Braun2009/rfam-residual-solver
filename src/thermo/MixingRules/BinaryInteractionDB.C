/*---------------------------------------------------------------------------*\
  realFluidAllMach - Binary Interaction Parameters Implementation
\*---------------------------------------------------------------------------*/

#include "BinaryInteractionDB.H"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace Foam
{

namespace
{

BinaryInteractionMetadata legacyMetadata(double value)
{
    BinaryInteractionMetadata metadata;
    metadata.zeroClass = value == 0.0
        ? BinaryInteractionZeroClass::UNVALIDATED_ZERO
        : BinaryInteractionZeroClass::NOT_ZERO;
    return metadata;
}

bool validProvenanceIdentifier(const std::string& value)
{
    bool hasNonWhitespace = false;
    for (const unsigned char character : value)
    {
        if (character < 0x20 || character > 0x7e) return false;
        if (character != 0x20) hasNonWhitespace = true;
    }
    return hasNonWhitespace;
}

BinaryInteractionReadinessStatus readinessStatus(
    const BinaryInteractionEntry* entry,
    const BinaryInteractionApplicability* requiredEnvelope)
{
    if (!entry) return BinaryInteractionReadinessStatus::MISSING;

    if (!std::isfinite(entry->value))
        return BinaryInteractionReadinessStatus::INVALID_METADATA;

    switch (entry->metadata.evidence)
    {
        case BinaryInteractionEvidenceClass::LEGACY_UNPROVENANCED:
            return BinaryInteractionReadinessStatus::LEGACY_UNPROVENANCED;
        case BinaryInteractionEvidenceClass::TEST_ONLY:
            return BinaryInteractionReadinessStatus::TEST_ONLY;
        case BinaryInteractionEvidenceClass::PROVISIONAL:
            return BinaryInteractionReadinessStatus::PROVISIONAL;
        case BinaryInteractionEvidenceClass::VALIDATED:
            break;
    }

    if (!validProvenanceIdentifier(entry->metadata.sourceId)
        || !validProvenanceIdentifier(entry->metadata.datasetId)
        || !entry->metadata.applicability.valid()
        || !std::isfinite(entry->metadata.absoluteUncertainty)
        || entry->metadata.absoluteUncertainty < 0.0)
        return BinaryInteractionReadinessStatus::INVALID_METADATA;

    if (entry->value == 0.0)
    {
        if (entry->metadata.zeroClass
            != BinaryInteractionZeroClass::VALIDATED_EXPLICIT_ZERO)
            return BinaryInteractionReadinessStatus::UNVALIDATED_ZERO;
    }
    else if (entry->metadata.zeroClass
        != BinaryInteractionZeroClass::NOT_ZERO)
    {
        return BinaryInteractionReadinessStatus::INVALID_METADATA;
    }

    if (requiredEnvelope
        && !entry->metadata.applicability.covers(*requiredEnvelope))
        return BinaryInteractionReadinessStatus::OUTSIDE_APPLICABILITY;

    return BinaryInteractionReadinessStatus::READY;
}

BinaryInteractionReadinessAudit auditProductionReadinessImpl(
    const BinaryInteractionDB& database,
    const std::vector<std::string>& potentiallyActiveSpecies,
    const BinaryInteractionApplicability* requiredEnvelope)
{
    BinaryInteractionReadinessAudit audit;
    audit.requiredEnvelopeSpecified = requiredEnvelope != nullptr;
    if (requiredEnvelope) audit.requiredEnvelope = *requiredEnvelope;

    if (requiredEnvelope && !requiredEnvelope->valid())
    {
        audit.requestError = "invalid required applicability envelope";
        return audit;
    }

    std::vector<std::string> species = potentiallyActiveSpecies;
    if (species.empty())
    {
        audit.requestError = "empty species set";
        return audit;
    }
    if (std::any_of(
            species.begin(), species.end(),
            [](const std::string& name) { return name.empty(); }))
    {
        audit.requestError = "empty species name";
        return audit;
    }
    std::sort(species.begin(), species.end());
    if (std::adjacent_find(species.begin(), species.end()) != species.end())
    {
        audit.requestError = "duplicate species name";
        return audit;
    }

    audit.requestValid = true;
    for (std::size_t i = 0; i < species.size(); ++i)
    {
        for (std::size_t j = i + 1; j < species.size(); ++j)
        {
            BinaryInteractionPairReadiness pair;
            pair.first = species[i];
            pair.second = species[j];
            pair.status = readinessStatus(
                database.find(pair.first, pair.second), requiredEnvelope);
            audit.pairs.push_back(std::move(pair));
        }
    }
    return audit;
}

} // End anonymous namespace

const char* binaryInteractionEvidenceClassName(
    BinaryInteractionEvidenceClass evidence)
{
    switch (evidence)
    {
        case BinaryInteractionEvidenceClass::LEGACY_UNPROVENANCED:
            return "LEGACY_UNPROVENANCED";
        case BinaryInteractionEvidenceClass::TEST_ONLY:
            return "TEST_ONLY";
        case BinaryInteractionEvidenceClass::PROVISIONAL:
            return "PROVISIONAL";
        case BinaryInteractionEvidenceClass::VALIDATED:
            return "VALIDATED";
    }
    return "UNKNOWN";
}

const char* binaryInteractionZeroClassName(
    BinaryInteractionZeroClass zeroClass)
{
    switch (zeroClass)
    {
        case BinaryInteractionZeroClass::NOT_ZERO:
            return "NOT_ZERO";
        case BinaryInteractionZeroClass::UNVALIDATED_ZERO:
            return "UNVALIDATED_ZERO";
        case BinaryInteractionZeroClass::VALIDATED_EXPLICIT_ZERO:
            return "VALIDATED_EXPLICIT_ZERO";
    }
    return "UNKNOWN";
}

bool BinaryInteractionApplicability::valid() const
{
    return std::isfinite(minimumTemperature)
        && std::isfinite(maximumTemperature)
        && std::isfinite(minimumPressure)
        && std::isfinite(maximumPressure)
        && minimumTemperature > 0.0
        && maximumTemperature >= minimumTemperature
        && minimumPressure > 0.0
        && maximumPressure >= minimumPressure;
}

bool BinaryInteractionApplicability::covers(
    const BinaryInteractionApplicability& required) const
{
    return valid() && required.valid()
        && minimumTemperature <= required.minimumTemperature
        && maximumTemperature >= required.maximumTemperature
        && minimumPressure <= required.minimumPressure
        && maximumPressure >= required.maximumPressure;
}

bool BinaryInteractionReadinessAudit::productionReady() const
{
    return requestValid && requiredEnvelopeSpecified && std::all_of(
        pairs.begin(), pairs.end(),
        [](const BinaryInteractionPairReadiness& pair)
        {
            return pair.status == BinaryInteractionReadinessStatus::READY;
        });
}

std::size_t BinaryInteractionReadinessAudit::count(
    BinaryInteractionReadinessStatus status) const
{
    return static_cast<std::size_t>(std::count_if(
        pairs.begin(), pairs.end(),
        [status](const BinaryInteractionPairReadiness& pair)
        {
            return pair.status == status;
        }));
}

std::vector<std::pair<std::string, std::string>>
BinaryInteractionReadinessAudit::pairsWithStatus(
    BinaryInteractionReadinessStatus status) const
{
    std::vector<std::pair<std::string, std::string>> result;
    for (const BinaryInteractionPairReadiness& pair : pairs)
    {
        if (pair.status == status)
            result.emplace_back(pair.first, pair.second);
    }
    return result;
}

std::string BinaryInteractionReadinessAudit::report() const
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << std::setprecision(17)
           << "schema=rfam.binary-interaction-readiness/v1\n"
           << "request_valid=" << (requestValid ? 1 : 0) << '\n'
           << "required_envelope="
           << (requiredEnvelopeSpecified ? 1 : 0);
    if (requiredEnvelopeSpecified)
    {
        stream << ',' << requiredEnvelope.minimumTemperature
               << ',' << requiredEnvelope.maximumTemperature
               << ',' << requiredEnvelope.minimumPressure
               << ',' << requiredEnvelope.maximumPressure;
    }
    stream << '\n';
    if (!requestError.empty())
        stream << "request_error=" << requestError << '\n';
    stream << "pair_count=" << pairs.size() << '\n';
    for (const BinaryInteractionPairReadiness& pair : pairs)
    {
        stream << "pair=" << pair.first << '|' << pair.second << '|'
               << binaryInteractionReadinessStatusName(pair.status) << '\n';
    }
    return stream.str();
}

const char* binaryInteractionReadinessStatusName(
    BinaryInteractionReadinessStatus status)
{
    switch (status)
    {
        case BinaryInteractionReadinessStatus::READY:
            return "READY";
        case BinaryInteractionReadinessStatus::MISSING:
            return "MISSING";
        case BinaryInteractionReadinessStatus::LEGACY_UNPROVENANCED:
            return "LEGACY_UNPROVENANCED";
        case BinaryInteractionReadinessStatus::TEST_ONLY:
            return "TEST_ONLY";
        case BinaryInteractionReadinessStatus::PROVISIONAL:
            return "PROVISIONAL";
        case BinaryInteractionReadinessStatus::INVALID_METADATA:
            return "INVALID_METADATA";
        case BinaryInteractionReadinessStatus::UNVALIDATED_ZERO:
            return "UNVALIDATED_ZERO";
        case BinaryInteractionReadinessStatus::OUTSIDE_APPLICABILITY:
            return "OUTSIDE_APPLICABILITY";
    }
    return "UNKNOWN";
}

BinaryInteractionDB::BinaryInteractionDB()
{}

std::pair<std::string, std::string> BinaryInteractionDB::makeKey(
    const std::string& name1,
    const std::string& name2
) const
{
    // Canonical ordering: alphabetically first comes first
    if (name1 < name2)
        return std::make_pair(name1, name2);
    else
        return std::make_pair(name2, name1);
}

void BinaryInteractionDB::add(
    const std::string& name1,
    const std::string& name2,
    double kij
)
{
    add(name1, name2, kij, legacyMetadata(kij));
}

void BinaryInteractionDB::add(
    const std::string& name1,
    const std::string& name2,
    double kij,
    const BinaryInteractionMetadata& metadata)
{
    BinaryInteractionEntry entry;
    entry.value = kij;
    entry.metadata = metadata;
    entries_[makeKey(name1, name2)] = std::move(entry);
}

double BinaryInteractionDB::get(
    const std::string& name1,
    const std::string& name2
) const
{
    const BinaryInteractionEntry* entry = find(name1, name2);
    if (entry) return entry->value;

    // Default: ideal mixing (kij = 0)
    return 0.0;
}

bool BinaryInteractionDB::has(
    const std::string& name1,
    const std::string& name2
) const
{
    return find(name1, name2) != nullptr;
}

const BinaryInteractionEntry* BinaryInteractionDB::find(
    const std::string& name1,
    const std::string& name2) const
{
    const auto it = entries_.find(makeKey(name1, name2));
    return it == entries_.end() ? nullptr : &it->second;
}

BinaryInteractionReadinessAudit
BinaryInteractionDB::auditProductionReadiness(
    const std::vector<std::string>& potentiallyActiveSpecies) const
{
    return auditProductionReadinessImpl(
        *this, potentiallyActiveSpecies, nullptr);
}

BinaryInteractionReadinessAudit
BinaryInteractionDB::auditProductionReadiness(
    const std::vector<std::string>& potentiallyActiveSpecies,
    const BinaryInteractionApplicability& requiredEnvelope) const
{
    return auditProductionReadinessImpl(
        *this, potentiallyActiveSpecies, &requiredEnvelope);
}

BinaryInteractionDB BinaryInteractionDB::createDefault()
{
    BinaryInteractionDB db;

    // These legacy hard-coded values preserve existing diagnostic numerics.
    // They have no source/range/uncertainty metadata and are not production
    // evidence, including the explicitly stored N2/O2 zero.
    db.add("N2", "O2",   0.0000);  // Ideal
    db.add("N2", "CH4",  0.0311);
    db.add("N2", "C2H6", 0.0515);
    db.add("N2", "C3H8", 0.0800);
    db.add("N2", "CO2",  -0.0170);

    db.add("O2", "CH4",  0.0230);
    db.add("O2", "C2H6", 0.0400);
    db.add("O2", "C3H8", 0.0550);
    db.add("O2", "H2",   0.0100);
    db.add("O2", "CO2",  0.0200);

    db.add("CH4", "C2H6", 0.0026);
    db.add("CH4", "C3H8", 0.0067);
    db.add("CH4", "CO2",  0.0919);
    db.add("CH4", "H2O",  0.4850);

    db.add("C2H6", "C3H8", 0.0011);
    db.add("C3H8", "CO2",  0.1300);

    db.add("H2", "CH4",   0.0340);
    db.add("H2", "CO",    0.1200);

    db.add("CO", "CO2",   0.0700);
    db.add("CO2", "H2O",  0.1900);

    return db;
}

} // End namespace Foam
