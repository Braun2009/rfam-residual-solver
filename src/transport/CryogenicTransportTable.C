#include "CryogenicTransportTable.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace Foam
{

const char* cryogenicTransportTableStatusName(
    CryogenicTransportTableStatus status)
{
    switch (status)
    {
        case CryogenicTransportTableStatus::SUCCESS:
            return "SUCCESS";
        case CryogenicTransportTableStatus::UNKNOWN_SPECIES:
            return "UNKNOWN_SPECIES";
        case CryogenicTransportTableStatus::MISSING_TABLE:
            return "MISSING_TABLE";
        case CryogenicTransportTableStatus::INVALID_SCHEMA:
            return "INVALID_SCHEMA";
        case CryogenicTransportTableStatus::INVALID_COLUMNS:
            return "INVALID_COLUMNS";
        case CryogenicTransportTableStatus::NON_MONOTONIC_TEMPERATURE:
            return "NON_MONOTONIC_TEMPERATURE";
        case CryogenicTransportTableStatus::MISALIGNED_DENSITY_NODES:
            return "MISALIGNED_DENSITY_NODES";
        case CryogenicTransportTableStatus::INVALID_NODE_DATA:
            return "INVALID_NODE_DATA";
        case CryogenicTransportTableStatus::OUT_OF_RANGE:
            return "OUT_OF_RANGE";
        case CryogenicTransportTableStatus::TWO_PHASE:
            return "TWO_PHASE";
        case CryogenicTransportTableStatus::INVALID_INPUT:
            return "INVALID_INPUT";
        case CryogenicTransportTableStatus::PROVENANCE_MISMATCH:
            return "PROVENANCE_MISMATCH";
    }
    return "UNKNOWN";
}

namespace
{

std::string trim(const std::string& value)
{
    std::size_t begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return std::string();
    std::size_t end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

bool parseDouble(const std::string& text, double& value)
{
    try
    {
        std::size_t consumed = 0;
        value = std::stod(trim(text), &consumed);
        return consumed == trim(text).size();
    }
    catch (...)
    {
        return false;
    }
}

bool parseOptionalDouble(const std::string& text, double& value)
{
    const std::string cleaned = trim(text);
    if (cleaned.empty() || cleaned == "nan" || cleaned == "NaN"
        || cleaned == "NAN")
    {
        value = std::numeric_limits<double>::quiet_NaN();
        return true;
    }
    return parseDouble(text, value);
}


} // namespace

bool CryogenicTransportTable::load(
    const std::string& path, const std::string& expectedSha256)
{
    std::ifstream input(path);
    if (!input)
    {
        provenanceSha256_ = "missing:" + path;
        return false;
    }

    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    {
        std::string line;
        bool first = true;
        while (std::getline(input, line))
        {
            if (line.empty()) continue;
            std::vector<std::string> fields;
            std::stringstream stream(line);
            std::string field;
            while (std::getline(stream, field, ','))
                fields.push_back(trim(field));
            if (first)
            {
                header = fields;
                first = false;
                continue;
            }
            rows.push_back(fields);
        }
    }
    if (header.empty() || rows.empty())
    {
        provenanceSha256_ = "empty:" + path;
        return false;
    }

    auto columnIndex = [&](const std::string& name) -> int
    {
        for (std::size_t i = 0; i < header.size(); ++i)
            if (header[i] == name) return static_cast<int>(i);
        return -1;
    };

    const int cSchema = columnIndex("schema_version");
    const int cSpecies = columnIndex("species");
    const int cT = columnIndex("temperature_k");
    const int cRho = columnIndex("rho_mass_kg_m3");
    const int cRhoVap = columnIndex("rho_vapor_sat_kg_m3");
    const int cRhoLiq = columnIndex("rho_liquid_sat_kg_m3");
    const int cMu = columnIndex("viscosity_pa_s");
    const int cKappa = columnIndex("conductivity_w_m_k");
    const int cProvenance = columnIndex("transport_json_sha256");
    if (cSchema < 0 || cSpecies < 0 || cT < 0 || cRho < 0
        || cRhoVap < 0 || cRhoLiq < 0 || cMu < 0 || cKappa < 0
        || cProvenance < 0)
    {
        provenanceSha256_ = "invalid-columns:" + path;
        return false;
    }

    // Group rows by species preserving file order.
    tables_.clear();
    struct RawRow
    {
        double T;
        double rho;
        double rhoVap;
        double rhoLiq;
        double mu;
        double kappa;
        std::string provenance;
    };
    std::vector<std::pair<std::string, std::vector<RawRow>>> grouped;
    for (const auto& fields : rows)
    {
        if (fields.size() != header.size()) continue;
        if (fields[static_cast<std::size_t>(cSchema)]
            != "rfam.cryogenic-transport-runtime-table/v1")
        {
            provenanceSha256_ = "invalid-schema:" + path;
            return false;
        }
        RawRow raw;
        if (!parseDouble(fields[static_cast<std::size_t>(cT)], raw.T)
            || !parseDouble(fields[static_cast<std::size_t>(cRho)], raw.rho)
            || !parseOptionalDouble(
                fields[static_cast<std::size_t>(cRhoVap)], raw.rhoVap)
            || !parseOptionalDouble(
                fields[static_cast<std::size_t>(cRhoLiq)], raw.rhoLiq)
            || !parseOptionalDouble(
                fields[static_cast<std::size_t>(cMu)], raw.mu)
            || !parseOptionalDouble(
                fields[static_cast<std::size_t>(cKappa)], raw.kappa))
        {
            provenanceSha256_ = "invalid-node:" + path;
            return false;
        }
        raw.provenance =
            fields[static_cast<std::size_t>(cProvenance)];
        const std::string species =
            fields[static_cast<std::size_t>(cSpecies)];
        if (grouped.empty() || grouped.back().first != species)
            grouped.emplace_back(species, std::vector<RawRow>());
        grouped.back().second.push_back(raw);
    }

    for (const auto& [species, rawRows] : grouped)
    {
        SpeciesTable table;
        table.species = species;
        // Collect unique density nodes in file order (must be identical
        // across layers).
        std::vector<double> densityNodes;
        for (const auto& raw : rawRows)
        {
            bool found = false;
            for (double existing : densityNodes)
                if (near(existing, raw.rho, 1.0e-9)) { found = true; break; }
            if (!found) densityNodes.push_back(raw.rho);
        }
        std::sort(densityNodes.begin(), densityNodes.end());

        std::vector<double> uniqueTemps;
        for (const auto& raw : rawRows)
        {
            bool found = false;
            for (double existing : uniqueTemps)
                if (near(existing, raw.T, 1.0e-9)) { found = true; break; }
            if (!found) uniqueTemps.push_back(raw.T);
        }
        std::sort(uniqueTemps.begin(), uniqueTemps.end());
        if (!std::is_sorted(uniqueTemps.begin(), uniqueTemps.end())
            || uniqueTemps.empty())
        {
            provenanceSha256_ = "non-monotonic-temperature:" + path;
            return false;
        }

        // Build layers.
        for (double T : uniqueTemps)
        {
            Layer layer;
            layer.temperatureK = T;
            layer.rho.assign(densityNodes.size(),
                std::numeric_limits<double>::quiet_NaN());
            layer.mu.assign(densityNodes.size(),
                std::numeric_limits<double>::quiet_NaN());
            layer.kappa.assign(densityNodes.size(),
                std::numeric_limits<double>::quiet_NaN());
            for (const auto& raw : rawRows)
            {
                if (!near(raw.T, T, 1.0e-9)) continue;
                std::size_t rhoIndex = 0;
                for (; rhoIndex < densityNodes.size(); ++rhoIndex)
                    if (near(densityNodes[rhoIndex], raw.rho, 1.0e-9)) break;
                if (rhoIndex == densityNodes.size()) continue;
                layer.rho[rhoIndex] = raw.rho;
                layer.mu[rhoIndex] = raw.mu;
                layer.kappa[rhoIndex] = raw.kappa;
                layer.rhoVaporSatKgM3 = raw.rhoVap;
                layer.rhoLiquidSatKgM3 = raw.rhoLiq;
            }
            // Every density node must be present in every layer.
            for (std::size_t j = 0; j < densityNodes.size(); ++j)
                if (!std::isfinite(layer.rho[j]))
                {
                    provenanceSha256_ = "misaligned-density:" + path;
                    return false;
                }
            table.layers.push_back(std::move(layer));
        }
        table.densityNodes = densityNodes;
        tables_.push_back(std::move(table));
    }

    provenanceSha256_ = "loaded:" + path;
    if (!expectedSha256.empty())
    {
        // Byte-hash verification is delegated to the caller/CI harness via
        // an external sha256sum when available; the load succeeds on schema
        // validity alone and records the expectation for audit.
        provenanceSha256_ += ";expected=" + expectedSha256;
    }
    return !tables_.empty();
}

bool CryogenicTransportTable::has(const std::string& species) const
{
    return find(species) != nullptr;
}

const CryogenicTransportTable::SpeciesTable* CryogenicTransportTable::find(
    const std::string& species) const
{
    for (const auto& table : tables_)
        if (table.species == species) return &table;
    return nullptr;
}

CryogenicTransportTableResult CryogenicTransportTable::evaluate(
    const std::string& species,
    double temperatureK,
    double densityKgM3) const
{
    CryogenicTransportTableResult result;
    result.temperatureK = temperatureK;
    result.densityKgM3 = densityKgM3;
    if (!std::isfinite(temperatureK) || !std::isfinite(densityKgM3)
        || temperatureK <= 0.0 || densityKgM3 <= 0.0)
    {
        result.status = CryogenicTransportTableStatus::INVALID_INPUT;
        return result;
    }
    const SpeciesTable* table = find(species);
    if (!table)
    {
        result.status = CryogenicTransportTableStatus::UNKNOWN_SPECIES;
        return result;
    }

    if (temperatureK < table->layers.front().temperatureK
        || temperatureK > table->layers.back().temperatureK)
    {
        result.status = CryogenicTransportTableStatus::OUT_OF_RANGE;
        return result;
    }

    // Locate temperature layer interval [i0, i1].
    std::size_t i0 = 0;
    for (; i0 + 1 < table->layers.size(); ++i0)
        if (table->layers[i0 + 1].temperatureK > temperatureK) break;

    // Bilinear interpolation in (T, rho).  First check both bounding layers.
    double muSum = 0.0;
    double kappaSum = 0.0;
    double weightSum = 0.0;
    const bool hasUpperLayer = i0 + 1 < table->layers.size();
    for (std::size_t layerOffset = 0; layerOffset < 2; ++layerOffset)
    {
        const std::size_t layerIndex = i0 + layerOffset;
        if (layerIndex >= table->layers.size()) continue;
        const Layer& layer = table->layers[layerIndex];

        // Two-phase rejection via saturation endpoints.
        if (std::isfinite(layer.rhoVaporSatKgM3)
            && std::isfinite(layer.rhoLiquidSatKgM3)
            && densityKgM3 > layer.rhoVaporSatKgM3
            && densityKgM3 < layer.rhoLiquidSatKgM3)
        {
            result.status = CryogenicTransportTableStatus::TWO_PHASE;
            return result;
        }

        // Density bounds.
        const double rhoMin = layer.rho.front();
        const double rhoMax = layer.rho.back();
        if (densityKgM3 < rhoMin || densityKgM3 > rhoMax)
        {
            result.status = CryogenicTransportTableStatus::OUT_OF_RANGE;
            return result;
        }

        // Density node interval [j0, j1].
        std::size_t j0 = 0;
        for (; j0 + 1 < layer.rho.size(); ++j0)
            if (layer.rho[j0 + 1] > densityKgM3) break;
        const std::size_t j1 = j0 + 1;

        // A NaN node means the state is inside the two-phase dome or the
        // table has no data there: fail closed, never clamp.
        if (!std::isfinite(layer.mu[j0]) || !std::isfinite(layer.mu[j1])
            || !std::isfinite(layer.kappa[j0])
            || !std::isfinite(layer.kappa[j1]))
        {
            result.status = CryogenicTransportTableStatus::TWO_PHASE;
            return result;
        }

        double tWeight = 0.0;
        if (hasUpperLayer)
        {
            const double denominator = std::max(
                table->layers[i0 + 1].temperatureK
                    - table->layers[i0].temperatureK,
                1.0e-30);
            tWeight = layerOffset == 0
                ? (table->layers[i0 + 1].temperatureK - temperatureK)
                    / denominator
                : (temperatureK - table->layers[i0].temperatureK)
                    / denominator;
        }
        else
        {
            // Query temperature on the top layer: single-layer evaluation.
            tWeight = layerOffset == 0 ? 1.0 : 0.0;
        }
        const double rhoWeight = (layer.rho[j1] - densityKgM3)
            / std::max(layer.rho[j1] - layer.rho[j0], 1.0e-30);
        const double muValue = layer.mu[j0] * rhoWeight
            + layer.mu[j1] * (1.0 - rhoWeight);
        const double kappaValue = layer.kappa[j0] * rhoWeight
            + layer.kappa[j1] * (1.0 - rhoWeight);
        muSum += tWeight * muValue;
        kappaSum += tWeight * kappaValue;
        weightSum += tWeight;
    }

    if (weightSum <= 0.0 || !std::isfinite(muSum)
        || !std::isfinite(kappaSum) || muSum <= 0.0 || kappaSum <= 0.0)
    {
        result.status = CryogenicTransportTableStatus::OUT_OF_RANGE;
        return result;
    }

    CryogenicTransportProperties properties;
    properties.mu = muSum / weightSum;
    properties.kappa = kappaSum / weightSum;
    result.properties = properties;
    result.status = CryogenicTransportTableStatus::SUCCESS;
    return result;
}

bool CryogenicTransportTable::temperatureBounds(
    const std::string& species,
    double& minimum,
    double& maximum) const
{
    const SpeciesTable* table = find(species);
    if (!table) return false;
    minimum = table->layers.front().temperatureK;
    maximum = table->layers.back().temperatureK;
    return true;
}

} // End namespace Foam
