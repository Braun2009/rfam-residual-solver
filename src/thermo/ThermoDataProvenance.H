#ifndef ThermoDataProvenance_H
#define ThermoDataProvenance_H

#include <cstddef>
#include <string>

namespace Foam
{
namespace thermoProvenance
{

void appendString(
    std::string& manifest,
    const std::string& key,
    const std::string& value);

void appendSize(
    std::string& manifest,
    const std::string& key,
    std::size_t value);

void appendSigned(
    std::string& manifest,
    const std::string& key,
    long long value);

void appendDouble(
    std::string& manifest,
    const std::string& key,
    double value);

std::string sha256Hex(const std::string& value);

} // End namespace thermoProvenance
} // End namespace Foam

#endif
