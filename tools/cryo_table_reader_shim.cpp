// CryogenicTransportTable reader shim for offline interpolation assessment.
//
// Usage:
//   cryo_table_reader <table.csv> <O2|C3H8> <T> <rho>
// Output: single-line JSON with {"status","mu","kappa"} on SUCCESS, or
//         {"status": "<name>"} otherwise.
//
// This is a diagnostic/validation tool, not a solver consumer.

#include <iostream>
#include <string>

#include "CryogenicTransportTable.H"

using namespace Foam;

int main(int argc, char** argv)
{
    if (argc < 5)
    {
        std::cerr << "usage: cryo_table_reader <table> <species> <T> <rho>\n";
        return 2;
    }
    const std::string tablePath = argv[1];
    const std::string species = argv[2];
    const double T = std::stod(argv[3]);
    const double rho = std::stod(argv[4]);

    CryogenicTransportTable table;
    if (!table.load(tablePath))
    {
        std::cout << "{\"status\":\"LOAD_FAILED\"}\n";
        return 0;
    }
    const CryogenicTransportTableResult result = table.evaluate(species, T, rho);
    if (result.success())
    {
        std::cout.precision(16);
        std::cout << "{\"status\":\"SUCCESS\",\"mu\":" << result.properties->mu
                  << ",\"kappa\":" << result.properties->kappa << "}\n";
    }
    else
    {
        std::cout << "{\"status\":\""
                  << cryogenicTransportTableStatusName(result.status) << "\"}\n";
    }
    return 0;
}
