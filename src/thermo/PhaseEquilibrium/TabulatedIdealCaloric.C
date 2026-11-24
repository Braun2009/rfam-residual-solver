#include "TabulatedIdealCaloric.H"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include "ThermoDataProvenance.H"

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

using IdealCaloricNode = TabulatedIdealCaloricNode;

struct IdealCaloricTable
{
    const char* species;
    const IdealCaloricNode* nodes;
    std::size_t size;
};

constexpr IdealCaloricNode o2Nodes[] =
{
    {8.00000000000000000e+01, -1.98799922488118580e+05, 9.09940514932144538e+02},
    {8.50000000000000000e+01, -1.94250295175434730e+05, 9.09911264832530037e+02},
    {9.00000000000000000e+01, -1.89700802165195579e+05, 9.09886629027717277e+02},
    {9.50000000000000000e+01, -1.85151422641482088e+05, 9.09865746490661763e+02},
    {1.00000000000000000e+02, -1.80602139531058172e+05, 9.09847970870067797e+02},
    {1.05000000000000000e+02, -1.76052938566684374e+05, 9.09832819167641105e+02},
    {1.10000000000000000e+02, -1.71503807561597874e+05, 9.09819937704991275e+02},
    {1.15000000000000000e+02, -1.66954735820258618e+05, 9.09809081082159878e+02},
    {1.20000000000000000e+02, -1.62405713627953315e+05, 9.09800101278621810e+02},
    {1.25000000000000000e+02, -1.57856731774041546e+05, 9.09792944776981926e+02},
    {1.30000000000000000e+02, -1.53307781073335063e+05, 9.09787655895163425e+02},
    {1.35000000000000000e+02, -1.48758851858888374e+05, 9.09784384605184641e+02},
    {1.40000000000000000e+02, -1.44209933428022079e+05, 9.09783397140408624e+02},
    {1.45000000000000000e+02, -1.39661013431799394e+05, 9.09785087738456014e+02},
    {1.50000000000000000e+02, -1.35112077206192567e+05, 9.09789989979145389e+02},
    {1.55000000000000000e+02, -1.30563107050434279e+05, 9.09798786367320304e+02},
    {1.60000000000000000e+02, -1.26014081464179239e+05, 9.09812315067935288e+02},
    {1.65000000000000000e+02, -1.21464974359825021e+05, 9.09831573001529591e+02},
    {1.70000000000000000e+02, -1.16915754269511497e+05, 9.09857714824316190e+02},
    {1.75000000000000000e+02, -1.12366383567923331e+05, 9.09892047622311338e+02},
    {1.80000000000000000e+02, -1.07816817732173949e+05, 9.09936021422293493e+02},
    {1.85000000000000000e+02, -1.03267004658945516e+05, 9.09991215849648938e+02},
    {1.90000000000000000e+02, -9.87168840569516469e+04, 9.10059323436865611e+02},
    {1.95000000000000000e+02, -9.41663869299524813e+04, 9.10142130205013814e+02},
    {2.00000000000000000e+02, -8.96154351622519316e+04, 9.10241494206968468e+02},
    {2.05000000000000000e+02, -8.50639412150962453e+04, 9.10359322741521169e+02},
    {2.10000000000000000e+02, -8.05118079388768820e+04, 9.10497548929777622e+02},
    {2.15000000000000000e+02, -7.59589285026897851e+04, 9.10658108297944750e+02},
    {2.20000000000000000e+02, -7.14051864397501922e+04, 9.10842915942221907e+02},
    {2.25000000000000000e+02, -6.68504558044798323e+04, 9.11053844769726993e+02},
    {2.30000000000000000e+02, -6.22946014348362223e+04, 9.11292705220909056e+02},
    {2.35000000000000000e+02, -5.77374793116499786e+04, 9.11561226789150624e+02},
    {2.40000000000000000e+02, -5.31789370053767052e+04, 9.11861041566512426e+02},
    {2.45000000000000000e+02, -4.86188141997303755e+04, 9.12193669963862703e+02},
    {2.50000000000000000e+02, -4.40569432811081642e+04, 9.12560508681070246e+02},
    {2.55000000000000000e+02, -3.94931499825018109e+04, 9.12962820939675794e+02},
    {2.60000000000000000e+02, -3.49272540706649015e+04, 9.13401728937034704e+02},
    {2.65000000000000000e+02, -3.03590700656239933e+04, 9.13878208437274566e+02},
    {2.70000000000000000e+02, -2.57884079821352789e+04, 9.14393085380125854e+02},
    {2.75000000000000000e+02, -2.12150740833486489e+04, 9.14947034363040871e+02},
    {2.80000000000000000e+02, -1.66388716377144447e+04, 9.15540578834152484e+02},
    {2.85000000000000000e+02, -1.20596016710076947e+04, 9.16174092822563466e+02},
    {2.90000000000000000e+02, -7.47706370622949908e+03, 9.16847804027229131e+02},
    {2.95000000000000000e+02, -2.89105648504151031e+03, 9.17561798085309647e+02},
    {2.98149999999999977e+02, 0.00000000000000000e+00, 9.18032280063897815e+02},
    {3.00000000000000000e+02, 1.69862133472319692e+03, 9.18316023844432721e+02},
};

constexpr IdealCaloricNode c3h8Nodes[] =
{
    {8.00000000000000000e+01, -2.71744578360415064e+05, 8.57617226424870069e+02},
    {8.25000000000000000e+01, -2.69588404233781039e+05, 8.67341324682845652e+02},
    {8.50000000000000000e+01, -2.67407794017461070e+05, 8.77158851095276873e+02},
    {8.75000000000000000e+01, -2.65202566427857324e+05, 8.87028803712855051e+02},
    {9.00000000000000000e+01, -2.62972636484813411e+05, 8.96915135660336773e+02},
    {9.25000000000000000e+01, -2.60718003217751218e+05, 9.06786661846853235e+02},
    {9.50000000000000000e+01, -2.58438737728946842e+05, 9.16616871490552512e+02},
    {9.75000000000000000e+01, -2.56134971814550401e+05, 9.26383673467330368e+02},
    {1.00000000000000000e+02, -2.53806887282781710e+05, 9.36069096508902248e+02},
    {1.02500000000000000e+02, -2.51454706059263088e+05, 9.45658961907795970e+02},
    {1.05000000000000000e+02, -2.49078681130105397e+05, 9.55142542649961456e+02},
    {1.07500000000000000e+02, -2.46679088342609291e+05, 9.64512219755813248e+02},
    {1.10000000000000000e+02, -2.44256219059842697e+05, 9.73763144009642019e+02},
    {1.12500000000000000e+02, -2.41810373647663393e+05, 9.82892909128148858e+02},
    {1.15000000000000000e+02, -2.39341855759858561e+05, 9.91901240693558520e+02},
    {1.17500000000000000e+02, -2.36850967378059984e+05, 1.00078970379206953e+03},
    {1.20000000000000000e+02, -2.34338004557164095e+05, 1.00956143119768217e+03},
    {1.22500000000000000e+02, -2.31803253823520616e+05, 1.01822087307565221e+03},
    {1.25000000000000000e+02, -2.29246989171589667e+05, 1.02677356850711340e+03},
    {1.27500000000000000e+02, -2.26669469604687765e+05, 1.03522593862169697e+03},
    {1.30000000000000000e+02, -2.24070937166484247e+05, 1.04358510073893376e+03},
    {1.32500000000000000e+02, -2.21451615411776584e+05, 1.05185870263758397e+03},
    {1.35000000000000000e+02, -2.18811708267541253e+05, 1.06005477587466339e+03},
    {1.37500000000000000e+02, -2.16151399238120299e+05, 1.06818160694624180e+03},
    {1.40000000000000000e+02, -2.13470850911540678e+05, 1.07624762500629413e+03},
    {1.42500000000000000e+02, -2.10770204727205215e+05, 1.08426130482668486e+03},
    {1.45000000000000000e+02, -2.08049580968508322e+05, 1.09223108368134626e+03},
    {1.47500000000000000e+02, -2.05309078947184433e+05, 1.10016529086316041e+03},
    {1.50000000000000000e+02, -2.02548777349378390e+05, 1.10807208858656941e+03},
    {1.55000000000000000e+02, -1.96968990036587813e+05, 1.12383498479852960e+03},
    {1.60000000000000000e+02, -1.91310456877466349e+05, 1.13958008886865377e+03},
    {1.65000000000000000e+02, -1.85573126214528573e+05, 1.15536276949999433e+03},
    {1.70000000000000000e+02, -1.79756683070479310e+05, 1.17123292679329575e+03},
    {1.75000000000000000e+02, -1.73860577151748294e+05, 1.18723477060075015e+03},
    {1.80000000000000000e+02, -1.67884051384869090e+05, 1.20340681697592277e+03},
    {1.85000000000000000e+02, -1.61826170062833233e+05, 1.21978203825178321e+03},
    {1.90000000000000000e+02, -1.55685845964195149e+05, 1.23638811656556913e+03},
    {1.95000000000000000e+02, -1.49461866027736978e+05, 1.25324776262100045e+03},
    {2.00000000000000000e+02, -1.43152915331228054e+05, 1.27037907125629658e+03},
    {2.05000000000000000e+02, -1.36757599244692188e+05, 1.28779589319136562e+03},
    {2.10000000000000000e+02, -1.30274463715830352e+05, 1.30550820842499184e+03},
    {2.15000000000000000e+02, -1.23702013705693302e+05, 1.32352249141732068e+03},
    {2.20000000000000000e+02, -1.17038729832847661e+05, 1.34184206168529840e+03},
    {2.25000000000000000e+02, -1.10283083309415029e+05, 1.36046741599419170e+03},
    {2.30000000000000000e+02, -1.03433549266616232e+05, 1.37939654015064866e+03},
    {2.35000000000000000e+02, -9.64886185740798246e+04, 1.39862519966322225e+03},
    {2.40000000000000000e+02, -8.94468082585902885e+04, 1.41814720937509128e+03},
    {2.45000000000000000e+02, -8.23066706259989878e+04, 1.43795468270320134e+03},
    {2.50000000000000000e+02, -7.50668011860025581e+04, 1.45803826142619891e+03},
    {2.55000000000000000e+02, -6.77258454743409529e+04, 1.47838732711798889e+03},
    {2.60000000000000000e+02, -6.02825048613223480e+04, 1.49899019537536651e+03},
    {2.65000000000000000e+02, -5.27355414298515534e+04, 1.51983429397446912e+03},
    {2.70000000000000000e+02, -4.50837820005854592e+04, 1.54090632603853282e+03},
    {2.75000000000000000e+02, -3.73261213765987195e+04, 1.56219241922751780e+03},
    {2.80000000000000000e+02, -2.94615248750908067e+04, 1.58367826188103550e+03},
    {2.85000000000000000e+02, -2.14890302092016209e+04, 1.60534922696779586e+03},
    {2.90000000000000000e+02, -1.34077487789223669e+04, 1.62719048462220394e+03},
    {2.95000000000000000e+02, -5.21686642635206226e+03, 1.64918710398416511e+03},
    {2.98149999999999977e+02, 0.00000000000000000e+00, 1.66311801759799960e+03},
    {3.00000000000000000e+02, 3.08435629288887139e+03, 1.67132414500254163e+03},
    {3.05000000000000000e+02, 1.14965842734586913e+04, 1.69358674081578602e+03},
    {3.10000000000000000e+02, 2.00204084136458114e+04, 1.71596017128435687e+03},
    {3.15000000000000000e+02, 2.86563465153956786e+04, 1.73842992821746316e+03},
    {3.20000000000000000e+02, 3.74048444907899247e+04, 1.76098177281028779e+03},
    {3.25000000000000000e+02, 4.62662778372695902e+04, 1.78360178578589785e+03},
    {3.30000000000000000e+02, 5.52409533463243861e+04, 1.80627641071739640e+03},
    {3.35000000000000000e+02, 6.43291110127052525e+04, 1.82899249098962127e+03},
    {3.40000000000000000e+02, 7.35309261134689441e+04, 1.85173730084485464e+03},
    {3.45000000000000000e+02, 8.28465114283661824e+04, 1.87449857094313347e+03},
    {3.50000000000000000e+02, 9.22759195751855150e+04, 1.89726450885406143e+03},
};

constexpr IdealCaloricTable tables[] =
{
    {"O2", o2Nodes, sizeof(o2Nodes)/sizeof(o2Nodes[0])},
    {"C3H8", c3h8Nodes, sizeof(c3h8Nodes)/sizeof(c3h8Nodes[0])},
};

const IdealCaloricTable* findTable(const std::string& species)
{
    for (const IdealCaloricTable& table : tables)
    {
        if (species == table.species) return &table;
    }
    return nullptr;
}

TabulatedIdealCaloricProperties interpolate(
    const IdealCaloricNode& lower,
    const IdealCaloricNode& upper,
    double temperature)
{
    const double dt = upper.temperature - lower.temperature;
    const double u = (temperature - lower.temperature)/dt;
    const double u2 = u*u;
    const double u3 = u2*u;

    const double h00 = 2.0*u3 - 3.0*u2 + 1.0;
    const double h10 = u3 - 2.0*u2 + u;
    const double h01 = -2.0*u3 + 3.0*u2;
    const double h11 = u3 - u2;

    const double dh00 = 6.0*u2 - 6.0*u;
    const double dh10 = 3.0*u2 - 4.0*u + 1.0;
    const double dh01 = -6.0*u2 + 6.0*u;
    const double dh11 = 3.0*u2 - 2.0*u;

    TabulatedIdealCaloricProperties properties;
    properties.temperature = temperature;
    properties.sensibleEnthalpy =
        h00*lower.sensibleEnthalpy
      + h10*dt*lower.cp
      + h01*upper.sensibleEnthalpy
      + h11*dt*upper.cp;
    properties.cp =
        (dh00*lower.sensibleEnthalpy
       + dh01*upper.sensibleEnthalpy)/dt
      + dh10*lower.cp
      + dh11*upper.cp;
    return properties;
}

double integrateCpOverTemperature(
    const IdealCaloricNode& lower,
    const IdealCaloricNode& upper,
    double startTemperature,
    double endTemperature)
{
    if (!(startTemperature > 0.0 && endTemperature > startTemperature
        && startTemperature >= lower.temperature
        && endTemperature <= upper.temperature))
        return std::numeric_limits<double>::quiet_NaN();

    const long double t0 = lower.temperature;
    const long double dt = upper.temperature - lower.temperature;
    const long double cp0 = lower.cp;
    const long double cp1 = upper.cp;
    const long double middleTemperature = t0 + 0.5L*dt;
    const long double cpMiddle = interpolate(
        lower, upper, static_cast<double>(middleTemperature)).cp;

    // cp(T)=a*(T-T0)^2+b*(T-T0)+c exactly for the cubic Hermite h(T).
    const long double d1 = cp1 - cp0;
    const long double dm = cpMiddle - cp0;
    const long double a = 2.0L*(d1 - 2.0L*dm)/(dt*dt);
    const long double b = (d1 - a*dt*dt)/dt;
    const long double c = cp0;

    // Convert to A*T^2+B*T+C and integrate cp/T analytically.
    const long double A = a;
    const long double B = b - 2.0L*a*t0;
    const long double C = a*t0*t0 - b*t0 + c;
    const long double begin = startTemperature;
    const long double end = endTemperature;
    const long double integral = 0.5L*A*(end*end - begin*begin)
        + B*(end - begin) + C*std::log(end/begin);
    return static_cast<double>(integral);
}

double sensibleEntropy(
    const IdealCaloricTable& table,
    double temperature)
{
    const double reference = TabulatedIdealCaloric::referenceTemperature;
    if (temperature == reference) return 0.0;
    const double lowerBound = std::min(temperature, reference);
    const double upperBound = std::max(temperature, reference);
    double integral = 0.0;
    double covered = lowerBound;

    for (std::size_t i = 1; i < table.size; ++i)
    {
        const IdealCaloricNode& lower = table.nodes[i - 1];
        const IdealCaloricNode& upper = table.nodes[i];
        const double start = std::max(lowerBound, lower.temperature);
        const double end = std::min(upperBound, upper.temperature);
        if (end <= start) continue;
        const double segment = integrateCpOverTemperature(
            lower, upper, start, end);
        if (!std::isfinite(segment))
            return std::numeric_limits<double>::quiet_NaN();
        integral += segment;
        covered = end;
    }
    if (std::abs(covered - upperBound)
        > 8.0*std::numeric_limits<double>::epsilon()
            *std::max(std::abs(upperBound), 1.0))
        return std::numeric_limits<double>::quiet_NaN();
    return temperature > reference ? integral : -integral;
}

bool validInterval(
    const IdealCaloricNode& lower,
    const IdealCaloricNode& upper)
{
    const double dt = upper.temperature - lower.temperature;
    if (!(std::isfinite(dt) && dt > 0.0) ||
        !(upper.sensibleEnthalpy > lower.sensibleEnthalpy))
    {
        return false;
    }

    // h(u)=a*u^3+b*u^2+c*u+d, so cp=dh/dT is a quadratic.
    // Checking its interior extremum proves positivity over the full segment.
    const double a =
        2.0*lower.sensibleEnthalpy
      - 2.0*upper.sensibleEnthalpy
      + dt*(lower.cp + upper.cp);
    const double b =
       -3.0*lower.sensibleEnthalpy
      + 3.0*upper.sensibleEnthalpy
      - dt*(2.0*lower.cp + upper.cp);

    if (a != 0.0)
    {
        const double extremum = -b/(3.0*a);
        if (extremum > 0.0 && extremum < 1.0)
        {
            const double temperature =
                lower.temperature + extremum*dt;
            const auto properties = interpolate(lower, upper, temperature);
            if (!(std::isfinite(properties.cp) && properties.cp > 0.0))
            {
                return false;
            }
        }
    }
    return true;
}

bool validateTable(const IdealCaloricTable& table)
{
    if (table.nodes == nullptr || table.size < 2) return false;
    bool foundReference = false;
    for (std::size_t i = 0; i < table.size; ++i)
    {
        const IdealCaloricNode& node = table.nodes[i];
        if (!std::isfinite(node.temperature) || node.temperature <= 0.0 ||
            !std::isfinite(node.sensibleEnthalpy) ||
            !std::isfinite(node.cp) || node.cp <= 0.0)
        {
            return false;
        }
        if (node.temperature == TabulatedIdealCaloric::referenceTemperature)
        {
            if (node.sensibleEnthalpy != 0.0) return false;
            foundReference = true;
        }
        if (i > 0 && !validInterval(table.nodes[i - 1], node)) return false;
    }
    return foundReference;
}

struct TableValidation
{
    bool valid[sizeof(tables)/sizeof(tables[0])];

    TableValidation()
    {
        for (std::size_t i = 0; i < sizeof(tables)/sizeof(tables[0]); ++i)
        {
            valid[i] = validateTable(tables[i]);
        }
    }
};

bool tableIsValid(const IdealCaloricTable& table)
{
    static const TableValidation validation;
    const std::size_t index = static_cast<std::size_t>(&table - tables);
    return index < sizeof(tables)/sizeof(tables[0]) && validation.valid[index];
}

} // End anonymous namespace

const char* tabulatedIdealCaloricStatusName(
    TabulatedIdealCaloricStatus status)
{
    switch (status)
    {
        case TabulatedIdealCaloricStatus::SUCCESS: return "SUCCESS";
        case TabulatedIdealCaloricStatus::UNKNOWN_SPECIES:
            return "UNKNOWN_SPECIES";
        case TabulatedIdealCaloricStatus::INVALID_TEMPERATURE:
            return "INVALID_TEMPERATURE";
        case TabulatedIdealCaloricStatus::OUT_OF_RANGE:
            return "OUT_OF_RANGE";
        case TabulatedIdealCaloricStatus::INVALID_FROZEN_TABLE:
            return "INVALID_FROZEN_TABLE";
        case TabulatedIdealCaloricStatus::NONFINITE_PROPERTY:
            return "NONFINITE_PROPERTY";
    }
    return "UNKNOWN_STATUS";
}

bool TabulatedIdealCaloric::has(const std::string& species) const
{
    return findTable(species) != nullptr;
}

bool TabulatedIdealCaloric::temperatureRange(
    const std::string& species,
    double& minimum,
    double& maximum) const
{
    const IdealCaloricTable* table = findTable(species);
    if (table == nullptr) return false;
    minimum = table->nodes[0].temperature;
    maximum = table->nodes[table->size - 1].temperature;
    return true;
}

std::vector<TabulatedIdealCaloricTableData>
TabulatedIdealCaloric::data() const
{
    std::vector<TabulatedIdealCaloricTableData> result;
    result.reserve(sizeof(tables)/sizeof(tables[0]));
    for (const IdealCaloricTable& table : tables)
    {
        TabulatedIdealCaloricTableData copy;
        copy.species = table.species;
        copy.nodes.assign(table.nodes, table.nodes + table.size);
        result.push_back(std::move(copy));
    }
    return result;
}

std::string TabulatedIdealCaloric::provenanceManifest(
    const std::vector<TabulatedIdealCaloricTableData>& tableData)
{
    std::string manifest = "schema=rfam.coolprop72-ideal-caloric/v1;";
    thermoProvenance::appendDouble(
        manifest, "referenceTemperature", referenceTemperature);
    thermoProvenance::appendSize(
        manifest, "table.count", tableData.size());
    for (std::size_t tableIndex = 0;
         tableIndex < tableData.size(); ++tableIndex)
    {
        const TabulatedIdealCaloricTableData& table = tableData[tableIndex];
        const std::string prefix =
            "table[" + std::to_string(tableIndex) + "].";
        thermoProvenance::appendString(
            manifest, prefix + "species", table.species);
        thermoProvenance::appendSize(
            manifest, prefix + "node.count", table.nodes.size());
        for (std::size_t nodeIndex = 0;
             nodeIndex < table.nodes.size(); ++nodeIndex)
        {
            const TabulatedIdealCaloricNode& node = table.nodes[nodeIndex];
            const std::string nodePrefix = prefix + "node["
                + std::to_string(nodeIndex) + "].";
            thermoProvenance::appendDouble(
                manifest, nodePrefix + "temperature", node.temperature);
            thermoProvenance::appendDouble(
                manifest,
                nodePrefix + "sensibleEnthalpy",
                node.sensibleEnthalpy);
            thermoProvenance::appendDouble(
                manifest, nodePrefix + "cp", node.cp);
        }
    }
    return manifest;
}

std::string TabulatedIdealCaloric::provenanceSha256(
    const std::vector<TabulatedIdealCaloricTableData>& tableData)
{
    return thermoProvenance::sha256Hex(provenanceManifest(tableData));
}

std::string TabulatedIdealCaloric::provenanceManifest() const
{
    return provenanceManifest(data());
}

std::string TabulatedIdealCaloric::provenanceSha256() const
{
    return provenanceSha256(data());
}

TabulatedIdealCaloricResult TabulatedIdealCaloric::evaluate(
    const std::string& species,
    double temperature) const
{
    TabulatedIdealCaloricResult result;
    result.species = species;

    const IdealCaloricTable* table = findTable(species);
    if (table == nullptr) return result;

    result.supportedTemperatureMinimum = table->nodes[0].temperature;
    result.supportedTemperatureMaximum =
        table->nodes[table->size - 1].temperature;

    if (!tableIsValid(*table))
    {
        result.status = TabulatedIdealCaloricStatus::INVALID_FROZEN_TABLE;
        return result;
    }
    if (!std::isfinite(temperature) || temperature <= 0.0)
    {
        result.status = TabulatedIdealCaloricStatus::INVALID_TEMPERATURE;
        return result;
    }
    if (temperature < result.supportedTemperatureMinimum ||
        temperature > result.supportedTemperatureMaximum)
    {
        result.status = TabulatedIdealCaloricStatus::OUT_OF_RANGE;
        return result;
    }

    const IdealCaloricNode* begin = table->nodes;
    const IdealCaloricNode* end = begin + table->size;
    const IdealCaloricNode* upper = std::upper_bound(
        begin, end, temperature,
        [](double value, const IdealCaloricNode& node)
        {
            return value < node.temperature;
        });
    if (upper == end) upper = end - 1;
    const IdealCaloricNode* lower =
        upper == begin ? begin : upper - 1;
    if (lower == upper) upper = begin + 1;

    TabulatedIdealCaloricProperties properties =
        interpolate(*lower, *upper, temperature);
    properties.sensibleEntropy = sensibleEntropy(*table, temperature);
    if (!std::isfinite(properties.sensibleEnthalpy) ||
        !std::isfinite(properties.cp) || properties.cp <= 0.0 ||
        !std::isfinite(properties.sensibleEntropy))
    {
        result.status = TabulatedIdealCaloricStatus::NONFINITE_PROPERTY;
        return result;
    }

    result.status = TabulatedIdealCaloricStatus::SUCCESS;
    result.state = properties;
    return result;
}

const TabulatedIdealCaloric& tabulatedIdealCaloric()
{
    static const TabulatedIdealCaloric instance;
    return instance;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
