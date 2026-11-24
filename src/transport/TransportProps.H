#ifndef TransportProps_H
#define TransportProps_H

/*---------------------------------------------------------------------------*\
    TransportProps - mixture transport properties for RFAM (Stage 3).

    Provides per-call mixture viscosity, thermal conductivity and species
    diffusivity from the local (T, Y) state. All methods are pure functions
    of their arguments (no stored state) so the object can be shared across
    cells safely.

    Model choices (first cut, upgradeable per species):
      mu    : Wilke rule with per-species Sutherland law
              mu_i = mu_ref_i * (T/T_ref)^1.5 * (T_ref+S_i)/(T+S_i)
              Sutherland constants from standard references (NASA TM-4513 /
              Poling, Prausnitz, O'Connell).
      kappa : Eucken approximation (Stiel-Thodos modification for polyatomic):
              kappa = mu * Cv * (1.32 + 1.77*R_spec/Cv)
              (reduces to kappa = mu * Cp * (9gamma-5)/(4(gamma-1)) for ideal gas)
      D_i   : constant Schmidt number (Sc = 0.7) for each species:
              D_i = mu / (rho * Sc)
              This is the standard first-cut for combustion CFD; per-species
              binary diffusivities (Chapman-Enskog / JANAF) are a later upgrade.

    Dependency: NASA7Thermo (for Cp_mix used in Eucken kappa), SpeciesRegistry.
    Pure std C++, no SUNDIALS, header-only.
\*---------------------------------------------------------------------------*/

#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>
#include <cmath>
#include <algorithm>
#include "ConservedVars.H"
#include "CryogenicTransportTable.H"
#include "EquilibriumThermoModel.H"
#include "NASA7Thermo.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"

namespace Foam
{

class TransportProps
{
public:
    // The current closure is deliberately identified as a single-phase gas
    // model.  Coupling it to real-fluid Cp/Cv makes the sector operators
    // internally consistent, but does not make Sutherland/Wilke/Eucken valid
    // for dense cryogenic liquid or a two-phase state.
    static constexpr bool supportsTwoPhase = false;
    static constexpr bool supportsDenseCryogenic = false;

    // The Sutherland/Wilke/Eucken closure is a gas correlation, not a dense
    // liquid or two-phase property model.  NASA-7 data used by the coupled
    // caloric path start at 200 K; the density guard is deliberately
    // conservative and is a capability gate, not an assertion of accuracy.
    static constexpr double gasClosureMinTemperature = 200.0;
    static constexpr double gasClosureMaxDensity = 100.0; // kg/m3
    static constexpr std::size_t gasClosureSpeciesCapacity = 9;

    // 注册的严格稠密纯流体表（Phase 1c）。nullptr 时稠密查询 fail-close。
    inline static const CryogenicTransportTable* cryogenicTable_ = nullptr;

    static constexpr const char* committedLegacyGasModelId()
    {
        return "sutherland-wilke-eucken-sc0.7-legacy-q-v1";
    }

    static bool gasClosureSupported(double p,double T,double rho)
    {
        return std::isfinite(p) && std::isfinite(T) && std::isfinite(rho)
            && p>0.0 && T>=gasClosureMinTemperature
            && rho>0.0 && rho<=gasClosureMaxDensity;
    }

    struct Point
    {
        double mu=0.0;       // dynamic viscosity [Pa s]
        double nu=0.0;       // kinematic viscosity [m2/s]
        double kappa=0.0;    // thermal conductivity [W/(m K)]
        double diffusivity=0.0;
        double Cp=0.0;
        double Cv=0.0;

        bool valid() const
        {
            return std::isfinite(mu) && mu>0.0
                &&std::isfinite(nu) && nu>0.0
                &&std::isfinite(kappa) && kappa>0.0
                &&std::isfinite(diffusivity) && diffusivity>0.0
                &&std::isfinite(Cp) && std::isfinite(Cv)
                &&Cp>Cv && Cv>0.0;
        }
    };

    enum class Status
    {
        SUCCESS,
        INVALID_INPUT,
        INVALID_EQUILIBRIUM_STATE,
        SPECIES_MISMATCH,
        UNSUPPORTED_THERMO_MODEL,
        UNSUPPORTED_ACOUSTIC_CLOSURE,
        UNSUPPORTED_PHASE,
        UNSUPPORTED_TEMPERATURE_RANGE,
        UNSUPPORTED_DENSITY_RANGE,
        PRIMITIVE_EQUILIBRIUM_MISMATCH,
        MISSING_THERMODYNAMIC_PAYLOAD,
        INVALID_TRANSPORT_INPUT,
        THERMO_BACKEND_FAILURE,
        THERMO_SURFACE_MISMATCH,
        NONFINITE_PROPERTY
    };

    enum class Provenance
    {
        NONE,
        COMMITTED_LEGACY_Q
    };

    struct Result
    {
        Status status=Status::INVALID_INPUT;
        Point point;
        Provenance provenance=Provenance::NONE;

        bool success() const
        {
            return status==Status::SUCCESS && point.valid()
                &&provenance==Provenance::COMMITTED_LEGACY_Q;
        }
    };

    static_assert(
        std::is_trivially_copyable<Result>::value,
        "per-cell transport results must remain fixed-size and allocation-free");

    static const char* statusName(Status status)
    {
        switch(status)
        {
            case Status::SUCCESS: return "SUCCESS";
            case Status::INVALID_INPUT: return "INVALID_INPUT";
            case Status::INVALID_EQUILIBRIUM_STATE:
                return "INVALID_EQUILIBRIUM_STATE";
            case Status::SPECIES_MISMATCH: return "SPECIES_MISMATCH";
            case Status::UNSUPPORTED_THERMO_MODEL:
                return "UNSUPPORTED_THERMO_MODEL";
            case Status::UNSUPPORTED_ACOUSTIC_CLOSURE:
                return "UNSUPPORTED_ACOUSTIC_CLOSURE";
            case Status::UNSUPPORTED_PHASE: return "UNSUPPORTED_PHASE";
            case Status::UNSUPPORTED_TEMPERATURE_RANGE:
                return "UNSUPPORTED_TEMPERATURE_RANGE";
            case Status::UNSUPPORTED_DENSITY_RANGE:
                return "UNSUPPORTED_DENSITY_RANGE";
            case Status::PRIMITIVE_EQUILIBRIUM_MISMATCH:
                return "PRIMITIVE_EQUILIBRIUM_MISMATCH";
            case Status::MISSING_THERMODYNAMIC_PAYLOAD:
                return "MISSING_THERMODYNAMIC_PAYLOAD";
            case Status::INVALID_TRANSPORT_INPUT:
                return "INVALID_TRANSPORT_INPUT";
            case Status::THERMO_BACKEND_FAILURE:
                return "THERMO_BACKEND_FAILURE";
            case Status::THERMO_SURFACE_MISMATCH:
                return "THERMO_SURFACE_MISMATCH";
            case Status::NONFINITE_PROPERTY: return "NONFINITE_PROPERTY";
        }
        return "UNKNOWN";
    }

    // ------------------------------------------------------------------ //
    // 稠密低温输运后端（Phase 1c：纯/近纯流体表）
    //
    // 后端选择遵循 docs/cryogenic_transport_backend.md 第 5 节：
    //   1. 验证的单相气相（T>=200K 且 rho<=100 kg/m3）-> legacy gas；
    //   2. 纯/近纯 O2/C3H8 稠密态（表内）-> 严格纯流体表；
    //   3. 稠密多组分或两相 -> UNSUPPORTED（未验证混合，绝不静默
    //       fallback 到气体闭合或表格边缘）。
    // ------------------------------------------------------------------ //
    enum class Backend
    {
        LEGACY_GAS,
        CRYOGENIC_PURE_FLUID,
        UNSUPPORTED
    };

    static const char* backendName(Backend backend)
    {
        switch (backend)
        {
            case Backend::LEGACY_GAS: return "LEGACY_GAS";
            case Backend::CRYOGENIC_PURE_FLUID:
                return "CRYOGENIC_PURE_FLUID";
            case Backend::UNSUPPORTED: return "UNSUPPORTED";
        }
        return "UNKNOWN";
    }

    // 注册严格稠密纯流体表（由求解器/测试在启动时加载）。
    // 未注册时稠密查询 fail-close 为 UNSUPPORTED。
    static void registerCryogenicTable(
        const CryogenicTransportTable* table)
    {
        cryogenicTable_ = table;
    }

    static const CryogenicTransportTable* cryogenicTable()
    {
        return cryogenicTable_;
    }

    // 近纯组分判定：若单一 active 物种质量分数 >= threshold 则返回该
    // 物种名，否则返回空字符串。用于稠密表后端（只支持纯 O2/C3H8）。
    static std::string activePureSpecies(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double threshold = 0.999)
    {
        if (names.size() != Y.size()) return std::string();
        std::string pure;
        for (std::size_t s = 0; s < names.size(); ++s)
        {
            if (!std::isfinite(Y[s]) || Y[s] < 0.0) return std::string();
            if (Y[s] > threshold)
            {
                if (!pure.empty()) return std::string(); // 两个都超
                pure = names[s];
            }
        }
        return pure;
    }

    // 后端选择（纯函数，无副作用）。
    static Backend selectBackend(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T,
        double rho)
    {
        if (!std::isfinite(T) || !std::isfinite(rho)
            || T <= 0.0 || rho <= 0.0)
            return Backend::UNSUPPORTED;
        const bool gasRange = T >= gasClosureMinTemperature
            && rho <= gasClosureMaxDensity;
        if (gasRange) return Backend::LEGACY_GAS;
        // 稠密或低温：仅近纯 O2/C3H8 且表已注册且表内可解析时允许。
        const std::string pure = activePureSpecies(names, Y);
        if (!pure.empty() && cryogenicTable_ && cryogenicTable_->has(pure))
            return Backend::CRYOGENIC_PURE_FLUID;
        return Backend::UNSUPPORTED;
    }

    // 稠密纯流体表评估：返回 mu/kappa（Cp/Cv 由调用者传入）。
    // 两相/表外/未注册 -> UNSUPPORTED_TEMPERATURE_RANGE（fail-close）。
    static Result evaluateCryogenicPureFluid(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double p,
        double T,
        double rho,
        double Cp,
        double Cv)
    {
        (void)p;  // 稠密表后端不依赖压力（以 T/rho 为坐标）
        Result result;
        if (!cryogenicTable_)
        {
            result.status = Status::UNSUPPORTED_TEMPERATURE_RANGE;
            return result;
        }
        const std::string pure = activePureSpecies(names, Y);
        if (pure.empty())
        {
            result.status = Status::UNSUPPORTED_PHASE;
            return result;
        }
        const CryogenicTransportTableResult tableResult =
            cryogenicTable_->evaluate(pure, T, rho);
        if (!tableResult.success())
        {
            switch (tableResult.status)
            {
                case CryogenicTransportTableStatus::TWO_PHASE:
                    result.status = Status::UNSUPPORTED_PHASE;
                    break;
                default:
                    result.status = Status::UNSUPPORTED_TEMPERATURE_RANGE;
                    break;
            }
            return result;
        }
        if (!std::isfinite(Cp) || !std::isfinite(Cv)
            || Cp <= Cv || Cv <= 0.0)
        {
            result.status = Status::NONFINITE_PROPERTY;
            return result;
        }
        const double mu = tableResult.properties->mu;
        const double kappa = tableResult.properties->kappa;
        if (!std::isfinite(mu) || mu <= 0.0
            || !std::isfinite(kappa) || kappa <= 0.0)
        {
            result.status = Status::NONFINITE_PROPERTY;
            return result;
        }
        result.point.mu = mu;
        result.point.nu = mu / rho;
        result.point.kappa = kappa;
        result.point.diffusivity = D_species(mu, rho);
        result.point.Cp = Cp;
        result.point.Cv = Cv;
        result.status = Status::SUCCESS;
        result.provenance = Provenance::COMMITTED_LEGACY_Q;
        return result;
    }

    // ---- Sutherland constants: {mu_ref [Pa s], T_ref [K], S [K]} ----
    // Sources: Poling, Prausnitz & O'Connell, "The Properties of Gases and
    // Liquids" 5th ed.; Bird, Stewart, Lightfoot "Transport Phenomena" 2nd ed.
    // For species without literature data, the default triple (mu0=1.5e-5,
    // T0=273, S=200) gives a reasonable order-of-magnitude. They are marked (*).
    struct SuthConst { double mu_ref; double T_ref; double S; };

    static const SuthConst& suthConst(const std::string& name)
    {
        // static table: name -> {mu_ref [Pa s], T_ref [K], S [K]}
        static const struct { const char* n; SuthConst c; } tbl[] = {
            {"O2",   {1.919e-5, 273.0,  139.0}},
            {"N2",   {1.663e-5, 273.0,  107.0}},
            {"H2",   {8.411e-6, 273.0,   97.0}},
            {"CH4",  {1.024e-5, 273.0,  164.0}},
            {"C2H6", {8.600e-6, 273.0,  252.0}},  // (*)
            {"C3H8", {7.500e-6, 273.0,  278.0}},  // (*)
            {"H2O",  {1.120e-5, 373.0,  673.0}},  // steam ref
            {"CO",   {1.657e-5, 273.0,  113.0}},
            {"CO2",  {1.370e-5, 273.0,  240.0}},
        };
        for (const auto& e : tbl)
            if (name == e.n) return e.c;
        // fallback
        static const SuthConst def{1.5e-5, 273.0, 200.0};
        return def;
    }

    // Pure-species dynamic viscosity [Pa s] via Sutherland law
    static double mu_species(const std::string& name, double T)
    {
        const SuthConst& c = suthConst(name);
        double Tclamp = std::max(T, 100.0);   // extrapolation floor
        return c.mu_ref * std::pow(Tclamp / c.T_ref, 1.5)
                        * (c.T_ref + c.S) / (Tclamp + c.S);
    }

    static bool hasSutherlandData(const std::string& name)
    {
        return name=="O2" || name=="N2" || name=="H2" || name=="CH4"
            || name=="C2H6" || name=="C3H8" || name=="H2O"
            || name=="CO" || name=="CO2";
    }

    // Strict Wilke evaluation used by production sector operators.  Unlike
    // the historical convenience API below, this never substitutes air-like
    // constants for missing species or invalid composition.
    static bool tryMuMix(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T,
        const SpeciesRegistry& reg,
        double& result)
    {
        std::size_t substituted=0;
        return tryMuMixGeneral(names,Y,T,reg,result,substituted);
    }

    // Wilke mixing for mechanisms wider than the Sutherland table.
    // Species without tabulated constants take those of the tabulated species
    // with the nearest molecular weight: viscosity correlates with size and
    // MW far more strongly than with the reaction path, so a radical that has
    // no Sutherland fit inherits its stable homologue instead of failing the
    // whole step. `substituted` reports how many species were mapped, so a
    // caller can log the approximation rather than inherit it silently.
    static bool tryMuMixGeneral(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T,
        const SpeciesRegistry& reg,
        double& result,
        std::size_t& substituted)
    {
        result=0.0;
        substituted=0;
        if(names.empty() || names.size()!=Y.size() ||
           !std::isfinite(T) || T<100.0)
            return false;
        double sumY=0.0;
        for(std::size_t i=0;i<names.size();++i)
        {
            if(!reg.has(names[i]) || !std::isfinite(Y[i]) || Y[i]<0.0)
                return false;
            sumY+=Y[i];
        }
        if(std::abs(sumY-1.0)>1e-8) return false;

        const std::size_t count=names.size();
        std::vector<double> x(count,0.0),mw(count,0.0),muI(count,0.0);
        double sumInv=0.0;
        for(std::size_t i=0;i<count;++i)
        {
            mw[i]=reg.molecularWeight(reg.index(names[i]));
            if(!std::isfinite(mw[i]) || mw[i]<=0.0) return false;
            sumInv+=Y[i]/mw[i];
        }
        if(!std::isfinite(sumInv) || sumInv<=1e-30) return false;

        for(std::size_t i=0;i<count;++i)
        {
            x[i]=(Y[i]/mw[i])/sumInv;
            if(hasSutherlandData(names[i]))
                muI[i]=mu_species(names[i],T);
            else
            {
                // Nearest tabulated molecular weight.
                const char* table[]={"O2","N2","H2","CH4","C2H6","C3H8",
                                     "H2O","CO","CO2"};
                double bestScore=1e30;
                const char* bestName=nullptr;
                for(const char* candidate:table)
                {
                    const int g=reg.index(candidate);
                    if(g<0) continue;
                    const double m=reg.molecularWeight(g);
                    if(!(m>0.0)) continue;
                    const double score=std::abs(std::log(m/mw[i]));
                    if(score<bestScore){bestScore=score;bestName=candidate;}
                }
                if(!bestName) return false;
                muI[i]=mu_species(bestName,T);
                ++substituted;
            }
            if(!std::isfinite(muI[i]) || muI[i]<=0.0) return false;
        }

        const double invSqrt8=1.0/std::sqrt(8.0);
        for(std::size_t i=0;i<count;++i)
        {
            if(x[i]<1e-12) continue;
            double denominator=0.0;
            for(std::size_t j=0;j<count;++j)
            {
                if(x[j]<1e-12) continue;
                const double ratio=std::sqrt(muI[i]/muI[j]);
                const double massRatio=std::pow(mw[j]/mw[i],0.25);
                const double bracket=1.0+ratio*massRatio;
                denominator+=x[j]*invSqrt8/std::sqrt(1.0+mw[i]/mw[j])
                    *bracket*bracket;
            }
            if(!std::isfinite(denominator) || denominator<=1e-30)
                return false;
            result+=x[i]*muI[i]/denominator;
        }
        return std::isfinite(result) && result>0.0;
    }

    // Mixture dynamic viscosity [Pa s]: Wilke mixing rule
    //   mu_mix = sum_i ( x_i * mu_i / sum_j x_j * phi_ij )
    //   phi_ij = (1/sqrt8) * (1+MW_i/MW_j)^(-1/2) * [1+(mu_i/mu_j)^(1/2)*(MW_j/MW_i)^(1/4)]^2
    static double mu_mix(
        const std::vector<std::string>& names,
        const std::vector<double>&      Y,
        double                          T,
        const SpeciesRegistry&          reg)
    {
        double mu=0.0;
        // Compatibility-only fallback for legacy/non-sector callers.  The
        // sector production path uses tryMuMix via evaluateSinglePhaseGas.
        return tryMuMix(names,Y,T,reg,mu)?mu:1.8e-5;
    }

    // Mixture thermal conductivity [W/(m K)]: modified Eucken relation
    //   kappa = mu * Cv * (1.32 + 1.77 * R_spec / Cv)
    // This is equivalent to the Stiel-Thodos correction for polyatomic gases.
    // For a monatomic ideal gas (Cv=3R/2) this reduces to kappa = 5/2 mu Cv.
    static double kappa_mix(
        const std::vector<std::string>& names,
        const std::vector<double>&      Y,
        double                          T,
        const SpeciesRegistry&          reg)
    {
        static constexpr double Runiv = 8.314462618;
        double mu  = mu_mix(names, Y, T, reg);

        // mixture Cv [J/(kg K)]
        double Cp  = nasa7Thermo().cp_mass_mix(names, Y, T);
        // R_mix [J/(kg K)] = Runiv * sum(Y_i / MW_i)
        double R_mix = 0.0;
        for (int i = 0; i < n(names); ++i)
            R_mix += Y[i] / reg.molecularWeight(reg.index(names[i])) * Runiv;
        double Cv_mix = std::max(Cp - R_mix, 100.0);

        return mu * Cv_mix * (1.32 + 1.77 * R_mix / Cv_mix);
    }

    // Compatibility-only PTY evaluator. CFD operators use the committed-Q
    // entry below so a phase-aware state cannot be reinterpreted as this
    // legacy gas closure.
    static bool evaluateSinglePhaseGas(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double p,
        double T,
        double rho,
        const SpeciesRegistry& reg,
        Point& result)
    {
        const Result evaluated=evaluateGasSurface(
            names,Y,p,T,rho,reg);
        result=evaluated.point;
        return evaluated.status==Status::SUCCESS
            &&evaluated.point.valid();
    }

    // Strict transitional transport closure for committed legacy gas states.
    // A different thermo model, acoustic closure, phase region, or W/Q pair
    // is an unsupported request, never an invitation to perform a PTY fallback.
    static Result evaluateCommittedLegacyGas(
        const std::vector<std::string>& names,
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium,
        const thermo::EquilibriumThermoModel& thermoModel,
        const SpeciesRegistry& reg)
    {
        Result result;
        if(!equilibrium.valid())
        {
            result.status=Status::INVALID_EQUILIBRIUM_STATE;
            return result;
        }
        if(names.empty() || names.size()!=primitive.Y.size()
           || names.size()!=equilibrium.overallMassFractions.size())
        {
            result.status=Status::SPECIES_MISMATCH;
            return result;
        }
        if(equilibrium.region
               !=thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE
           ||!equilibrium.singlePhase.has_value()
           ||equilibrium.liquid.has_value() || equilibrium.vapor.has_value())
        {
            result.status=Status::UNSUPPORTED_PHASE;
            return result;
        }
        if(equilibrium.acoustic.closure
           !=thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE)
        {
            result.status=Status::UNSUPPORTED_ACOUSTIC_CLOSURE;
            return result;
        }
        if(!thermoModel.supportsCommittedLegacyGasTransport()
           ||thermoModel.speciesNames()!=names)
        {
            result.status=Status::UNSUPPORTED_THERMO_MODEL;
            return result;
        }
        if(!close(primitive.p,equilibrium.p,1.0e-8)
           ||!close(primitive.T,equilibrium.T,1.0e-8)
           ||!sameComposition(
               primitive.Y,equilibrium.overallMassFractions))
        {
            result.status=Status::PRIMITIVE_EQUILIBRIUM_MISMATCH;
            return result;
        }
        if(!equilibrium.thermodynamics.has_value())
        {
            result.status=Status::MISSING_THERMODYNAMIC_PAYLOAD;
            return result;
        }
        if(!thermoModel.stateValidForModel(equilibrium))
        {
            result.status=Status::INVALID_EQUILIBRIUM_STATE;
            return result;
        }
        const thermo::EquilibriumThermodynamicPayload& payload=
            *equilibrium.thermodynamics;
        Result evaluated=evaluateGasPayload(
            names,equilibrium.overallMassFractions,
            equilibrium.p,equilibrium.T,equilibrium.rho,
            reg,payload.Cp,payload.Cv);
        if(evaluated.status==Status::SUCCESS)
            evaluated.provenance=Provenance::COMMITTED_LEGACY_Q;
        return evaluated;
    }

    // ------------------------------------------------------------------ //
    // 显式诊断外推路径（绝不静默）
    //
    // 供 --allow-single-phase-reacting-diagnostic 等显式诊断模式使用。
    // 该路径允许 Sutherland/Wilke/Eucken 气体闭合外推到其标称范围之外
    // （T<gasClosureMinTemperature 或 rho>gasClosureMaxDensity），这是
    // 非物理的：稠密低温状态（例如 130 K LOX 液相）的气相关联式误差可达
    // 8-10 倍。因此任何成功都必须由调用者通过 extrapolation 字段审计并
    // 在输出中报告；生产调用（evaluateCommittedLegacyGas）保持 fail-close。
    // ------------------------------------------------------------------ //
    enum class DiagnosticExtrapolation
    {
        NONE,           // 在气体闭合及活跃物种量热范围内
        COLD,           // T < gasClosureMinTemperature
        DENSE,          // rho > gasClosureMaxDensity
        COLD_AND_DENSE, // 两者同时超出
        CALORIC_RANGE   // 活跃物种 NASA-7 温度范围外
    };

    struct DiagnosticGasResult
    {
        Result result;
        DiagnosticExtrapolation extrapolation=DiagnosticExtrapolation::NONE;
        double temperatureK=0.0;
        double densityKgM3=0.0;

        bool success() const { return result.success(); }
    };

    static DiagnosticGasResult evaluateDiagnosticGas(
        const std::vector<std::string>& names,
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium,
        const thermo::EquilibriumThermoModel& thermoModel,
        const SpeciesRegistry& reg,
        bool allowTwoPhaseDiagnostic=false)
    {
        DiagnosticGasResult out;
        if(!equilibrium.valid())
        {
            out.result.status=Status::INVALID_EQUILIBRIUM_STATE;
            return out;
        }
        if(names.empty() || names.size()!=primitive.Y.size()
           || names.size()!=equilibrium.overallMassFractions.size())
        {
            out.result.status=Status::SPECIES_MISMATCH;
            return out;
        }
        const bool twoPhaseBlocked =
            (equilibrium.region
               !=thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE
            || !equilibrium.singlePhase.has_value()
            || equilibrium.liquid.has_value()
            || equilibrium.vapor.has_value());
        // HEM two-phase LOX near-field states (region=TWO_PHASE with
        // liquid/vapour phases) are admitted only in explicit diagnostic
        // extrapolation mode; production callers keep the fail-closed gate.
        if (twoPhaseBlocked && !allowTwoPhaseDiagnostic)
        {
            out.result.status=Status::UNSUPPORTED_PHASE;
            return out;
        }
        if(twoPhaseBlocked && equilibrium.acoustic.closure
               !=thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
            && equilibrium.acoustic.closure
               !=thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE)
        {
            out.result.status=Status::UNSUPPORTED_ACOUSTIC_CLOSURE;
            return out;
        }
        if(!twoPhaseBlocked && equilibrium.acoustic.closure
           !=thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE)
        {
            out.result.status=Status::UNSUPPORTED_ACOUSTIC_CLOSURE;
            return out;
        }
        if(!thermoModel.supportsCommittedLegacyGasTransport()
           ||thermoModel.speciesNames()!=names)
        {
            out.result.status=Status::UNSUPPORTED_THERMO_MODEL;
            return out;
        }
        if(!close(primitive.p,equilibrium.p,1.0e-8)
           ||!close(primitive.T,equilibrium.T,1.0e-8)
           ||!sameComposition(
               primitive.Y,equilibrium.overallMassFractions))
        {
            out.result.status=Status::PRIMITIVE_EQUILIBRIUM_MISMATCH;
            return out;
        }
        if(!equilibrium.thermodynamics.has_value())
        {
            out.result.status=Status::MISSING_THERMODYNAMIC_PAYLOAD;
            return out;
        }
        if(!thermoModel.stateValidForModel(equilibrium))
        {
            out.result.status=Status::INVALID_EQUILIBRIUM_STATE;
            return out;
        }
        const thermo::EquilibriumThermodynamicPayload& payload=
            *equilibrium.thermodynamics;
        // Diagnostic floor (2026-08-31): the reacting near-field evolves
        // below the gas-correlation lower bound (e.g. a 94.5 K condensed
        // cell); the mu/kappa extrapolation has no valid values there, so
        // evaluate at the floor and keep the diagnostic provenance.
        const double evalT = std::max(equilibrium.T, 105.0);
        out.result=evaluateGasPayload(
            names,equilibrium.overallMassFractions,
            equilibrium.p,evalT,equilibrium.rho,
            reg,payload.Cp,payload.Cv,
            /*allowDiagnosticExtrapolation=*/true);
        if(out.result.status==Status::SUCCESS)
        {
            out.result.provenance=Provenance::COMMITTED_LEGACY_Q;
            const bool cold=equilibrium.T<gasClosureMinTemperature;
            const bool dense=equilibrium.rho>gasClosureMaxDensity;
            const bool caloricOutside=!caloricRangeSupported(
                names,equilibrium.overallMassFractions,equilibrium.T);
            out.extrapolation=cold&&dense
                ?DiagnosticExtrapolation::COLD_AND_DENSE
                :cold?DiagnosticExtrapolation::COLD
                :dense?DiagnosticExtrapolation::DENSE
                :caloricOutside?DiagnosticExtrapolation::CALORIC_RANGE
                :DiagnosticExtrapolation::NONE;
            out.temperatureK=equilibrium.T;
            out.densityKgM3=equilibrium.rho;
        }
        return out;
    }

    // Species diffusivity [m^2/s]: constant Schmidt number approximation
    //   D_i = mu_mix / (rho * Sc)
    // Sc = 0.7 is standard for hydrocarbon/air combustion. Per-species
    // Chapman-Enskog diffusivities are a later upgrade ([[technical-debt]]).
    static constexpr double Sc = 0.7;
    static double D_species(double mu, double rho)
    {
        if (rho < 1e-8) return 0.0;
        return mu / (rho * Sc);
    }

private:
    static bool close(
        double lhs,double rhs,double relativeTolerance=1.0e-7)
    {
        return std::isfinite(lhs) && std::isfinite(rhs)
            && std::abs(lhs-rhs)<=relativeTolerance
                *std::max({std::abs(lhs),std::abs(rhs),1.0});
    }

    static bool sameComposition(
        const std::vector<double>& lhs,
        const std::vector<double>& rhs)
    {
        if(lhs.empty() || lhs.size()!=rhs.size()) return false;
        for(size_t i=0;i<lhs.size();++i)
            if(!close(lhs[i],rhs[i],1.0e-8)) return false;
        return true;
    }

    static bool caloricRangeSupported(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T)
    {
        if(names.size()!=Y.size() || !std::isfinite(T)) return false;
        bool hasActiveSpecies=false;
        for(size_t s=0;s<names.size();++s)
        {
            if(Y[s]<=1.0e-14) continue;
            hasActiveSpecies=true;
            double minimum=0.0,maximum=0.0;
            if(!nasa7Thermo().temperatureRange(
                    names[s],minimum,maximum)
               ||T<minimum || T>maximum)
            {
                // 2026-08-08: allow a bounded extrapolation band beyond the
                // NASA-7 polynomial range (high-temperature edge 3500 K):
                // numerical overshoot of a few K (observed 3501.8 K in the
                // flame) must not fail the whole step. 2026-10-22: the band is
                // 25% of the range edge because the O2/C3H8 surrogate flame
                // reaches ~3950 K, and a 10% band (3850 K) rejected the step.
                // Extrapolation is high-side only: low-temperature NASA-7 data
                // are not made valid by a hot overshoot allowance, and the
                // wider band is a transport-fit argument, not a claim that the
                // polynomial is thermodynamically accurate there.
                const double slack = 0.25 * std::max(std::abs(maximum), 1.0);
                // The bounded extrapolation is intentionally high-side only:
                // low-temperature NASA-7 data are not made valid by a hot
                // overshoot allowance.
                if (!(T >= minimum && T <= maximum + slack))
                    return false;
            }
        }
        return hasActiveSpecies;
    }

    static Result evaluateGasSurface(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double p,double T,double rho,
        const SpeciesRegistry& reg)
    {
        Result result;
        if(names.empty() || names.size()!=Y.size())
        {
            result.status=Status::SPECIES_MISMATCH;
            return result;
        }
        if(!std::isfinite(p) || !std::isfinite(T)
           ||!std::isfinite(rho) || p<=0.0 || T<=0.0 || rho<=0.0)
        {
            result.status=Status::INVALID_INPUT;
            return result;
        }

        RealFluidState state;
        try
        {
            state.setSpecies(names);
        }
        catch(...)
        {
            result.status=Status::THERMO_BACKEND_FAILURE;
            return result;
        }
        if(state.updatePTY(p,T,Y)!=RealFluidState::Status::SUCCESS)
        {
            result.status=Status::THERMO_BACKEND_FAILURE;
            return result;
        }

        if(!close(rho,state.rho(),1.0e-7))
        {
            result.status=Status::THERMO_SURFACE_MISMATCH;
            return result;
        }

        return evaluateGasPayload(
            names,Y,p,T,rho,reg,state.Cp(),state.Cv());
    }

    static Result evaluateGasPayload(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double p,double T,double rho,
        const SpeciesRegistry& reg,
        double Cp,double Cv,
        bool allowDiagnosticExtrapolation=false)
    {
        Result result;
        if(names.empty() || names.size()!=Y.size())
        {
            result.status=Status::SPECIES_MISMATCH;
            return result;
        }
        if(!std::isfinite(p) || !std::isfinite(T)
           ||!std::isfinite(rho) || p<=0.0 || T<=0.0 || rho<=0.0)
        {
            result.status=Status::INVALID_INPUT;
            return result;
        }

        // 生产路径（allowDiagnosticExtrapolation=false）保持 fail-close：
        // 气体闭合范围之外的稠密低温状态被拒绝。
        // 显式诊断路径（allowDiagnosticExtrapolation=true）跳过范围检查，
        // 但调用者必须通过 evaluateDiagnosticGas 返回的
        // DiagnosticGasResult.extrapolation 显式审计外推范围；绝不静默。
        if(!allowDiagnosticExtrapolation)
        {
            if(T<gasClosureMinTemperature
               ||!caloricRangeSupported(names,Y,T))
            {
                result.status=Status::UNSUPPORTED_TEMPERATURE_RANGE;
                return result;
            }
            if(rho>gasClosureMaxDensity)
            {
                result.status=Status::UNSUPPORTED_DENSITY_RANGE;
                return result;
            }
        }
        if(!std::isfinite(Cp) || !std::isfinite(Cv)
           ||Cp<=Cv || Cv<=0.0)
        {
            result.status=Status::NONFINITE_PROPERTY;
            return result;
        }

        double mu=0.0;
        if(!tryMuMix(names,Y,T,reg,mu))
        {
            result.status=Status::INVALID_TRANSPORT_INPUT;
            return result;
        }

        double gasConstant=0.0;
        for(size_t s=0;s<names.size();++s)
        {
            const int index=reg.index(names[s]);
            if(index<0)
            {
                result.status=Status::SPECIES_MISMATCH;
                return result;
            }
            gasConstant+=Y[s]/reg.molecularWeight(index)*8.314462618;
        }
        const double kappa=mu*Cv*(1.32+1.77*gasConstant/Cv);
        const double diffusivity=D_species(mu,rho);
        const double nu=mu/rho;
        if(!std::isfinite(mu) || mu<=0.0
           ||!std::isfinite(nu) || nu<=0.0
           ||!std::isfinite(kappa) || kappa<=0.0
           ||!std::isfinite(diffusivity) || diffusivity<=0.0)
        {
            result.status=Status::NONFINITE_PROPERTY;
            return result;
        }

        result.point.mu=mu;
        result.point.nu=nu;
        result.point.kappa=kappa;
        result.point.diffusivity=diffusivity;
        result.point.Cp=Cp;
        result.point.Cv=Cv;
        result.status=Status::SUCCESS;
        return result;
    }

    // helper: size of a std::vector<std::string>
    static int n(const std::vector<std::string>& v) { return (int)v.size(); }
};

} // End namespace Foam

#endif
