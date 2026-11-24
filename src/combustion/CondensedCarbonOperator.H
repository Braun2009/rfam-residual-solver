#ifndef CondensedCarbonOperator_H
#define CondensedCarbonOperator_H

#include "CondensedCarbonModel.H"
#include "Operator.H"
#include "SectorMeshMetrics.H"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace Foam
{

std::string condensedCarbonSectorMeshFingerprint(
    const SectorMeshMetrics& mesh);

struct CondensedCarbonInventoryAudit
{
    double fluidMass = 0.0;
    double sootMass = 0.0;
    double wallCarbonMass = 0.0;
    double materialMass = 0.0;
    std::array<double,4> fluidElementMoles = {0.0,0.0,0.0,0.0};
    std::array<double,4> totalElementMoles = {0.0,0.0,0.0,0.0};
    double fluidResolvedEnergy = 0.0;
    double formationEnergyLedger = 0.0;
    double sootDiagnosticKineticEnergy = 0.0;
    double modeledEnergyLedger = 0.0;
    double sootBoundaryMassNetIntoDomain = 0.0;
    double sootBoundaryNumberNetIntoDomain = 0.0;
    bool physicalEnergyComplete = false;
    bool wallDepositionEnergyModeled = false;
};

class SectorCondensedCarbonOperator : public Operator
{
    const SectorMeshMetrics& mesh_;
    const SpeciesRegistry& registry_;
    CondensedCarbonOptions options_;
    const CondensedCarbonKinetics* kinetics_;
    const WallCarbonDepositionModel* wallModel_;
    std::string kineticsModelId_;
    std::string wallModelId_;
    std::string meshFingerprint_;
    std::string modelFingerprint_;

    mutable std::vector<CondensedCarbonState> state_;
    mutable std::vector<WallCarbonInventory> wall_;
    mutable bool initialized_ = false;
    mutable double sootBoundaryMassNetIntoDomain_ = 0.0;
    mutable double sootBoundaryNumberNetIntoDomain_ = 0.0;
    mutable SootApplicabilityReport lastApplicability_;

    mutable std::vector<CondensedCarbonState> stateStepStart_;
    mutable std::vector<WallCarbonInventory> wallStepStart_;
    mutable bool initializedStepStart_ = false;
    mutable double boundaryMassStepStart_ = 0.0;
    mutable double boundaryNumberStepStart_ = 0.0;
    mutable SootApplicabilityReport applicabilityStepStart_;
    mutable bool stepTransactionActive_ = false;

    int index(int i,int j,int k) const noexcept;
    int wallIndex(int i,int k) const noexcept;
    void initialize() const;
    void validateDomainLayout(const Domain& domain) const;
    void applyFieldBoundaries(
        const Domain& domain,
        std::vector<CondensedCarbonState>& field) const;
    double scalarOutwardFlux(
        const Domain& domain,
        const std::vector<CondensedCarbonState>& field,
        int i,int j,int k,int ni,int nj,int nk,
        const FaceAreaVector& area,
        bool numberField) const;
    double particleDiameter(const CondensedCarbonState& state) const;
    SootApplicabilityReport applicability(
        const Domain& domain,
        int i,int j,int k,
        const ConservedVars& fluid,
        const CondensedCarbonState& carbon) const;

public:
    SectorCondensedCarbonOperator(
        const SectorMeshMetrics& mesh,
        const SpeciesRegistry& registry,
        CondensedCarbonOptions options,
        const CondensedCarbonKinetics* kinetics = nullptr,
        const WallCarbonDepositionModel* wallModel = nullptr);

    bool ownsPersistentState() const override { return true; }
    void beginStepTransaction(const Domain& domain) const override;
    void acceptStepTransaction() const noexcept override;
    void rollbackStepTransaction() const noexcept override;
    void averageWithStepStartState() const override;
    void apply(Domain& domain,double dt) const override;
    double maxStableDt(const Domain& domain,double cfl) const override;
    std::string name() const override
    { return "Sector condensed-carbon v1"; }

    bool meshMatches(const Domain& domain) const noexcept;
    bool stepTransactionActive() const noexcept
    { return stepTransactionActive_; }
    bool persistentStateInitialized() const noexcept { return initialized_; }
    bool numberDensityEnabled() const noexcept
    { return options_.momentModel==CondensedCarbonMomentModel::TWO_MOMENT; }
    CondensedCarbonMomentModel momentModel() const noexcept
    { return options_.momentModel; }
    const CondensedCarbonOptions& options() const noexcept { return options_; }
    const std::string& modelFingerprint() const noexcept
    { return modelFingerprint_; }
    const std::string& meshFingerprint() const noexcept
    { return meshFingerprint_; }
    const std::string& kineticsModelId() const noexcept
    { return kineticsModelId_; }
    const std::string& wallModelId() const noexcept
    { return wallModelId_; }

    double rhoSoot(int i,int j,int k) const;
    double numberDensity(int i,int j,int k) const;
    double wallCarbonMassPerArea(int i,int k) const;
    const std::vector<CondensedCarbonState>& stateField() const;
    const std::vector<WallCarbonInventory>& wallField() const;
    double sootBoundaryMassNetIntoDomain() const noexcept
    { return sootBoundaryMassNetIntoDomain_; }
    double sootBoundaryNumberNetIntoDomain() const noexcept
    { return sootBoundaryNumberNetIntoDomain_; }
    const SootApplicabilityReport& lastApplicabilityReport() const noexcept
    { return lastApplicability_; }

    bool fieldsAreRestorable(
        const std::vector<CondensedCarbonState>& state,
        const std::vector<WallCarbonInventory>& wall,
        double boundaryMass,
        double boundaryNumber) const noexcept;
    void restoreFields(
        std::vector<CondensedCarbonState> state,
        std::vector<WallCarbonInventory> wall,
        double boundaryMass,
        double boundaryNumber);

    CondensedCarbonInventoryAudit audit(const Domain& domain) const;
};

} // End namespace Foam

#endif
