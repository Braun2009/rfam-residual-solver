#!/usr/bin/env python3
"""Validate the frozen RFAM low-concentration regularization contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any


class ContractError(RuntimeError):
    pass


EXPECTED_CONTRACT_PATH = Path(
    "tests/reference/r2_low_concentration_regularization_contract_v1.json"
)
EXPECTED_CONTRACT_SHA256 = (
    "ae8d748d5cf1f3b510719222bf14e219109278fcc1711c6f97810086dc18a0ba"
)

EXPECTED_ARTIFACTS = [
    ("probe_source", "tools/r2_low_concentration_regularization_probe.cpp", 16071,
     "7727dc47f71dedc1bd973243b89ad668a1dad4ad2c2c30b388c54c3bd7710ae1"),
    ("arrhenius_rate", "src/combustion/FiniteRate/ArrheniusRate.H", 1814,
     "d768ea9ef577abd3851837e394e47e343b3b957d4cd7fc2cb7210737a125f20c"),
    ("mass_action_thermo", "src/combustion/FiniteRate/MassActionThermo.H", 4511,
     "7cc012ea2b0bf9e6cee067b40a7590e3e5f3a27819b04fecca3edb8f2fe94d7b"),
    ("nasa7_thermo", "src/thermo/NASA7Thermo.H", 14091,
     "40ca86c9b52ee7332299ae8401680cc8ac84ad0dd035aab1183d4df9b5e5fafc"),
    ("species_registry", "src/thermo/SpeciesRegistry.H", 5789,
     "22b3c0e96d5bff8b6ffa6b69ecd7568a4f87422f5188af85de31c50136b27ba7"),
    ("thermo_provenance_header", "src/thermo/ThermoDataProvenance.H", 660,
     "eaac29b0f0cea543b2a0ae20adb6e570faf16f849ede709ab0e6140f856af2a9"),
    ("real_fluid_cvode_wrapper",
     "src/combustion/FiniteRate/RealFluidCVODEWrapper.H", 36186,
     "88ddb4dad7132f5143cf7781a281a599a73d092d8a9321a9c39ad3a028f8b0f2"),
    ("rhs_composition_preparation",
     "src/combustion/FiniteRate/RhsCompositionPreparation.H", 3578,
     "ebf7511a90b6fc2a9745d2008da8442973b7fedde99edb54f580091a86fb806a"),
    ("chemistry_model", "src/combustion/FiniteRate/ChemistryModel.H", 2390,
     "f6fab066c321b3de47148f08705cc640f10a0290928344933239e5d1e47876b1"),
    ("equilibrium_thermo_model_header",
     "src/thermo/EquilibriumThermoModel.H", 10136,
     "d6a04013e33b6a34ae23735203be71505ab16b471f9bf01239efc1659141f2e2"),
    ("state_contract_header", "src/variables/StateContract.H", 2088,
     "4ad881962f93df46aab9024c0a7968dc62bfcbbf367bacfd9800f7747eaad186"),
    ("real_fluid_state_header", "src/thermo/RealFluidState.H", 3920,
     "c76d35eb5749957444470fe748e3c9055dcef457a39174fd66aaff35bdb3e499"),
    ("peng_robinson_header",
     "src/thermo/EquationOfState/PengRobinsonEOS.H", 2966,
     "4dcceec1a615c59c53b86f14c95a6118d94cef8d8e45db121d27c57c906c3e1f"),
    ("species_critical_db_header",
     "src/thermo/MixingRules/SpeciesCriticalDB.H", 973,
     "f82fbd6671873d0120baf081d4ec651901151f8050bb03b22e2ef2bf8ea6efde"),
    ("binary_interaction_db_header",
     "src/thermo/MixingRules/BinaryInteractionDB.H", 4058,
     "32e49334aa6a7e33ccd782879b6c510a2892311e45683c5834853a695ac9b908"),
    ("vdw_mixing_header",
     "src/thermo/MixingRules/VdWOneFluidMixing.H", 3068,
     "9919c59f501a4ec0097dbdac8c9e6fc1df587f44218d8aedbfe43f5fb769d664"),
    ("species_mapper", "src/combustion/SpeciesMapper.H", 2948,
     "5edb5efb615aac4d011205e87e6f0d6b58b6d3eea63dbdea913aeea0f29afa97"),
    ("two_step_source", "src/combustion/FiniteRate/TwoStepCO.C", 3538,
     "013d7af8f055a0ee834a095fa2946abd45214cf769ca5ccd43c3d83d16fd6465"),
    ("two_step_header", "src/combustion/FiniteRate/TwoStepCO.H", 3955,
     "9636dc4b2280af37d0231880a31672557556d5c9ee955bd4c46d373cad46d8c5"),
    ("dynamic_chemistry", "src/combustion/FiniteRate/DynamicChemistry.H", 15314,
     "04ad57fc8153982d5cdbe4744dcb873351caeece102e7e3bc3b6f5e39a3c33b6"),
    ("rfam_two_step_yaml", "tests/twostep_co.yaml", 1086,
     "f192aee746695a89d7e754b23d3d97a36da0e16fd510824db769e9590b92ae3c"),
    ("trajectory_contract",
     "tests/reference/r2_same_mechanism_trajectory_contract_v1.json", 7259,
     "301649763727d10947473fef504956439a15bb374f70a303bb8d7bda64d25af4"),
    ("probe_executable",
     "build-phase-foundation/rfam_r2_low_concentration_regularization_probe",
     42304,
     "5541eb0a9a82dfc0f1aeb0ec4d9ca69c2288c7e58631fd3fd3d9093e1e130a21"),
]

EXPECTED_APPROVAL_KEYS = {
    "regularization_threshold_selected_for_production",
    "threshold_convergence_validated",
    "detailed_mechanism_validation_complete",
    "experimental_holdout_validation_complete",
    "reacting_hot_start_approved",
    "production_chemistry_ready",
}

EXPECTED_COMPILE_COMMAND = [
    "/usr/bin/g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
    "-Isrc", "-Isrc/thermo", "-Isrc/framework", "-Isrc/combustion",
    "-Isrc/combustion/FiniteRate", "-Isrc/variables",
    "-I/mnt/f/opt/sundials/include", "-Isrc/thermo/EquationOfState",
    "-Isrc/thermo/MixingRules", "-Isrc/thermo/RootSelection",
    "-Isrc/numerics/flux", "-Isrc/numerics/preconditioning",
    "-Isrc/numerics/reconstruction", "-Isrc/mesh", "-Isrc/transport",
    "-Isrc/turbulence", "-Isrc/combustion/FPV",
    "tools/r2_low_concentration_regularization_probe.cpp", "-o",
    "build-phase-foundation/rfam_r2_low_concentration_regularization_probe",
]

EXPECTED_EXPONENTS = [0.1, 0.5, 1.0, 1.65]
EXPECTED_RATIOS = [
    0.0, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2,
    0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999999, 1.0, 1.000001, 10.0,
]

EXPECTED_RHS_GLOBAL_NAMES = [
    "O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2",
]
EXPECTED_RHS_CHEMISTRY_NAMES = ["C3H8", "O2", "CO", "CO2", "H2O", "N2"]
EXPECTED_RHS_CHEMISTRY_TO_GLOBAL = [5, 0, 7, 8, 6, 1]
EXPECTED_RHS_GLOBAL_MW = [
    0.0319988, 0.0280134, 0.00201588, 0.01604246, 0.03006904,
    0.04409562, 0.01801528, 0.0280101, 0.0440095,
]
EXPECTED_RHS_TEMPLATES = {
    "aligned": [0.20, 0.20, 0.00, 0.00, 0.00, 0.20, 0.10, 0.10, 0.20],
    "trace": [0.20, 0.20, 0.01, 0.02, 0.03, 0.20, 0.10, 0.10, 0.14],
}
EXPECTED_RHS_CASES = [
    ("aligned-rho1-ratio0.001", 1.0, 0.001, "aligned"),
    ("aligned-rho1-ratio0.5", 1.0, 0.5, "aligned"),
    ("trace-rho1-ratio0.5", 1.0, 0.5, "trace"),
    ("trace-rho0.1-ratio0.5", 0.1, 0.5, "trace"),
    ("trace-rho10-ratio0.5", 10.0, 0.5, "trace"),
    ("trace-rho1-ratio0.9", 1.0, 0.9, "trace"),
]


def duplicate_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def load_json_bytes(data: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(data.decode("utf-8"), object_pairs_hook=duplicate_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"json-decode-failed:{label}") from error
    if not isinstance(value, dict):
        raise ContractError(f"json-root-not-object:{label}")
    return value


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ContractError(f"artifact-open-failed:{path}") from error


def exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != expected:
        raise ContractError(f"object-shape-mismatch:{label}")
    return value


def finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"number-invalid:{label}")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"number-nonfinite:{label}")
    return result


def close(actual: float, expected: float, relative: float, label: str) -> None:
    scale = max(abs(expected), 1.0e-300)
    if abs(actual - expected) > relative * scale:
        raise ContractError(f"numeric-mismatch:{label}")


def state_expected(concentration: float, transition: float) -> float:
    if concentration <= 0.0:
        return 0.0
    if concentration >= transition:
        return concentration
    scaled = concentration / transition
    return transition * scaled * scaled * (2.0 - scaled)


def power_expected(concentration: float, exponent: float, transition: float) -> float:
    if exponent == 0.0:
        return 1.0
    if concentration <= 0.0:
        return 0.0
    if concentration >= transition or exponent >= 1.0:
        return concentration ** exponent
    scaled = concentration / transition
    return transition ** exponent * (
        (3.0 - exponent) * scaled * scaled
        + (exponent - 2.0) * scaled * scaled * scaled
    )


def validate_rhs_probe_contract(contract: dict[str, Any]) -> dict[str, Any]:
    section = exact_keys(contract["rhs_composition_probe"], {
        "global_species_names", "chemistry_species_names",
        "chemistry_to_global", "global_molecular_weights_kg_mol",
        "templates", "oxygen_remainder_fraction",
        "nitrogen_remainder_fraction", "cases", "scope",
    }, "rhs-composition-probe")
    if section["global_species_names"] != EXPECTED_RHS_GLOBAL_NAMES:
        raise ContractError("rhs-global-species-mismatch")
    if section["chemistry_species_names"] != EXPECTED_RHS_CHEMISTRY_NAMES:
        raise ContractError("rhs-chemistry-species-mismatch")
    if section["chemistry_to_global"] != EXPECTED_RHS_CHEMISTRY_TO_GLOBAL:
        raise ContractError("rhs-mapping-mismatch")
    if section["global_molecular_weights_kg_mol"] != EXPECTED_RHS_GLOBAL_MW:
        raise ContractError("rhs-global-molecular-weights-mismatch")
    if finite(section["oxygen_remainder_fraction"], "rhs-oxygen-fraction") != 0.23:
        raise ContractError("rhs-oxygen-fraction-mismatch")
    if finite(section["nitrogen_remainder_fraction"], "rhs-nitrogen-fraction") != 0.77:
        raise ContractError("rhs-nitrogen-fraction-mismatch")
    if section["scope"] != (
        "exact RhsCompositionPreparation::prepare path called by rhsFn, "
        "fixed diagnostic backgrounds only"
    ):
        raise ContractError("rhs-scope-mismatch")
    templates = exact_keys(section["templates"], {"aligned", "trace"}, "rhs-templates")
    if templates != EXPECTED_RHS_TEMPLATES:
        raise ContractError("rhs-template-mismatch")
    cases = section["cases"]
    if not isinstance(cases, list) or len(cases) != len(EXPECTED_RHS_CASES):
        raise ContractError("rhs-case-count-mismatch")
    for item, expected in zip(cases, EXPECTED_RHS_CASES):
        item = exact_keys(item, {
            "id", "rho_kg_m3", "fuel_concentration_ratio", "template",
        }, "rhs-case")
        if item["id"] != expected[0] or item["template"] != expected[3]:
            raise ContractError("rhs-case-identity-mismatch")
        if finite(item["rho_kg_m3"], "rhs-case-rho") != expected[1]:
            raise ContractError("rhs-case-rho-mismatch")
        if finite(item["fuel_concentration_ratio"], "rhs-case-ratio") != expected[2]:
            raise ContractError("rhs-case-ratio-mismatch")
    return section


def validate_contract(root: Path, contract_path: Path) -> dict[str, Any]:
    expected_contract = (root / EXPECTED_CONTRACT_PATH).resolve()
    if contract_path.resolve() != expected_contract:
        raise ContractError("contract-path-mismatch")
    if sha256(contract_path) != EXPECTED_CONTRACT_SHA256:
        raise ContractError("contract-sha256-mismatch")
    contract = load_json_bytes(contract_path.read_bytes(), "contract")
    exact_keys(contract, {
        "schema", "contract_version", "status", "scope", "approvals",
        "artifacts", "toolchain", "regularization", "probe_grid",
        "rhs_composition_probe",
        "trajectory_diagnostic_bound", "acceptance", "holds", "stop_conditions",
    }, "contract")
    if contract["schema"] != "rfam.r2-low-concentration-regularization-contract/v1":
        raise ContractError("schema-mismatch")
    if isinstance(contract["contract_version"], bool) or contract["contract_version"] != 1:
        raise ContractError("contract-version-mismatch")
    if contract["status"] != "blocked-numerical-characterization-only":
        raise ContractError("status-mismatch")
    if contract["scope"] != (
        "current-rfam-exact-rhs-composition-and-low-concentration-smoothing-"
        "fixed-backgrounds-no-fit-no-physical-validation"
    ):
        raise ContractError("scope-mismatch")

    approvals = exact_keys(contract["approvals"], EXPECTED_APPROVAL_KEYS, "approvals")
    if any(value is not False for value in approvals.values()):
        raise ContractError("approval-enabled")

    artifacts = contract["artifacts"]
    if not isinstance(artifacts, list) or len(artifacts) != len(EXPECTED_ARTIFACTS):
        raise ContractError("artifact-count-mismatch")
    for item, expected in zip(artifacts, EXPECTED_ARTIFACTS):
        item = exact_keys(item, {"id", "path", "size", "sha256"}, "artifact")
        artifact_id, relative_path, expected_size, expected_hash = expected
        if (item["id"], item["path"], item["size"], item["sha256"]) != expected:
            raise ContractError(f"artifact-contract-mismatch:{artifact_id}")
        path = (root / relative_path).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise ContractError(f"artifact-path-escapes-root:{artifact_id}") from error
        if not path.is_file() or path.stat().st_size != expected_size:
            raise ContractError(f"artifact-size-mismatch:{artifact_id}")
        if sha256(path) != expected_hash:
            raise ContractError(f"artifact-sha256-mismatch:{artifact_id}")

    toolchain = exact_keys(
        contract["toolchain"],
        {"compiler_path", "compiler_sha256", "compiler_version", "compile_command"},
        "toolchain",
    )
    if toolchain["compiler_path"] != "/usr/bin/g++":
        raise ContractError("compiler-path-mismatch")
    compiler = Path(toolchain["compiler_path"])
    if not compiler.is_file() or sha256(compiler) != (
        "0c0b987719385f8e242819dfef6a3e97b4b2540ce4035dbc11565df2966bd69e"
    ):
        raise ContractError("compiler-identity-mismatch")
    if toolchain["compiler_sha256"] != sha256(compiler):
        raise ContractError("compiler-contract-mismatch")
    if toolchain["compiler_version"] != "g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0":
        raise ContractError("compiler-version-contract-mismatch")
    version = subprocess.run(
        [str(compiler), "--version"], check=False, text=True, capture_output=True)
    if version.returncode != 0 or version.stdout.splitlines()[0] != toolchain["compiler_version"]:
        raise ContractError("compiler-version-runtime-mismatch")
    if toolchain["compile_command"] != EXPECTED_COMPILE_COMMAND:
        raise ContractError("compile-command-mismatch")

    regularization = exact_keys(contract["regularization"], {
        "transition_concentration_mol_m3", "state_smoothing",
        "irreversible_fractional_power_smoothing", "composition_order",
        "current_two_step_policy",
    }, "regularization")
    if finite(regularization["transition_concentration_mol_m3"], "transition") != 1.0e-6:
        raise ContractError("transition-mismatch")
    state_smoothing = exact_keys(regularization["state_smoothing"], {
        "source", "below_transition_formula", "at_or_above_transition",
        "nonpositive",
    }, "state-smoothing")
    if state_smoothing != {
        "source": "RhsCompositionPreparation.regularizedNonnegative",
        "below_transition_formula": "c0*x^2*(2-x), x=c/c0",
        "at_or_above_transition": "identity",
        "nonpositive": "zero",
    }:
        raise ContractError("state-smoothing-mismatch")
    fractional_smoothing = exact_keys(
        regularization["irreversible_fractional_power_smoothing"], {
            "source", "below_transition_formula", "at_or_above_transition",
            "exponent_at_or_above_one", "nonpositive",
        }, "fractional-smoothing")
    if fractional_smoothing != {
        "source": "MassActionThermo.regularizedConcentrationPower",
        "below_transition_formula": (
            "c0^a*((3-a)*x^2+(a-2)*x^3), x=c/c0"
        ),
        "at_or_above_transition": "exact-power",
        "exponent_at_or_above_one": "exact-power",
        "nonpositive": "zero",
    }:
        raise ContractError("fractional-smoothing-mismatch")
    if regularization["composition_order"] != [
        "state mass fraction is smoothed at the concentration-equivalent threshold",
        "all mass fractions are renormalized",
        "chemistry receives the smoothed state",
        "irreversible fractional concentration powers are smoothed again",
    ]:
        raise ContractError("composition-order-mismatch")
    policy = exact_keys(regularization["current_two_step_policy"], {
        "R1_C3H8_exponent", "R1_O2_exponent", "R1_reversible",
        "R2_CO_exponent", "R2_O2_exponent", "R2_reversible",
        "R2_fractional_power_policy",
    }, "current-two-step-policy")
    if policy != {
        "R1_C3H8_exponent": 0.1, "R1_O2_exponent": 1.65,
        "R1_reversible": False, "R2_CO_exponent": 1.0,
        "R2_O2_exponent": 0.5, "R2_reversible": True,
        "R2_fractional_power_policy": "exact-power-no-regularization",
    }:
        raise ContractError("two-step-policy-mismatch")

    grid = exact_keys(contract["probe_grid"], {
        "exponents", "concentration_ratios", "power_rows", "state_rows",
        "pre_normalization_combined_rows", "rhs_composition_rows",
    }, "probe-grid")
    if grid != {
        "exponents": EXPECTED_EXPONENTS,
        "concentration_ratios": EXPECTED_RATIOS,
        "power_rows": 64,
        "state_rows": 16,
        "pre_normalization_combined_rows": 32,
        "rhs_composition_rows": 6,
    }:
        raise ContractError("probe-grid-mismatch")

    acceptance = exact_keys(contract["acceptance"], {
        "finite_nonnegative", "monotone", "continuous_at_zero_and_transition",
        "first_derivative_continuous_at_zero_and_transition",
        "exact_at_or_above_transition", "exact_for_exponent_at_or_above_one",
        "smoothed_value_not_above_exact_fractional_power",
        "pre_normalization_combined_value_not_above_single_power_smoothing",
        "rhs_composition_matches_independent_reconstruction",
        "review_counterexample_reproduced",
        "transition_derivative_relative_tolerance",
        "zero_derivative_normalized_maximum", "numeric_relative_tolerance",
    }, "acceptance")
    for key in (
        "finite_nonnegative", "monotone", "continuous_at_zero_and_transition",
        "first_derivative_continuous_at_zero_and_transition",
        "exact_at_or_above_transition", "exact_for_exponent_at_or_above_one",
        "smoothed_value_not_above_exact_fractional_power",
        "pre_normalization_combined_value_not_above_single_power_smoothing",
        "rhs_composition_matches_independent_reconstruction",
        "review_counterexample_reproduced",
    ):
        if acceptance[key] is not True:
            raise ContractError(f"acceptance-weakened:{key}")
    if (
        finite(acceptance["transition_derivative_relative_tolerance"], "derivative-tol")
        != 5.0e-6
        or finite(acceptance["zero_derivative_normalized_maximum"], "zero-tol")
        != 3.0e-8
        or finite(acceptance["numeric_relative_tolerance"], "numeric-tol")
        != 5.0e-13
    ):
        raise ContractError("acceptance-tolerance-mismatch")

    validate_rhs_probe_contract(contract)

    return contract


def derive_oxygen_bound(root: Path, diagnostic: dict[str, Any]) -> float:
    trajectory = load_json_bytes(
        (root / "tests/reference/r2_same_mechanism_trajectory_contract_v1.json").read_bytes(),
        "trajectory-contract",
    )
    molecular_weights = [
        finite(value, f"molecular-weight-{index}")
        for index, value in enumerate(diagnostic["molecular_weights_kg_mol"])
    ]
    gas_constant = finite(
        diagnostic["universal_gas_constant_j_mol_k"], "gas-constant")
    maximum = 0.0
    for case in trajectory["cases"]:
        fractions = [finite(value, "trajectory-Y") for value in case["mass_fractions"]]
        mixture_mw = 1.0 / sum(
            fraction / molecular_weight
            for fraction, molecular_weight in zip(fractions, molecular_weights)
        )
        density = finite(case["pressure_pa"], "trajectory-pressure") * mixture_mw / (
            gas_constant * finite(case["temperature_k"], "trajectory-temperature")
        )
        oxygen = density * fractions[1] / molecular_weights[1]
        maximum = max(maximum, oxygen)
    return maximum


def validate_evidence(
    root: Path,
    contract: dict[str, Any],
    evidence: dict[str, Any],
) -> int:
    checks = 0
    exact_keys(evidence, {
        "schema", "status", "transition_mol_m3", "rows", "state_rows",
        "pre_normalization_combined_rows", "rhs_composition_rows",
        "derivatives", "state_derivative", "source_bound",
        "r2_fractional_power_policy", "production_ready",
    }, "evidence")
    if evidence["schema"] != "rfam.r2-low-concentration-regularization/v1":
        raise ContractError("evidence-schema-mismatch")
    if evidence["status"] != "numerical-characterization-only-production-blocked":
        raise ContractError("evidence-status-mismatch")
    if evidence["r2_fractional_power_policy"] != (
        "reversible-exact-power-no-regularization"
    ):
        raise ContractError("evidence-r2-policy-mismatch")
    if evidence["production_ready"] is not False:
        raise ContractError("evidence-production-enabled")
    transition = finite(evidence["transition_mol_m3"], "evidence-transition")
    if transition != 1.0e-6:
        raise ContractError("evidence-transition-mismatch")
    checks += 5

    tolerance = 5.0e-13
    rows = evidence["rows"]
    if not isinstance(rows, list) or len(rows) != 64:
        raise ContractError("power-row-count-mismatch")
    previous: dict[float, float] = {}
    for row, (exponent, ratio) in zip(
        rows,
        ((exponent, ratio) for exponent in EXPECTED_EXPONENTS for ratio in EXPECTED_RATIOS),
    ):
        row = exact_keys(row, {
            "exponent", "concentration_ratio", "concentration_mol_m3",
            "regularized", "exact", "regularized_to_exact",
        }, "power-row")
        actual_exponent = finite(row["exponent"], "row-exponent")
        actual_ratio = finite(row["concentration_ratio"], "row-ratio")
        if actual_exponent != exponent or actual_ratio != ratio:
            raise ContractError("power-row-order-mismatch")
        concentration = transition * ratio
        actual = finite(row["regularized"], "row-regularized")
        exact = concentration ** exponent if concentration > 0.0 else 0.0
        expected = power_expected(concentration, exponent, transition)
        close(finite(row["concentration_mol_m3"], "row-concentration"),
              concentration, tolerance, "row-concentration")
        close(actual, expected, tolerance, "row-regularized")
        close(finite(row["exact"], "row-exact"), exact, tolerance, "row-exact")
        if ratio == 0.0:
            if row["regularized_to_exact"] is not None:
                raise ContractError("zero-power-ratio-not-null")
        else:
            close(finite(row["regularized_to_exact"], "row-ratio-value"),
                  expected / exact, tolerance, "row-ratio-value")
        if actual < 0.0 or not math.isfinite(actual):
            raise ContractError("power-row-invalid")
        if exponent in previous and actual < previous[exponent]:
            raise ContractError("power-row-not-monotone")
        previous[exponent] = actual
        if ratio < 1.0 and exponent < 1.0 and actual > exact:
            raise ContractError("power-smoothing-above-exact")
        if (ratio >= 1.0 or exponent >= 1.0) and actual != exact:
            close(actual, exact, tolerance, "power-exact-branch")
    checks += len(rows)

    state_rows = evidence["state_rows"]
    if not isinstance(state_rows, list) or len(state_rows) != 16:
        raise ContractError("state-row-count-mismatch")
    previous_state = -1.0
    for row, ratio in zip(state_rows, EXPECTED_RATIOS):
        row = exact_keys(row, {
            "concentration_ratio", "input_concentration_mol_m3",
            "regularized_concentration_mol_m3", "regularized_to_input",
        }, "state-row")
        if finite(row["concentration_ratio"], "state-ratio") != ratio:
            raise ContractError("state-row-order-mismatch")
        concentration = transition * ratio
        expected = state_expected(concentration, transition)
        actual = finite(row["regularized_concentration_mol_m3"], "state-value")
        close(finite(row["input_concentration_mol_m3"], "state-input"),
              concentration, tolerance, "state-input")
        close(actual, expected, tolerance, "state-value")
        if actual < previous_state:
            raise ContractError("state-row-not-monotone")
        previous_state = actual
        if ratio == 0.0:
            if row["regularized_to_input"] is not None:
                raise ContractError("zero-state-ratio-not-null")
        else:
            close(finite(row["regularized_to_input"], "state-value-ratio"),
                  expected / concentration, tolerance, "state-value-ratio")
    checks += len(state_rows)

    combined_rows = evidence["pre_normalization_combined_rows"]
    if not isinstance(combined_rows, list) or len(combined_rows) != 32:
        raise ContractError("pre-normalization-combined-row-count-mismatch")
    for row, (exponent, ratio) in zip(
        combined_rows,
        ((exponent, ratio) for exponent in EXPECTED_EXPONENTS[:2]
         for ratio in EXPECTED_RATIOS),
    ):
        row = exact_keys(row, {
            "exponent", "concentration_ratio", "state_regularized_mol_m3",
            "combined_factor", "exact_factor", "combined_to_exact",
        }, "combined-row")
        if (finite(row["exponent"], "combined-exponent") != exponent
                or finite(row["concentration_ratio"], "combined-ratio") != ratio):
            raise ContractError("combined-row-order-mismatch")
        concentration = transition * ratio
        state_value = state_expected(concentration, transition)
        combined = power_expected(state_value, exponent, transition)
        single = power_expected(concentration, exponent, transition)
        exact = concentration ** exponent if concentration > 0.0 else 0.0
        close(finite(row["state_regularized_mol_m3"], "combined-state"),
              state_value, tolerance, "combined-state")
        close(finite(row["combined_factor"], "combined-factor"),
              combined, tolerance, "combined-factor")
        close(finite(row["exact_factor"], "combined-exact"),
              exact, tolerance, "combined-exact")
        if combined > single:
            raise ContractError("combined-above-single-smoothing")
        if ratio == 0.0:
            if row["combined_to_exact"] is not None:
                raise ContractError("zero-combined-ratio-not-null")
        else:
            close(finite(row["combined_to_exact"], "combined-value-ratio"),
                  combined / exact, tolerance, "combined-value-ratio")
    checks += len(combined_rows)

    rhs_rows = evidence["rhs_composition_rows"]
    if not isinstance(rhs_rows, list) or len(rhs_rows) != 6:
        raise ContractError("rhs-composition-row-count-mismatch")
    global_mw = [finite(value, "rhs-global-molecular-weight")
                 for value in EXPECTED_RHS_GLOBAL_MW]
    chemistry_map = EXPECTED_RHS_CHEMISTRY_TO_GLOBAL
    chemistry_mw = [global_mw[index] for index in chemistry_map]
    rhs_expected_by_id = {
        case[0]: case for case in EXPECTED_RHS_CASES
    }
    counterexample_relative = None
    for row in rhs_rows:
        row = exact_keys(row, {
            "id", "rho_kg_m3", "fuel_concentration_ratio",
            "global_template_sum", "input_chemistry_mass_fraction_sum",
            "pre_normalization_mass_fraction_sum",
            "normalized_global_mass_fractions",
            "normalized_chemistry_mass_fractions",
            "pre_normalized_fuel_concentration_mol_m3",
            "normalized_fuel_concentration_mol_m3",
            "shortcut_r1_fuel_factor", "rhs_r1_fuel_factor",
        }, "rhs-composition-row")
        row_id = row["id"]
        if row_id not in rhs_expected_by_id:
            raise ContractError("rhs-composition-row-identity-mismatch")
        _, rho, ratio, template_name = rhs_expected_by_id[row_id]
        if finite(row["rho_kg_m3"], "rhs-row-rho") != rho:
            raise ContractError("rhs-row-rho-mismatch")
        if finite(row["fuel_concentration_ratio"], "rhs-row-ratio") != ratio:
            raise ContractError("rhs-row-ratio-mismatch")
        template = EXPECTED_RHS_TEMPLATES[template_name]
        if finite(row["global_template_sum"], "rhs-template-sum") != sum(template):
            raise ContractError("rhs-template-sum-mismatch")
        transition = 1.0e-6
        fuel_transition = transition * chemistry_mw[0] / rho
        fuel_mass = fuel_transition * ratio
        remaining = 1.0 - fuel_mass
        chemistry_input = [fuel_mass, 0.23 * remaining, 0.0, 0.0, 0.0,
                           0.77 * remaining]
        close(finite(row["input_chemistry_mass_fraction_sum"], "rhs-input-sum"),
              sum(chemistry_input), 5.0e-13, "rhs-input-sum")
        global_before = [state_expected(value, transition * mw / rho)
                         for value, mw in zip(template, global_mw)]
        for chem_index, global_index in enumerate(chemistry_map):
            global_before[global_index] = state_expected(
                chemistry_input[chem_index], transition * chemistry_mw[chem_index] / rho)
        pre_sum = sum(global_before)
        close(finite(row["pre_normalization_mass_fraction_sum"], "rhs-pre-sum"),
              pre_sum, 5.0e-13, "rhs-pre-sum")
        normalized_global = [value / pre_sum for value in global_before]
        normalized_chemistry = [normalized_global[index] for index in chemistry_map]
        actual_global = row["normalized_global_mass_fractions"]
        actual_chemistry = row["normalized_chemistry_mass_fractions"]
        if not isinstance(actual_global, list) or len(actual_global) != 9:
            raise ContractError("rhs-global-vector-shape-mismatch")
        if not isinstance(actual_chemistry, list) or len(actual_chemistry) != 6:
            raise ContractError("rhs-chemistry-vector-shape-mismatch")
        for actual, expected in zip(actual_global, normalized_global):
            close(finite(actual, "rhs-global-value"), expected,
                  5.0e-13, "rhs-global-value")
        for actual, expected in zip(actual_chemistry, normalized_chemistry):
            close(finite(actual, "rhs-chemistry-value"), expected,
                  5.0e-13, "rhs-chemistry-value")
        pre_fuel_concentration = rho * global_before[chemistry_map[0]] / chemistry_mw[0]
        normalized_fuel_concentration = (
            rho * normalized_chemistry[0] / chemistry_mw[0])
        shortcut_factor = power_expected(pre_fuel_concentration, 0.1, transition)
        rhs_factor = power_expected(normalized_fuel_concentration, 0.1, transition)
        close(finite(row["pre_normalized_fuel_concentration_mol_m3"], "rhs-pre-fuel"),
              pre_fuel_concentration, 5.0e-13, "rhs-pre-fuel")
        close(finite(row["normalized_fuel_concentration_mol_m3"], "rhs-fuel"),
              normalized_fuel_concentration, 5.0e-13, "rhs-fuel")
        close(finite(row["shortcut_r1_fuel_factor"], "rhs-shortcut-factor"),
              shortcut_factor, 5.0e-13, "rhs-shortcut-factor")
        close(finite(row["rhs_r1_fuel_factor"], "rhs-factor"),
              rhs_factor, 5.0e-13, "rhs-factor")
        if row_id == "aligned-rho1-ratio0.5":
            counterexample_relative = abs(rhs_factor - shortcut_factor) / max(
                abs(shortcut_factor), 1.0e-300)
        checks += 1
    if counterexample_relative is None or counterexample_relative <= 5.0e-9:
        raise ContractError("rhs-counterexample-not-reproduced")
    checks += 9

    derivatives = evidence["derivatives"]
    if not isinstance(derivatives, list) or len(derivatives) != 2:
        raise ContractError("derivative-count-mismatch")
    derivative_tolerance = 5.0e-6
    for row, exponent in zip(derivatives, EXPECTED_EXPONENTS[:2]):
        row = exact_keys(row, {
            "exponent", "transition_left", "transition_right",
            "transition_expected", "zero_right_normalized",
        }, "derivative")
        if finite(row["exponent"], "derivative-exponent") != exponent:
            raise ContractError("derivative-order-mismatch")
        expected = exponent * transition ** (exponent - 1.0)
        close(finite(row["transition_expected"], "derivative-expected"),
              expected, tolerance, "derivative-expected")
        close(finite(row["transition_left"], "derivative-left"),
              expected, derivative_tolerance, "derivative-left")
        close(finite(row["transition_right"], "derivative-right"),
              expected, derivative_tolerance, "derivative-right")
        if finite(row["zero_right_normalized"], "derivative-zero") > 3.0e-8:
            raise ContractError("derivative-zero-too-large")
    checks += len(derivatives)

    state_derivative = exact_keys(evidence["state_derivative"], {
        "transition_left", "transition_right", "transition_expected",
        "zero_right_normalized",
    }, "state-derivative")
    if finite(state_derivative["transition_expected"], "state-derivative-expected") != 1.0:
        raise ContractError("state-derivative-expected-mismatch")
    close(finite(state_derivative["transition_left"], "state-derivative-left"),
          1.0, derivative_tolerance, "state-derivative-left")
    close(finite(state_derivative["transition_right"], "state-derivative-right"),
          1.0, derivative_tolerance, "state-derivative-right")
    if finite(state_derivative["zero_right_normalized"], "state-derivative-zero") > 3.0e-8:
        raise ContractError("state-derivative-zero-too-large")
    checks += 4

    diagnostic = exact_keys(contract["trajectory_diagnostic_bound"], {
        "temperature_bound_k", "oxygen_bound_rule", "oxygen_bound_mol_m3",
        "universal_gas_constant_j_mol_k", "molecular_weights_kg_mol", "r1",
        "expected", "scope",
    }, "diagnostic-bound")
    oxygen_bound = derive_oxygen_bound(root, diagnostic)
    close(oxygen_bound, finite(diagnostic["oxygen_bound_mol_m3"], "oxygen-bound"),
          tolerance, "derived-oxygen-bound")
    source = exact_keys(evidence["source_bound"], {
        "temperature_k", "oxygen_mol_m3", "k1", "maximum_fuel_factor",
        "progress_mol_m3_s", "fuel_mass_source_kg_m3_s", "heat_release_w_m3",
        "fuel_inventory_kg_m3", "depletion_time_s",
    }, "source-bound")
    expected_source = exact_keys(diagnostic["expected"], {
        "k1", "maximum_fuel_factor", "progress_mol_m3_s",
        "fuel_mass_source_kg_m3_s", "heat_release_w_m3",
        "fuel_inventory_kg_m3", "depletion_time_s",
    }, "expected-source")
    close(finite(source["temperature_k"], "source-temperature"),
          finite(diagnostic["temperature_bound_k"], "bound-temperature"),
          tolerance, "source-temperature")
    close(finite(source["oxygen_mol_m3"], "source-oxygen"), oxygen_bound,
          tolerance, "source-oxygen")
    for key, value in expected_source.items():
        close(finite(source[key], f"source-{key}"), finite(value, f"expected-{key}"),
              tolerance, f"source-{key}")
    if finite(source["depletion_time_s"], "source-depletion") >= 1.0e-15:
        raise ContractError("diagnostic-stiffness-hold-missing")
    checks += 10
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()
    try:
        root = args.root.resolve()
        if not root.is_dir():
            raise ContractError("root-not-directory")
        contract = validate_contract(root, args.contract.resolve())
        expected_executable = (root / EXPECTED_ARTIFACTS[-1][1]).resolve()
        if args.executable.resolve() != expected_executable:
            raise ContractError("executable-path-mismatch")
        diagnostic = contract["trajectory_diagnostic_bound"]
        oxygen_bound = derive_oxygen_bound(root, diagnostic)
        result = subprocess.run(
            [
                str(expected_executable),
                format(finite(diagnostic["temperature_bound_k"], "temperature"), ".17g"),
                format(oxygen_bound, ".17g"),
            ],
            cwd=root,
            check=False,
            text=True,
            capture_output=True,
            timeout=30,
        )
        if result.returncode != 0:
            raise ContractError("probe-return-code")
        if result.stderr:
            raise ContractError("probe-stderr-not-empty")
        lines = result.stdout.splitlines()
        if len(lines) != 2 or lines[1] != "R2_LOW_CONCENTRATION_PROBE_OK":
            raise ContractError("probe-output-shape-mismatch")
        evidence = load_json_bytes(lines[0].encode("utf-8"), "probe-evidence")
        checks = validate_evidence(root, contract, evidence)
    except ContractError as error:
        print(f"R2_LOW_CONCENTRATION_REGULARIZATION_FAIL {error}", file=sys.stderr)
        return 2
    except Exception as error:
        print(
            "R2_LOW_CONCENTRATION_REGULARIZATION_FAIL unexpected:"
            + type(error).__name__,
            file=sys.stderr,
        )
        return 2

    print(
        "R2_LOW_CONCENTRATION_REGULARIZATION_BLOCKED "
        f"checks={checks} power_rows=64 state_rows=16 "
        "pre_normalization_combined_rows=32 rhs_composition_rows=6 "
        "threshold_convergence_validated=0 physical_validation_approved=0 "
        "production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
