#!/usr/bin/env python3
"""Validate the candidate-only O2/C3H8 topology uncertainty tube."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Iterable


CONTRACT_SCHEMA = "rfam.o2-c3h8-topology-uncertainty-tube-contract/v1"
EVIDENCE_SCHEMA = "rfam.o2-c3h8-topology-uncertainty-tube-evidence/v1"
FAIL_MARKER = "RFAM_O2_C3H8_TOPOLOGY_UNCERTAINTY_TUBE_FAIL"

STATUS_THREE_PHASE = "THREE_PHASE_REQUIRED"
STATUS_OUTSIDE = "OUT_OF_ENVELOPE"
STATUS_CLEAR = "SCREEN_ONLY_NO_MODEL_RISK_DETECTED"

CONTRACT_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_topology_uncertainty_tube_contract_v1.json"
)
BASE_VALIDATOR_PATH = (
    "realFluidAllMach/tools/validate_o2_c3h8_topology_envelope.py"
)
BASE_CONTRACT_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_topology_envelope_contract_v1.json"
)
TRAJECTORY_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_predicted_trajectory_512_v1.csv"
)
RISK_PATH = (
    "realFluidAllMach/tests/reference/"
    "o2_c3h8_topology_model_risk_v1.json"
)

EXPECTED_AUTHORITIES = {
    "base_topology_envelope_validator": {
        "id": "base-topology-envelope-validator",
        "path": BASE_VALIDATOR_PATH,
        "sha256": "bb4284f0eae6b7da1d85594d22391a431f3bf7608e81cf85a421a49995eec284",
        "byte_count": 38931,
    },
    "base_topology_envelope_contract": {
        "id": "base-topology-envelope-contract",
        "path": BASE_CONTRACT_PATH,
        "sha256": "cf5c8855b3afe7fd1e5c4d67366e046d4b7b3515a381be4e78af4281e18763ca",
        "byte_count": 8890,
    },
    "trajectory_reference": {
        "id": "predicted-trajectory",
        "path": TRAJECTORY_PATH,
        "sha256": "a6d498887225a685b2f7b5d40248913ad2e7ba93158e12bbe64963a01783cdc0",
        "byte_count": 37328,
    },
    "model_risk_registry": {
        "id": "model-risk-registry",
        "path": RISK_PATH,
        "sha256": "3ef2f001f56af5cbdb8fa851992736ec517cb17c5ed1656f32b85d9946d072b2",
        "byte_count": 3361,
    },
}

EXPECTED_COMPONENTS = {
    "trial_kij_interval": {
        "basis": "per-anchor-trial-kij-offset-not-fit-posterior",
        "nominal_anchor_kij_range": [-0.05, 0.15],
        "kij_offset": [-0.005, 0.005],
        "kij_fitted": False,
        "physical_validation": False,
    },
    "inlet_condition_error": {
        "basis": "declared-screening-input-error-not-measurement-confidence",
        "pressure_pa_offset": [-50000.0, 50000.0],
        "temperature_k_offset": [-0.2, 0.2],
        "overall_z_o2_offset": [-0.002, 0.002],
    },
    "numerical_interpolation_error": {
        "basis": "absolute-bounds-added-by-minkowski-sum",
        "pressure_pa_offset": [-1000.0, 1000.0],
        "temperature_k_offset": [-0.01, 0.01],
        "overall_z_o2_offset": [-0.00001, 0.00001],
        "kij_offset": [-0.00001, 0.00001],
    },
    "topology_detector_margin": {
        "basis": "risk-box-expansion-from-bound-model-risk-registry",
        "registry_field": "approach_guard_offsets",
        "pressure_pa_offset": [-50000.0, 50000.0],
        "temperature_k_offset": [-0.05, 0.05],
        "overall_z_o2_offset": [0.0, 0.0],
        "kij_offset": [-0.001, 0.001],
    },
}

EXPECTED_INTERVAL_POLICY = {
    "anchor_bound": (
        "nominal-plus-trial-kij-plus-inlet-error-plus-"
        "numerical-interpolation-error"
    ),
    "risk_bound": "registered-risk-box-plus-topology-detector-margin",
    "segment_bound": "axis-aligned-hull-of-adjacent-effective-anchor-boxes",
    "intersection": "closed-interval-all-axes",
    "admission_rule": (
        "if-any-realization-contacts-any-risk-box-return-three-phase-required"
    ),
    "sampling_role": "diagnostic-only-never-admission-authority",
    "uncertainty_widening_must_be_monotone": True,
}

EXPECTED_SUPPORTED = {
    "pressure_pa": [4000000.0, 6000000.0],
    "temperature_k": [110.0, 300.0],
    "kij": [-0.05, 0.15],
    "overall_z_o2": [0.0, 1.0],
}

EXPECTED_RISK_POLICY = {
    "pair": ["O2", "C3H8"],
    "intersection_status": STATUS_THREE_PHASE,
    "outside_or_uncommittable_status": STATUS_OUTSIDE,
    "clear_status": STATUS_CLEAR,
    "status_precedence": [STATUS_THREE_PHASE, STATUS_OUTSIDE, STATUS_CLEAR],
    "physical_vlle_claim": False,
    "z_localization_policy": "conservative-all-compositions-until-validated",
}

EXPECTED_APPROVALS = {
    "trial_kij_only": True,
    "physical_validation": False,
    "kij_fitted": False,
    "three_phase_owner_implemented": False,
    "hem_enabled": False,
    "fluent_execution_approved": False,
    "hot_start_approved": False,
    "production_ready": False,
}

EXPECTED_LIMITATIONS = [
    "the kij interval is a trial screening offset and is not a fitted posterior or confidence interval",
    "inlet and numerical intervals are declared screen inputs and are not validated measurement or model-form uncertainty",
    "the detector margin expands model-risk boxes and does not prove physical VLLE",
    "the possible-three-phase registry is PR model-risk evidence and overall composition remains conservatively unlocalized",
    "axis-aligned Minkowski bounds discard correlation and may over-reject but are not a joint probability model",
    "a clear finite tube does not prove global topology uniqueness or physical validity",
    "no TP or PH runtime gate, three-phase owner, HEM capability, Fluent execution, hot-start, or production flag is enabled",
]

AXES = ("pressure_pa", "temperature_k", "kij", "overall_z_o2")


class TubeError(ValueError):
    """A malformed authority or uncertainty tube fails closed."""


@dataclass(frozen=True)
class Interval:
    low: float
    high: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.low) or not math.isfinite(self.high):
            raise TubeError("nonfinite-interval")
        if self.low > self.high:
            raise TubeError("reversed-interval")

    def shifted(self, offset: "Interval") -> "Interval":
        return Interval(self.low + offset.low, self.high + offset.high)

    def hull(self, other: "Interval") -> "Interval":
        return Interval(min(self.low, other.low), max(self.high, other.high))

    def overlaps(self, other: "Interval") -> bool:
        return self.low <= other.high and other.low <= self.high

    def contained_by(self, other: "Interval") -> bool:
        return other.low <= self.low and self.high <= other.high

    def contains(self, value: float) -> bool:
        return self.low <= value <= self.high

    def samples(self) -> tuple[float, ...]:
        if self.low == self.high:
            return (self.low,)
        return (self.low, 0.5 * (self.low + self.high), self.high)


@dataclass(frozen=True)
class Anchor:
    identifier: str
    group: str
    order: float
    pressure: float
    temperature: float
    kij: float
    z_o2: float
    committable: bool = True


@dataclass(frozen=True)
class RiskRecord:
    identifier: str
    pressure: float
    temperature: float
    kij: float
    z_interval: Interval


Box = dict[str, Interval]


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise TubeError(f"duplicate-json-key:{key}")
        result[key] = value
    return result


def _read_bytes(path: Path, label: str) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise TubeError(f"{label}-read-failed:{exc}") from exc


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    payload = _read_bytes(path, label)
    if b"\r" in payload:
        raise TubeError(f"{label}-not-canonical-lf")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TubeError(f"{label}-not-ascii") from exc
    try:
        value = json.loads(text, object_pairs_hook=_reject_duplicates)
    except TubeError:
        raise
    except (json.JSONDecodeError, ValueError) as exc:
        raise TubeError(f"{label}-invalid-json:{exc}") from exc
    if not isinstance(value, dict):
        raise TubeError(f"{label}-root-not-object")
    return value, payload


def _strict_equal(actual: Any, expected: Any) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            _strict_equal(actual[key], expected[key]) for key in expected
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _strict_equal(left, right) for left, right in zip(actual, expected)
        )
    return actual == expected


def _exact(actual: Any, expected: Any, label: str) -> None:
    if not _strict_equal(actual, expected):
        raise TubeError(f"value-mismatch:{label}")


def _object(value: Any, label: str, keys: Iterable[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TubeError(f"type-mismatch:{label}:object")
    expected = set(keys)
    actual = set(value)
    if actual != expected:
        raise TubeError(
            f"key-mismatch:{label}:missing={sorted(expected - actual)}:"
            f"extra={sorted(actual - expected)}"
        )
    return value


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TubeError(f"type-mismatch:{label}:number")
    result = float(value)
    if not math.isfinite(result):
        raise TubeError(f"nonfinite:{label}")
    return result


def _interval(value: Any, label: str, contains_zero: bool = False) -> Interval:
    if not isinstance(value, list) or len(value) != 2:
        raise TubeError(f"type-mismatch:{label}:two-number-array")
    result = Interval(
        _number(value[0], f"{label}[0]"),
        _number(value[1], f"{label}[1]"),
    )
    if contains_zero and not result.contains(0.0):
        raise TubeError(f"offset-interval-must-contain-zero:{label}")
    return result


def _workspace_path(root: Path, relative: str, label: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise TubeError(f"invalid-relative-path:{label}")
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise TubeError(f"path-escapes-root:{label}")
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise TubeError(f"path-escapes-root:{label}") from exc
    return resolved


def _validate_live_authority(root: Path, name: str) -> bytes:
    authority = EXPECTED_AUTHORITIES[name]
    path = _workspace_path(root, authority["path"], name)
    payload = _read_bytes(path, name)
    if len(payload) != authority["byte_count"]:
        raise TubeError(f"{name}-byte-count-mismatch")
    if _sha256(payload) != authority["sha256"]:
        raise TubeError(f"{name}-sha256-mismatch")
    return payload


def load_base_validator(root: Path) -> ModuleType:
    # Hash before import so changed Python bytes cannot execute as authority.
    _validate_live_authority(root, "base_topology_envelope_validator")
    path = _workspace_path(root, BASE_VALIDATOR_PATH, "base-validator")
    module_name = "_rfam_bound_topology_envelope_" + _sha256(
        str(path).encode("utf-8")
    )[:16]
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise TubeError("base-validator-import-spec-failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        sys.modules.pop(module_name, None)
        raise TubeError(f"base-validator-import-failed:{exc}") from exc
    return module


def _parse_components(document: dict[str, Any]) -> dict[str, Any]:
    components = document["tube_components"]
    trial = components["trial_kij_interval"]
    inlet = components["inlet_condition_error"]
    numerical = components["numerical_interpolation_error"]
    detector = components["topology_detector_margin"]
    zero = Interval(0.0, 0.0)
    return {
        "trial_nominal_range": _interval(
            trial["nominal_anchor_kij_range"], "trial.nominal-kij-range"
        ),
        "trial_kij_offset": _interval(
            trial["kij_offset"], "trial.kij-offset", True
        ),
        "inlet": {
            "pressure_pa": _interval(
                inlet["pressure_pa_offset"], "inlet.pressure-offset", True
            ),
            "temperature_k": _interval(
                inlet["temperature_k_offset"], "inlet.temperature-offset", True
            ),
            "kij": zero,
            "overall_z_o2": _interval(
                inlet["overall_z_o2_offset"], "inlet.z-offset", True
            ),
        },
        "numerical": {
            "pressure_pa": _interval(
                numerical["pressure_pa_offset"],
                "numerical.pressure-offset",
                True,
            ),
            "temperature_k": _interval(
                numerical["temperature_k_offset"],
                "numerical.temperature-offset",
                True,
            ),
            "kij": _interval(
                numerical["kij_offset"], "numerical.kij-offset", True
            ),
            "overall_z_o2": _interval(
                numerical["overall_z_o2_offset"],
                "numerical.z-offset",
                True,
            ),
        },
        "detector": {
            "pressure_pa": _interval(
                detector["pressure_pa_offset"],
                "detector.pressure-offset",
                True,
            ),
            "temperature_k": _interval(
                detector["temperature_k_offset"],
                "detector.temperature-offset",
                True,
            ),
            "kij": _interval(
                detector["kij_offset"], "detector.kij-offset", True
            ),
            "overall_z_o2": _interval(
                detector["overall_z_o2_offset"], "detector.z-offset", True
            ),
        },
    }


def supported_set(value: dict[str, Any]) -> Box:
    return {
        "pressure_pa": _interval(value["pressure_pa"], "supported.pressure"),
        "temperature_k": _interval(
            value["temperature_k"], "supported.temperature"
        ),
        "kij": _interval(value["kij"], "supported.kij"),
        "overall_z_o2": _interval(value["overall_z_o2"], "supported.z"),
    }


def validate_contract(
    root: Path, path: Path
) -> tuple[dict[str, Any], list[Anchor], list[RiskRecord]]:
    expected_path = _workspace_path(root, CONTRACT_PATH, "contract")
    if path.resolve() != expected_path:
        raise TubeError("contract-path-mismatch")
    document, _ = _load_json(path, "contract")
    root_object = _object(
        document,
        "contract",
        {
            "schema",
            "contract_version",
            "status",
            "scope",
            "authorities",
            "tube_components",
            "interval_policy",
            "supported_screen_envelope",
            "risk_policy",
            "approvals",
            "limitations",
        },
    )
    _exact(root_object["schema"], CONTRACT_SCHEMA, "contract.schema")
    _exact(root_object["contract_version"], 1, "contract.contract-version")
    _exact(
        root_object["status"],
        "offline-trial-model-risk-gate-blocked",
        "contract.status",
    )
    _exact(
        root_object["scope"],
        (
            "current-rfam-pr-o2-c3h8-predicted-trajectory-explicit-"
            "uncertainty-tube-no-physical-vlle"
        ),
        "contract.scope",
    )
    _exact(root_object["authorities"], EXPECTED_AUTHORITIES, "contract.authorities")
    _exact(root_object["tube_components"], EXPECTED_COMPONENTS, "contract.tube-components")
    _exact(root_object["interval_policy"], EXPECTED_INTERVAL_POLICY, "contract.interval-policy")
    _exact(root_object["supported_screen_envelope"], EXPECTED_SUPPORTED, "contract.supported")
    _exact(root_object["risk_policy"], EXPECTED_RISK_POLICY, "contract.risk-policy")
    _exact(root_object["approvals"], EXPECTED_APPROVALS, "contract.approvals")
    _exact(root_object["limitations"], EXPECTED_LIMITATIONS, "contract.limitations")

    for name in EXPECTED_AUTHORITIES:
        _validate_live_authority(root, name)

    base = load_base_validator(root)
    base_contract_path = _workspace_path(root, BASE_CONTRACT_PATH, "base-contract")
    try:
        base_document, base_risks, trajectory_payload, _ = base.validate_contract(
            root, base_contract_path
        )
        base_anchors = base.parse_trajectory(trajectory_payload)
    except Exception as exc:
        raise TubeError(f"base-envelope-validation-failed:{exc}") from exc

    # The new names split the accepted base arithmetic without changing values.
    base_uncertainty = base_document["uncertainty_intervals"]
    _exact(
        base_uncertainty,
        {
            "basis": "declared-screening-intervals-not-physical-uncertainty",
            "kij_offset": root_object["tube_components"]["trial_kij_interval"]["kij_offset"],
            "temperature_k_offset": root_object["tube_components"]["inlet_condition_error"]["temperature_k_offset"],
            "pressure_pa_offset": root_object["tube_components"]["inlet_condition_error"]["pressure_pa_offset"],
            "overall_z_o2_offset": root_object["tube_components"]["inlet_condition_error"]["overall_z_o2_offset"],
        },
        "base-contract.uncertainty-component-equivalence",
    )
    _exact(
        base_document["numerical_interpolation_error_intervals"],
        {
            "basis": "absolute-bounds-added-by-minkowski-sum",
            "kij_offset": root_object["tube_components"]["numerical_interpolation_error"]["kij_offset"],
            "temperature_k_offset": root_object["tube_components"]["numerical_interpolation_error"]["temperature_k_offset"],
            "pressure_pa_offset": root_object["tube_components"]["numerical_interpolation_error"]["pressure_pa_offset"],
            "overall_z_o2_offset": root_object["tube_components"]["numerical_interpolation_error"]["overall_z_o2_offset"],
        },
        "base-contract.numerical-component-equivalence",
    )
    _exact(
        base_document["supported_screen_envelope"],
        root_object["supported_screen_envelope"],
        "base-contract.supported-equivalence",
    )
    risk_document, _ = _load_json(
        _workspace_path(root, RISK_PATH, "risk-registry"), "risk-registry"
    )
    detector = root_object["tube_components"]["topology_detector_margin"]
    _exact(
        risk_document["approach_guard_offsets"],
        {
            "pressure_pa_offset": detector["pressure_pa_offset"],
            "temperature_k_offset": detector["temperature_k_offset"],
            "kij_offset": detector["kij_offset"],
            "overall_z_o2_offset": detector["overall_z_o2_offset"],
        },
        "risk-registry.detector-margin-equivalence",
    )

    anchors = [
        Anchor(
            identifier=item.identifier,
            group=item.group,
            order=float(item.order),
            pressure=float(item.pressure),
            temperature=float(item.temperature),
            kij=float(item.kij),
            z_o2=float(item.z_o2),
            committable=bool(item.committable),
        )
        for item in base_anchors
    ]
    risks = [
        RiskRecord(
            identifier=item.identifier,
            pressure=float(item.pressure),
            temperature=float(item.temperature),
            kij=float(item.kij),
            z_interval=Interval(
                float(item.z_interval.low), float(item.z_interval.high)
            ),
        )
        for item in base_risks
    ]
    _parse_components(root_object)
    return root_object, anchors, risks


def zero_offsets() -> Box:
    zero = Interval(0.0, 0.0)
    return {name: zero for name in AXES}


def anchor_box(
    anchor: Anchor,
    trial_kij_offset: Interval,
    inlet_error: Box,
    numerical_error: Box,
) -> Box:
    nominal = {
        "pressure_pa": Interval(anchor.pressure, anchor.pressure),
        "temperature_k": Interval(anchor.temperature, anchor.temperature),
        "kij": Interval(anchor.kij, anchor.kij),
        "overall_z_o2": Interval(anchor.z_o2, anchor.z_o2),
    }
    result = {
        name: nominal[name].shifted(inlet_error[name]).shifted(
            numerical_error[name]
        )
        for name in AXES
    }
    result["kij"] = result["kij"].shifted(trial_kij_offset)
    return result


def risk_box(record: RiskRecord, detector_margin: Box) -> Box:
    nominal = {
        "pressure_pa": Interval(record.pressure, record.pressure),
        "temperature_k": Interval(record.temperature, record.temperature),
        "kij": Interval(record.kij, record.kij),
        "overall_z_o2": record.z_interval,
    }
    return {
        name: nominal[name].shifted(detector_margin[name]) for name in AXES
    }


def hull_box(left: Box, right: Box) -> Box:
    return {name: left[name].hull(right[name]) for name in AXES}


def boxes_overlap(left: Box, right: Box) -> bool:
    return all(left[name].overlaps(right[name]) for name in AXES)


def box_contained(box: Box, supported: Box) -> bool:
    return all(box[name].contained_by(supported[name]) for name in AXES)


def point_in_box(point: dict[str, float], box: Box) -> bool:
    return all(box[name].contains(point[name]) for name in AXES)


def sample_box(box: Box) -> Iterable[dict[str, float]]:
    axes = [box[name].samples() for name in AXES]
    for values in itertools.product(*axes):
        yield dict(zip(AXES, values))


def _validate_offset_set(value: Box, label: str) -> None:
    if set(value) != set(AXES):
        raise TubeError(f"axis-mismatch:{label}")
    for name in AXES:
        if not value[name].contains(0.0):
            raise TubeError(f"offset-does-not-contain-zero:{label}.{name}")


def evaluate_tube(
    anchors: list[Anchor],
    risks: list[RiskRecord],
    trial_nominal_range: Interval,
    trial_kij_offset: Interval,
    inlet_error: Box,
    numerical_error: Box,
    detector_margin: Box,
    supported: Box,
) -> dict[str, Any]:
    if not anchors:
        raise TubeError("empty-trajectory")
    if not risks:
        raise TubeError("empty-risk-registry")
    if not trial_kij_offset.contains(0.0):
        raise TubeError("trial-kij-offset-does-not-contain-zero")
    for label, value in (
        ("inlet", inlet_error),
        ("numerical", numerical_error),
        ("detector", detector_margin),
    ):
        _validate_offset_set(value, label)
    if set(supported) != set(AXES):
        raise TubeError("axis-mismatch:supported")
    if len({item.identifier for item in anchors}) != len(anchors):
        raise TubeError("duplicate-anchor-id")
    if len({item.identifier for item in risks}) != len(risks):
        raise TubeError("duplicate-risk-id")

    anchor_boxes: list[tuple[Anchor, Box]] = []
    for anchor in anchors:
        numeric = (
            anchor.order,
            anchor.pressure,
            anchor.temperature,
            anchor.kij,
            anchor.z_o2,
        )
        if not all(math.isfinite(value) for value in numeric):
            raise TubeError(f"nonfinite-anchor:{anchor.identifier}")
        if not (0.0 <= anchor.z_o2 <= 1.0):
            raise TubeError(f"anchor-composition-bound:{anchor.identifier}")
        if not trial_nominal_range.contains(anchor.kij):
            raise TubeError(f"anchor-kij-outside-trial-range:{anchor.identifier}")
        anchor_boxes.append(
            (
                anchor,
                anchor_box(
                    anchor, trial_kij_offset, inlet_error, numerical_error
                ),
            )
        )

    risk_boxes: list[tuple[str, Box]] = []
    for record in risks:
        numeric = (record.pressure, record.temperature, record.kij)
        if not all(math.isfinite(value) for value in numeric):
            raise TubeError(f"nonfinite-risk:{record.identifier}")
        risk_boxes.append((record.identifier, risk_box(record, detector_margin)))

    groups: dict[str, list[tuple[Anchor, Box]]] = {}
    for anchor, box in anchor_boxes:
        groups.setdefault(anchor.group, []).append((anchor, box))
    segment_boxes: list[tuple[str, Box]] = []
    for group_name, values in groups.items():
        ordered = sorted(values, key=lambda item: item[0].order)
        orders = [item[0].order for item in ordered]
        if len(orders) != len(set(orders)):
            raise TubeError(f"duplicate-group-order:{group_name}")
        for index in range(len(ordered) - 1):
            segment_boxes.append(
                (
                    f"{group_name}|segment={index}",
                    hull_box(ordered[index][1], ordered[index + 1][1]),
                )
            )

    all_boxes: list[tuple[str, str, Box]] = [
        ("anchor", anchor.identifier, box) for anchor, box in anchor_boxes
    ]
    all_boxes.extend(
        ("segment", identifier, box) for identifier, box in segment_boxes
    )
    risk_hits: list[str] = []
    outside: list[str] = []
    sampled_risk_hits = 0
    sample_count = 0
    for kind, identifier, box in all_boxes:
        if not box_contained(box, supported):
            outside.append(f"{kind}:{identifier}")
        for risk_id, candidate in risk_boxes:
            if boxes_overlap(box, candidate):
                risk_hits.append(f"{kind}:{identifier}->{risk_id}")
        for point in sample_box(box):
            sample_count += 1
            if any(point_in_box(point, candidate) for _, candidate in risk_boxes):
                sampled_risk_hits += 1

    uncommittable = [item.identifier for item in anchors if not item.committable]
    if risk_hits:
        status = STATUS_THREE_PHASE
    elif outside or uncommittable:
        status = STATUS_OUTSIDE
    else:
        status = STATUS_CLEAR
    return {
        "schema": EVIDENCE_SCHEMA,
        "status": status,
        "anchors": len(anchor_boxes),
        "segments": len(segment_boxes),
        "model_risk_records": len(risk_boxes),
        "risk_box_hits": len(risk_hits),
        "sample_points": sample_count,
        "sampled_risk_hits": sampled_risk_hits,
        "out_of_envelope_boxes": len(outside),
        "uncommittable_anchors": len(uncommittable),
        "first_risk_hits": risk_hits[:8],
        "first_out_of_envelope": outside[:8],
        "interval_admission": True,
        "sampling_is_admission_authority": False,
        "trial_kij_only": True,
        "physical_validation": False,
        "kij_fitted": False,
        "three_phase_owner_implemented": False,
        "hem_enabled": False,
        "fluent_execution_approved": False,
        "hot_start_approved": False,
        "production_ready": False,
    }


def evaluate_contract(
    document: dict[str, Any],
    anchors: list[Anchor],
    risks: list[RiskRecord],
) -> dict[str, Any]:
    components = _parse_components(document)
    return evaluate_tube(
        anchors,
        risks,
        components["trial_nominal_range"],
        components["trial_kij_offset"],
        components["inlet"],
        components["numerical"],
        components["detector"],
        supported_set(document["supported_screen_envelope"]),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--contract", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    default_root = Path(__file__).resolve().parents[2]
    root = (args.root if args.root is not None else default_root).resolve()
    contract = (
        args.contract.resolve()
        if args.contract is not None
        else _workspace_path(root, CONTRACT_PATH, "contract")
    )
    try:
        document, anchors, risks = validate_contract(root, contract)
        evidence = evaluate_contract(document, anchors, risks)
    except Exception as exc:
        print(f"{FAIL_MARKER} reason={exc}")
        return 2

    print(json.dumps(evidence, sort_keys=True, separators=(",", ":")))
    print(
        "RFAM_O2_C3H8_TOPOLOGY_UNCERTAINTY_TUBE_"
        f"{evidence['status']} anchors={evidence['anchors']} "
        f"segments={evidence['segments']} "
        f"risk_hits={evidence['risk_box_hits']} "
        f"outside={evidence['out_of_envelope_boxes']} "
        "trial_kij_only=1 physical_validation=0 production_ready=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
