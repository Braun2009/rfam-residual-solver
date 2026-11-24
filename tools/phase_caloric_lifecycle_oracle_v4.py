#!/usr/bin/env python3
"""Pure identity/lifecycle oracle for the blocked phase-caloric v4 design.

This module owns no thermodynamic state and performs no phase recovery.  It is
only a deterministic model of the proposed token lifetime and publication
ordering so the design can be tested before a runtime owner is authorized.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import re
from typing import Any, Dict, Mapping, Optional, Tuple


UINT64_MAX = (1 << 64) - 1
SHA256_RE = re.compile(r"[0-9a-f]{64}")
OWNER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:/-]{0,127}")
PHASE_IDENTITIES = ("VAPOR", "SUPERCRITICAL")
ROOT_IDENTITIES = ("PR_VAPOR_ROOT", "PR_UNIQUE_SUPERCRITICAL_ROOT")
ENERGY_CONVENTION = "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT"


class LifecycleError(ValueError):
    """Fail-closed lifecycle-oracle rejection."""


def _uint64(value: Any, label: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise LifecycleError(f"invalid-uint64:{label}")
    if value < (1 if positive else 0) or value > UINT64_MAX:
        raise LifecycleError(f"uint64-out-of-range:{label}")
    return value


def _increment(value: int, label: str) -> int:
    _uint64(value, label)
    if value == UINT64_MAX:
        raise LifecycleError(f"uint64-overflow:{label}")
    return value + 1


def _sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise LifecycleError(f"invalid-sha256:{label}")
    if value == "0" * 64:
        raise LifecycleError(f"zero-sha256:{label}")
    return value


def _owner(value: Any) -> str:
    if not isinstance(value, str) or OWNER_RE.fullmatch(value) is None:
        raise LifecycleError("invalid-parent-owner-id")
    return value


def _canonical_value(value: Any, label: str) -> Any:
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        _uint64(value, label)
        return value
    if isinstance(value, list):
        return [
            _canonical_value(item, f"{label}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, tuple):
        return [
            _canonical_value(item, f"{label}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, dict):
        result: Dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise LifecycleError(f"invalid-canonical-key:{label}")
            result[key] = _canonical_value(item, f"{label}.{key}")
        return result
    # Floats are intentionally forbidden. Physical binary64 values enter a
    # fingerprint as exact f64be hexadecimal strings defined by the contract.
    raise LifecycleError(f"unsupported-canonical-type:{label}")


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise LifecycleError(f"duplicate-canonical-json-key:{key}")
        result[key] = value
    return result


def _reject_nonfinite_json(token: str) -> None:
    raise LifecycleError(f"nonfinite-canonical-json-number:{token}")


def canonical_json_bytes(value: Mapping[str, Any]) -> bytes:
    if not isinstance(value, dict):
        raise LifecycleError("canonical-payload-not-object")
    canonical = _canonical_value(dict(value), "payload")
    return json.dumps(
        canonical,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")


def parse_canonical_json_bytes(payload: bytes) -> Dict[str, Any]:
    if not isinstance(payload, bytes):
        raise LifecycleError("canonical-json-payload-not-bytes")
    try:
        value = json.loads(
            payload.decode("ascii"),
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_nonfinite_json,
        )
    except LifecycleError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LifecycleError("canonical-json-invalid") from exc
    if not isinstance(value, dict):
        raise LifecycleError("canonical-payload-not-object")
    if canonical_json_bytes(value) != payload:
        raise LifecycleError("canonical-json-bytes-not-normalized")
    return value


def canonical_fingerprint(domain: str, value: Mapping[str, Any]) -> str:
    if not isinstance(domain, str) or not domain.startswith("rfam."):
        raise LifecycleError("invalid-fingerprint-domain")
    try:
        prefix = domain.encode("ascii")
    except UnicodeEncodeError as exc:
        raise LifecycleError("non-ascii-fingerprint-domain") from exc
    digest = hashlib.sha256()
    digest.update(prefix)
    digest.update(b"\0")
    digest.update(canonical_json_bytes(value))
    return digest.hexdigest()


@dataclass(frozen=True)
class ParentBinding:
    parent_owner_id: str
    committed_parent_q_fingerprint: str
    parent_phase_payload_fingerprint: str
    parent_phase_identity: str
    root_identity: str
    model_id: str
    phase_caloric_schema_id: str
    phase_caloric_options_fingerprint: str
    ordered_species_fingerprint: str
    ordered_species_count: int
    energy_convention: str
    process_runtime_epoch: int
    parent_generation: int

    def __post_init__(self) -> None:
        _owner(self.parent_owner_id)
        _sha256(
            self.committed_parent_q_fingerprint,
            "committed-parent-q-fingerprint",
        )
        _sha256(
            self.parent_phase_payload_fingerprint,
            "parent-phase-payload-fingerprint",
        )
        _sha256(
            self.phase_caloric_options_fingerprint,
            "phase-caloric-options-fingerprint",
        )
        _sha256(
            self.ordered_species_fingerprint,
            "ordered-species-fingerprint",
        )
        _uint64(
            self.ordered_species_count,
            "ordered-species-count",
            positive=True,
        )
        if self.parent_phase_identity not in PHASE_IDENTITIES:
            raise LifecycleError("invalid-parent-phase-identity")
        if self.root_identity not in ROOT_IDENTITIES:
            raise LifecycleError("invalid-root-identity")
        if (
            self.parent_phase_identity == "VAPOR"
            and self.root_identity != "PR_VAPOR_ROOT"
        ) or (
            self.parent_phase_identity == "SUPERCRITICAL"
            and self.root_identity != "PR_UNIQUE_SUPERCRITICAL_ROOT"
        ):
            raise LifecycleError("phase-root-identity-mismatch")
        if not isinstance(self.model_id, str) or not self.model_id.strip():
            raise LifecycleError("invalid-model-id")
        if self.phase_caloric_schema_id != "rfam.phase-caloric/v4":
            raise LifecycleError("invalid-phase-caloric-schema-id")
        if self.energy_convention != ENERGY_CONVENTION:
            raise LifecycleError("invalid-energy-convention")
        _uint64(self.process_runtime_epoch, "process-runtime-epoch", positive=True)
        _uint64(self.parent_generation, "parent-generation", positive=True)

    def identity_payload(self) -> Dict[str, Any]:
        return {
            "committed_parent_q_fingerprint": self.committed_parent_q_fingerprint,
            "energy_convention": self.energy_convention,
            "model_id": self.model_id,
            "ordered_species_count": self.ordered_species_count,
            "ordered_species_fingerprint": self.ordered_species_fingerprint,
            "parent_generation": self.parent_generation,
            "parent_owner_id": self.parent_owner_id,
            "parent_phase_identity": self.parent_phase_identity,
            "parent_phase_payload_fingerprint": self.parent_phase_payload_fingerprint,
            "phase_caloric_options_fingerprint": self.phase_caloric_options_fingerprint,
            "phase_caloric_schema_id": self.phase_caloric_schema_id,
            "process_runtime_epoch": self.process_runtime_epoch,
            "root_identity": self.root_identity,
        }

    def identity_digest(self) -> str:
        return canonical_fingerprint(
            "rfam.phase-caloric.parent-binding/v4",
            self.identity_payload(),
        )


@dataclass(frozen=True)
class LiveToken:
    parent_owner_id: str
    transaction_id: str
    process_runtime_epoch: int
    parent_generation: int
    transaction_sequence: int
    parent_binding_digest: str
    root_identity: str
    token_digest: str


@dataclass(frozen=True)
class RestartHint:
    parent_owner_id: str
    source_process_runtime_epoch: int
    source_parent_generation: int
    committed_parent_q_fingerprint: str
    parent_phase_payload_fingerprint: str
    parent_phase_identity: str
    root_identity: str
    source_token_digest: str
    hint_digest: str


@dataclass
class _Transaction:
    transaction_id: str
    ordinal: int
    binding_digest: str
    current_sequence: int = 0
    current_token_digest: Optional[str] = None
    active: bool = True


class PhaseCaloricLifecycleOracle:
    """Deterministic, per-parent lifecycle model with no thermo behavior."""

    def __init__(self, process_runtime_epoch: int) -> None:
        self.process_runtime_epoch = _uint64(
            process_runtime_epoch, "process-runtime-epoch", positive=True
        )
        self._parents: Dict[str, ParentBinding] = {}
        self._transactions: Dict[str, _Transaction] = {}
        self._transaction_ordinals: Dict[str, int] = {}

    def register_parent(self, binding: ParentBinding) -> None:
        if binding.process_runtime_epoch != self.process_runtime_epoch:
            raise LifecycleError("parent-runtime-epoch-mismatch")
        if binding.parent_owner_id in self._parents:
            raise LifecycleError("parent-already-registered")
        self._parents[binding.parent_owner_id] = binding
        self._transaction_ordinals[binding.parent_owner_id] = 0

    def parent(self, parent_owner_id: str) -> ParentBinding:
        owner = _owner(parent_owner_id)
        try:
            return self._parents[owner]
        except KeyError as exc:
            raise LifecycleError("unknown-parent-owner") from exc

    def begin_transaction(self, parent_owner_id: str) -> str:
        binding = self.parent(parent_owner_id)
        existing = self._transactions.get(binding.parent_owner_id)
        if existing is not None and existing.active:
            raise LifecycleError("parent-transaction-already-active")
        ordinal = _increment(
            self._transaction_ordinals[binding.parent_owner_id],
            "transaction-ordinal",
        )
        self._transaction_ordinals[binding.parent_owner_id] = ordinal
        transaction_id = canonical_fingerprint(
            "rfam.phase-caloric.transaction/v4",
            {
                "parent_binding_digest": binding.identity_digest(),
                "parent_owner_id": binding.parent_owner_id,
                "transaction_ordinal": ordinal,
            },
        )
        self._transactions[binding.parent_owner_id] = _Transaction(
            transaction_id=transaction_id,
            ordinal=ordinal,
            binding_digest=binding.identity_digest(),
        )
        return transaction_id

    def _transaction(self, parent_owner_id: str) -> Tuple[ParentBinding, _Transaction]:
        binding = self.parent(parent_owner_id)
        transaction = self._transactions.get(binding.parent_owner_id)
        if transaction is None or not transaction.active:
            raise LifecycleError("no-active-parent-transaction")
        if transaction.binding_digest != binding.identity_digest():
            raise LifecycleError("transaction-parent-binding-drift")
        return binding, transaction

    @staticmethod
    def _token_digest(
        binding: ParentBinding,
        transaction: _Transaction,
        sequence: int,
    ) -> str:
        return canonical_fingerprint(
            "rfam.phase-caloric.live-root-token/v4",
            {
                "parent_binding_digest": binding.identity_digest(),
                "parent_generation": binding.parent_generation,
                "parent_owner_id": binding.parent_owner_id,
                "process_runtime_epoch": binding.process_runtime_epoch,
                "root_identity": binding.root_identity,
                "transaction_id": transaction.transaction_id,
                "transaction_sequence": sequence,
            },
        )

    def _issue(self, binding: ParentBinding, transaction: _Transaction) -> LiveToken:
        sequence = _increment(
            transaction.current_sequence,
            f"transaction-sequence:{binding.parent_owner_id}",
        )
        digest = self._token_digest(binding, transaction, sequence)
        transaction.current_sequence = sequence
        transaction.current_token_digest = digest
        return LiveToken(
            parent_owner_id=binding.parent_owner_id,
            transaction_id=transaction.transaction_id,
            process_runtime_epoch=binding.process_runtime_epoch,
            parent_generation=binding.parent_generation,
            transaction_sequence=sequence,
            parent_binding_digest=binding.identity_digest(),
            root_identity=binding.root_identity,
            token_digest=digest,
        )

    def issue_initial_token(self, parent_owner_id: str) -> LiveToken:
        binding, transaction = self._transaction(parent_owner_id)
        if transaction.current_sequence != 0:
            raise LifecycleError("initial-token-already-issued")
        return self._issue(binding, transaction)

    def continue_from_token(self, token: LiveToken) -> LiveToken:
        if not self.accepts_live_token(token):
            raise LifecycleError("stale-or-invalid-live-token")
        binding, transaction = self._transaction(token.parent_owner_id)
        return self._issue(binding, transaction)

    def accepts_live_token(self, token: Any) -> bool:
        if not isinstance(token, LiveToken):
            return False
        binding = self._parents.get(token.parent_owner_id)
        transaction = self._transactions.get(token.parent_owner_id)
        if binding is None or transaction is None or not transaction.active:
            return False
        if (
            token.process_runtime_epoch != self.process_runtime_epoch
            or token.process_runtime_epoch != binding.process_runtime_epoch
            or token.parent_generation != binding.parent_generation
            or token.parent_binding_digest != binding.identity_digest()
            or token.root_identity != binding.root_identity
            or token.transaction_id != transaction.transaction_id
            or token.transaction_sequence != transaction.current_sequence
            or token.token_digest != transaction.current_token_digest
        ):
            return False
        try:
            expected = self._token_digest(
                binding, transaction, token.transaction_sequence
            )
        except LifecycleError:
            return False
        return token.token_digest == expected

    def _invalidate_and_advance(
        self,
        parent_owner_id: str,
        committed_parent_q_fingerprint: str,
        parent_phase_payload_fingerprint: str,
        parent_phase_identity: str,
        root_identity: str,
    ) -> ParentBinding:
        binding, transaction = self._transaction(parent_owner_id)
        # Prepare and validate the complete replacement before mutating either
        # the transaction or its published parent binding.
        generation = _increment(binding.parent_generation, "parent-generation")
        replacement = replace(
            binding,
            committed_parent_q_fingerprint=_sha256(
                committed_parent_q_fingerprint,
                "committed-parent-q-fingerprint",
            ),
            parent_phase_payload_fingerprint=_sha256(
                parent_phase_payload_fingerprint,
                "parent-phase-payload-fingerprint",
            ),
            parent_phase_identity=parent_phase_identity,
            root_identity=root_identity,
            parent_generation=generation,
        )
        transaction.active = False
        transaction.current_token_digest = None
        self._parents[binding.parent_owner_id] = replacement
        return replacement

    def rollback(
        self,
        parent_owner_id: str,
        *,
        restored_committed_parent_q_fingerprint: str,
        restored_parent_phase_payload_fingerprint: str,
    ) -> ParentBinding:
        binding = self.parent(parent_owner_id)
        return self._invalidate_and_advance(
            binding.parent_owner_id,
            restored_committed_parent_q_fingerprint,
            restored_parent_phase_payload_fingerprint,
            binding.parent_phase_identity,
            binding.root_identity,
        )

    def final_reflash_publish(
        self,
        parent_owner_id: str,
        *,
        committed_parent_q_fingerprint: str,
        parent_phase_payload_fingerprint: str,
        parent_phase_identity: str,
        root_identity: str,
    ) -> ParentBinding:
        return self._invalidate_and_advance(
            parent_owner_id,
            committed_parent_q_fingerprint,
            parent_phase_payload_fingerprint,
            parent_phase_identity,
            root_identity,
        )

    def make_restart_hint(self, token: LiveToken) -> RestartHint:
        if not self.accepts_live_token(token):
            raise LifecycleError("cannot-persist-stale-token")
        binding = self.parent(token.parent_owner_id)
        payload = {
            "committed_parent_q_fingerprint": binding.committed_parent_q_fingerprint,
            "parent_generation": binding.parent_generation,
            "parent_owner_id": binding.parent_owner_id,
            "parent_phase_identity": binding.parent_phase_identity,
            "parent_phase_payload_fingerprint": binding.parent_phase_payload_fingerprint,
            "root_identity": binding.root_identity,
            "source_process_runtime_epoch": binding.process_runtime_epoch,
            "source_token_digest": token.token_digest,
        }
        return RestartHint(
            parent_owner_id=binding.parent_owner_id,
            source_process_runtime_epoch=binding.process_runtime_epoch,
            source_parent_generation=binding.parent_generation,
            committed_parent_q_fingerprint=binding.committed_parent_q_fingerprint,
            parent_phase_payload_fingerprint=binding.parent_phase_payload_fingerprint,
            parent_phase_identity=binding.parent_phase_identity,
            root_identity=binding.root_identity,
            source_token_digest=token.token_digest,
            hint_digest=canonical_fingerprint(
                "rfam.phase-caloric.restart-hint/v4", payload
            ),
        )

    def restart(self) -> "PhaseCaloricLifecycleOracle":
        restarted = PhaseCaloricLifecycleOracle(
            _increment(self.process_runtime_epoch, "process-runtime-epoch")
        )
        for binding in self._parents.values():
            restarted.register_parent(
                replace(
                    binding,
                    process_runtime_epoch=restarted.process_runtime_epoch,
                    parent_generation=_increment(
                        binding.parent_generation, "parent-generation"
                    ),
                )
            )
        return restarted

    def revalidate_restart_hint(
        self,
        hint: RestartHint,
    ) -> LiveToken:
        """Reuse numerical guidance only after exact current-parent replay.

        A hint is deliberately not consumed.  It may guide later generations
        whose complete parent fingerprints and phase/root identities replay
        exactly, but every reuse opens a fresh transaction and mints a new
        sequence-one LiveToken in the current runtime epoch.
        """
        if not isinstance(hint, RestartHint):
            raise LifecycleError("invalid-restart-hint-type")
        binding = self.parent(hint.parent_owner_id)
        if (
            hint.source_process_runtime_epoch >= self.process_runtime_epoch
            or hint.source_parent_generation >= binding.parent_generation
            or hint.committed_parent_q_fingerprint
                != binding.committed_parent_q_fingerprint
            or hint.parent_phase_payload_fingerprint
                != binding.parent_phase_payload_fingerprint
            or hint.parent_phase_identity != binding.parent_phase_identity
            or hint.root_identity != binding.root_identity
        ):
            raise LifecycleError("restart-hint-parent-revalidation-failed")
        expected_hint = canonical_fingerprint(
            "rfam.phase-caloric.restart-hint/v4",
            {
                "committed_parent_q_fingerprint": hint.committed_parent_q_fingerprint,
                "parent_generation": hint.source_parent_generation,
                "parent_owner_id": hint.parent_owner_id,
                "parent_phase_identity": hint.parent_phase_identity,
                "parent_phase_payload_fingerprint": hint.parent_phase_payload_fingerprint,
                "root_identity": hint.root_identity,
                "source_process_runtime_epoch": hint.source_process_runtime_epoch,
                "source_token_digest": hint.source_token_digest,
            },
        )
        if hint.hint_digest != expected_hint:
            raise LifecycleError("restart-hint-digest-mismatch")
        self.begin_transaction(binding.parent_owner_id)
        return self.issue_initial_token(binding.parent_owner_id)
