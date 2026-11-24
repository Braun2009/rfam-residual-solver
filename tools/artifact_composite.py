#!/usr/bin/env python3
"""Generate and verify deterministic artifact composite manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


SCHEMA = "rfam.artifact-composite/v1"
RESULT_SCHEMA = "rfam.artifact-composite-result/v1"
CONTRACT_VERSION = 1
ALGORITHM = "sha256(path-nul-size-u64be-bytes)"
PATH_ENCODING = "ascii"
MAX_U64 = (1 << 64) - 1
CHUNK_SIZE = 1024 * 1024

ROOT_KEYS = {
    "schema",
    "contract_version",
    "algorithm",
    "path_encoding",
    "artifact_count",
    "artifacts",
    "composite_sha256",
}
ARTIFACT_KEYS = {"id", "path", "size", "sha256"}
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
WINDOWS_RESERVED = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *("COM%d" % index for index in range(1, 10)),
    *("LPT%d" % index for index in range(1, 10)),
}
WINDOWS_FORBIDDEN = set('<>:"|?*')


class CompositeError(ValueError):
    """A deterministic manifest or artifact violated the contract."""


class _StrictArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise CompositeError("argument-error:" + message)


def _reject_duplicate_pairs(pairs: List[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CompositeError("duplicate-json-key:" + key)
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise CompositeError("non-finite-json-number:" + value)


def _strict_object(value: Any, label: str, keys: set) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise CompositeError("not-object:" + label)
    missing = keys - set(value)
    unknown = set(value) - keys
    if missing:
        raise CompositeError(
            "missing-keys:%s:%s" % (label, ",".join(sorted(missing)))
        )
    if unknown:
        raise CompositeError(
            "unknown-keys:%s:%s" % (label, ",".join(sorted(unknown)))
        )
    return value


def _strict_int(value: Any, label: str, minimum: int, maximum: int) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise CompositeError("invalid-integer:" + label)
    return value


def _valid_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise CompositeError("invalid-sha256:" + label)
    if value == "0" * 64:
        raise CompositeError("invalid-sha256:" + label)
    return value


def _canonical_path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise CompositeError("invalid-path:" + label)
    try:
        value.encode(PATH_ENCODING)
    except UnicodeEncodeError as exc:
        raise CompositeError("non-ascii-path:" + label) from exc
    if value.startswith("/") or "\\" in value:
        raise CompositeError("non-relative-posix-path:" + label)
    if any(ord(character) < 0x20 or ord(character) == 0x7F for character in value):
        raise CompositeError("control-character-in-path:" + label)

    parts = value.split("/")
    if any(not part for part in parts):
        raise CompositeError("empty-path-segment:" + label)
    for part in parts:
        if part in (".", ".."):
            raise CompositeError("non-canonical-path-segment:" + label)
        if part != part.strip() or part.endswith("."):
            raise CompositeError("non-portable-path-segment:" + label)
        if any(character in WINDOWS_FORBIDDEN for character in part):
            raise CompositeError("windows-forbidden-path-character:" + label)
        basename = part.split(".", 1)[0].upper()
        if basename in WINDOWS_RESERVED:
            raise CompositeError("windows-reserved-path-segment:" + label)
    return "/".join(parts)


def artifact_id(path: str) -> str:
    """Return the non-configurable semantic identity for a canonical path."""

    return "artifact:" + path


def _is_link_or_reparse(info: os.stat_result) -> bool:
    if stat.S_ISLNK(info.st_mode):
        return True
    attributes = getattr(info, "st_file_attributes", 0)
    reparse_flag = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return bool(attributes & reparse_flag)


def _absolute_lexical(path: Path) -> Path:
    return Path(os.path.abspath(str(path)))


def _assert_existing_chain_without_links(path: Path, label: str) -> None:
    absolute = _absolute_lexical(path)
    anchor = Path(absolute.anchor)
    current = anchor
    parts = absolute.parts[1:] if absolute.anchor else absolute.parts
    for part in parts:
        current = current / part
        try:
            info = os.lstat(str(current))
        except OSError as exc:
            raise CompositeError("path-component-unavailable:%s:%s" % (label, part)) from exc
        if _is_link_or_reparse(info):
            raise CompositeError("symlink-or-reparse-point:" + label)


def _assert_exact_relative_spelling(root: Path, relative: str) -> None:
    current = root
    for part in relative.split("/"):
        try:
            with os.scandir(str(current)) as entries:
                matches = [
                    entry.name
                    for entry in entries
                    if entry.name.casefold() == part.casefold()
                ]
        except OSError as exc:
            raise CompositeError(
                "path-component-unavailable:artifact:%s:%s" % (relative, part)
            ) from exc
        if not matches:
            raise CompositeError(
                "path-component-unavailable:artifact:%s:%s" % (relative, part)
            )
        if len(matches) != 1 or matches[0] != part:
            raise CompositeError("path-case-mismatch:" + relative)
        current = current / part


def _validated_root(root: Path) -> Tuple[Path, Path]:
    lexical = _absolute_lexical(root)
    _assert_existing_chain_without_links(lexical, "root")
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as exc:
        raise CompositeError("root-unavailable") from exc
    if not resolved.is_dir():
        raise CompositeError("root-not-directory")
    return lexical, resolved


def _artifact_path(root_lexical: Path, root_resolved: Path, relative: str) -> Path:
    candidate = root_lexical.joinpath(*relative.split("/"))
    _assert_existing_chain_without_links(candidate, "artifact:" + relative)
    _assert_exact_relative_spelling(root_lexical, relative)
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise CompositeError("path-escapes-root:" + relative) from exc
    except OSError as exc:
        raise CompositeError("artifact-unavailable:" + relative) from exc
    return candidate


def _same_file_identity(left: os.stat_result, right: os.stat_result) -> bool:
    if left.st_dev != right.st_dev:
        return False
    if left.st_ino and right.st_ino and left.st_ino != right.st_ino:
        return False
    return True


def _scan_artifact(
    root_lexical: Path,
    root_resolved: Path,
    relative: str,
    composite: "hashlib._Hash",
) -> Tuple[int, str, Path, os.stat_result]:
    candidate = _artifact_path(root_lexical, root_resolved, relative)
    try:
        before = os.lstat(str(candidate))
    except OSError as exc:
        raise CompositeError("artifact-stat-failed:" + relative) from exc
    if _is_link_or_reparse(before) or not stat.S_ISREG(before.st_mode):
        raise CompositeError("artifact-not-regular-file:" + relative)
    size = _strict_int(before.st_size, "artifact-size:" + relative, 0, MAX_U64)

    composite.update(relative.encode(PATH_ENCODING))
    composite.update(b"\0")
    composite.update(size.to_bytes(8, "big"))
    artifact_hash = hashlib.sha256()
    byte_count = 0
    try:
        with candidate.open("rb") as stream:
            opened = os.fstat(stream.fileno())
            if not stat.S_ISREG(opened.st_mode) or not _same_file_identity(before, opened):
                raise CompositeError("artifact-changed-before-read:" + relative)
            while True:
                chunk = stream.read(CHUNK_SIZE)
                if not chunk:
                    break
                byte_count += len(chunk)
                artifact_hash.update(chunk)
                composite.update(chunk)
    except CompositeError:
        raise
    except OSError as exc:
        raise CompositeError("artifact-read-failed:" + relative) from exc

    try:
        after = os.lstat(str(candidate))
    except OSError as exc:
        raise CompositeError("artifact-stat-failed-after-read:" + relative) from exc
    if (
        _is_link_or_reparse(after)
        or not _same_file_identity(before, after)
        or byte_count != size
        or after.st_size != before.st_size
        or getattr(after, "st_mtime_ns", None) != getattr(before, "st_mtime_ns", None)
    ):
        raise CompositeError("artifact-changed-during-read:" + relative)
    return size, artifact_hash.hexdigest(), candidate, after


def _physical_identity(path: Path, info: os.stat_result) -> Tuple[Any, ...]:
    if info.st_ino:
        return ("inode", info.st_dev, info.st_ino)
    return ("resolved-path", str(path.resolve(strict=True)).casefold())


def _verify_final_snapshots(
    snapshots: Sequence[Tuple[str, Path, os.stat_result]]
) -> None:
    for relative, candidate, before in snapshots:
        try:
            after = os.lstat(str(candidate))
        except OSError as exc:
            raise CompositeError("artifact-stat-failed-final:" + relative) from exc
        if (
            _is_link_or_reparse(after)
            or not stat.S_ISREG(after.st_mode)
            or not _same_file_identity(before, after)
            or after.st_size != before.st_size
            or getattr(after, "st_mtime_ns", None)
            != getattr(before, "st_mtime_ns", None)
        ):
            raise CompositeError("artifact-changed-after-scan:" + relative)


def _prepare_paths(paths: Iterable[str]) -> List[str]:
    canonical = [
        _canonical_path(path, "artifact[%d].path" % index)
        for index, path in enumerate(paths)
    ]
    if not canonical:
        raise CompositeError("artifact-list-empty")
    seen = set()
    for path in canonical:
        folded = path.casefold()
        if folded in seen:
            raise CompositeError("duplicate-or-case-aliased-path:" + path)
        seen.add(folded)
    return sorted(canonical)


def _canonical_manifest_bytes(document: Dict[str, Any]) -> bytes:
    return (
        json.dumps(document, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("ascii")


def _build_document(root: Path, paths: Iterable[str]) -> Dict[str, Any]:
    canonical_paths = _prepare_paths(paths)
    root_lexical, root_resolved = _validated_root(root)
    composite = hashlib.sha256()
    artifacts = []
    snapshots = []
    physical_identities = set()
    for relative in canonical_paths:
        size, digest, candidate, snapshot = _scan_artifact(
            root_lexical, root_resolved, relative, composite
        )
        physical_identity = _physical_identity(candidate, snapshot)
        if physical_identity in physical_identities:
            raise CompositeError("duplicate-physical-file:" + relative)
        physical_identities.add(physical_identity)
        snapshots.append((relative, candidate, snapshot))
        artifacts.append(
            {
                "id": artifact_id(relative),
                "path": relative,
                "sha256": digest,
                "size": size,
            }
        )
    _verify_final_snapshots(snapshots)
    return {
        "algorithm": ALGORITHM,
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
        "composite_sha256": composite.hexdigest(),
        "contract_version": CONTRACT_VERSION,
        "path_encoding": PATH_ENCODING,
        "schema": SCHEMA,
    }


def _reject_manifest_as_artifact(
    root: Path, manifest: Path, paths: Sequence[str]
) -> None:
    root_lexical, root_resolved = _validated_root(root)
    manifest_lexical = _absolute_lexical(manifest)
    try:
        relative = manifest_lexical.relative_to(root_lexical).as_posix()
    except ValueError:
        return
    try:
        relative = _canonical_path(relative, "manifest-relative-path")
    except CompositeError:
        return
    if relative.casefold() in {path.casefold() for path in paths}:
        raise CompositeError("manifest-cannot-be-an-artifact")
    # Resolve the root here so aliases through a different lexical spelling do
    # not evade the check when the output already exists.
    if manifest_lexical.exists():
        try:
            resolved = manifest_lexical.resolve(strict=True)
            resolved.relative_to(root_resolved)
        except (OSError, ValueError):
            return
        for path in paths:
            candidate = _artifact_path(root_lexical, root_resolved, path)
            if resolved == candidate.resolve(strict=True):
                raise CompositeError("manifest-cannot-be-an-artifact")


def _write_atomic(path: Path, payload: bytes) -> None:
    target = _absolute_lexical(path)
    parent = target.parent
    try:
        parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise CompositeError("manifest-parent-create-failed") from exc
    _assert_existing_chain_without_links(parent, "manifest-parent")
    if target.exists() or target.is_symlink():
        try:
            info = os.lstat(str(target))
        except OSError as exc:
            raise CompositeError("manifest-stat-failed") from exc
        if _is_link_or_reparse(info) or not stat.S_ISREG(info.st_mode):
            raise CompositeError("manifest-target-not-regular-file")

    descriptor = -1
    temporary_name = ""
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix="." + target.name + ".", suffix=".tmp", dir=str(parent)
        )
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, str(target))
        temporary_name = ""
        if target.read_bytes() != payload:
            raise CompositeError("manifest-post-write-byte-mismatch")
    except CompositeError:
        raise
    except OSError as exc:
        raise CompositeError("manifest-write-failed") from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary_name:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass


def generate_manifest(root: Path, manifest: Path, paths: Iterable[str]) -> str:
    """Generate a canonical manifest and return its composite SHA-256."""

    canonical_paths = _prepare_paths(paths)
    _reject_manifest_as_artifact(root, manifest, canonical_paths)
    document = _build_document(root, canonical_paths)
    _write_atomic(manifest, _canonical_manifest_bytes(document))
    return document["composite_sha256"]


def _load_manifest(path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    try:
        manifest = _absolute_lexical(path)
        _assert_existing_chain_without_links(manifest, "manifest")
        info = os.lstat(str(manifest))
        if _is_link_or_reparse(info) or not stat.S_ISREG(info.st_mode):
            raise CompositeError("manifest-not-regular-file")
        payload = manifest.read_bytes()
        text = payload.decode("ascii")
        document = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=_reject_json_constant,
        )
    except CompositeError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CompositeError("manifest-json-invalid") from exc

    root_object = _strict_object(document, "root", ROOT_KEYS)
    if root_object["schema"] != SCHEMA:
        raise CompositeError("schema-mismatch")
    if (
        type(root_object["contract_version"]) is not int
        or root_object["contract_version"] != CONTRACT_VERSION
    ):
        raise CompositeError("contract-version-mismatch")
    if root_object["algorithm"] != ALGORITHM:
        raise CompositeError("algorithm-mismatch")
    if root_object["path_encoding"] != PATH_ENCODING:
        raise CompositeError("path-encoding-mismatch")
    count = _strict_int(root_object["artifact_count"], "artifact_count", 1, MAX_U64)
    if not isinstance(root_object["artifacts"], list):
        raise CompositeError("artifacts-not-array")
    if len(root_object["artifacts"]) != count:
        raise CompositeError("artifact-count-mismatch")

    artifacts: List[Dict[str, Any]] = []
    seen_paths = set()
    seen_ids = set()
    previous_path = None
    for index, value in enumerate(root_object["artifacts"]):
        label = "artifacts[%d]" % index
        item = _strict_object(value, label, ARTIFACT_KEYS)
        relative = _canonical_path(item["path"], label + ".path")
        expected_id = artifact_id(relative)
        if item["id"] != expected_id:
            raise CompositeError("artifact-id-path-binding-mismatch:" + relative)
        if previous_path is not None and relative <= previous_path:
            raise CompositeError("artifacts-not-strictly-path-sorted")
        previous_path = relative
        folded_path = relative.casefold()
        folded_id = expected_id.casefold()
        if folded_path in seen_paths:
            raise CompositeError("duplicate-or-case-aliased-path:" + relative)
        if folded_id in seen_ids:
            raise CompositeError("duplicate-or-case-aliased-id:" + expected_id)
        seen_paths.add(folded_path)
        seen_ids.add(folded_id)
        size = _strict_int(item["size"], label + ".size", 0, MAX_U64)
        digest = _valid_sha256(item["sha256"], label + ".sha256")
        artifacts.append(
            {"id": expected_id, "path": relative, "size": size, "sha256": digest}
        )

    _valid_sha256(root_object["composite_sha256"], "composite_sha256")
    if _canonical_manifest_bytes(root_object) != payload:
        raise CompositeError("manifest-not-canonical-json")
    return root_object, artifacts


def verify_manifest(root: Path, manifest: Path, expected_composite: str) -> str:
    """Verify a manifest against live files and one externally frozen SHA."""

    expected = _valid_sha256(expected_composite, "expected-composite")
    document, artifacts = _load_manifest(manifest)
    embedded = document["composite_sha256"]
    if embedded != expected:
        raise CompositeError("expected-composite-mismatch")

    root_lexical, root_resolved = _validated_root(root)
    composite = hashlib.sha256()
    snapshots = []
    physical_identities = set()
    for artifact in artifacts:
        relative = artifact["path"]
        size, digest, candidate, snapshot = _scan_artifact(
            root_lexical, root_resolved, relative, composite
        )
        physical_identity = _physical_identity(candidate, snapshot)
        if physical_identity in physical_identities:
            raise CompositeError("duplicate-physical-file:" + relative)
        physical_identities.add(physical_identity)
        snapshots.append((relative, candidate, snapshot))
        if size != artifact["size"]:
            raise CompositeError("artifact-size-mismatch:" + relative)
        if digest != artifact["sha256"]:
            raise CompositeError("artifact-sha256-mismatch:" + relative)
    _verify_final_snapshots(snapshots)
    actual = composite.hexdigest()
    if actual != embedded:
        raise CompositeError("computed-composite-mismatch")
    return actual


def _parser() -> argparse.ArgumentParser:
    parser = _StrictArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate = subparsers.add_parser(
        "generate", help="write a canonical composite manifest"
    )
    generate.add_argument("--root", type=Path, default=Path("."))
    generate.add_argument("--manifest", type=Path, required=True)
    generate.add_argument("--artifact", action="append", required=True)
    generate.add_argument("--print-composite", action="store_true")

    verify = subparsers.add_parser(
        "verify", help="verify live artifacts against a frozen composite"
    )
    verify.add_argument("--root", type=Path, default=Path("."))
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--expected-composite", required=True)
    verify.add_argument("--print-composite", action="store_true")
    return parser


def _emit_failure(reason: str) -> None:
    print(
        json.dumps(
            {"reason": reason, "schema": RESULT_SCHEMA, "status": "FAIL"},
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        ),
        file=sys.stderr,
    )


def main(argv: Sequence[str] = None) -> int:
    try:
        args = _parser().parse_args(argv)
        if args.command == "generate":
            digest = generate_manifest(args.root, args.manifest, args.artifact)
        else:
            digest = verify_manifest(
                args.root, args.manifest, args.expected_composite
            )
        if args.print_composite:
            print(digest)
        return 0
    except CompositeError as exc:
        _emit_failure(str(exc))
        return 2
    except Exception as exc:
        _emit_failure("unexpected:%s:%s" % (type(exc).__name__, exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
