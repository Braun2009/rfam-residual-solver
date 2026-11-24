#!/usr/bin/env python3
"""Extract the six Lam 2010 supplement tables into canonical CSV."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


SOURCE_BYTE_COUNT = 115712
SOURCE_SHA256 = "154b690336f9bd2aa6a56c417386d620d5dc360e5e996a8b41900238bd17f977"
CSV_FIELDS = [
    "table_id",
    "row_in_table",
    "facility_mode",
    "t5_k",
    "p5_atm",
    "dp5_dt_percent_per_ms",
    "dp5_dt_status",
    "idt_ms",
]
EXPECTED_ROWS = [17, 18, 19, 18, 10, 13]
EXPECTED_MODES = [
    "conventional-pressure-rise",
    "near-constant-volume",
    "conventional-pressure-rise",
    "near-constant-volume",
    "conventional-pressure-rise",
    "near-constant-volume",
]
EXPECTED_COLUMNS = [4, 3, 4, 3, 4, 3]
EXPECTED_TITLE = (
    "Shock Tube Ignition Delay Time Measurements in Propane/O2/Argon "
    "Mixtures at Near-Constant-Volume Conditions"
)
EXPECTED_EXPERIMENT_DESCRIPTION = (
    "The following tables include experimental ignition delay time data for a "
    "mixture of 0.8% C3H8/ 8% O2/ Ar at nominal pressures of 6, 24 and 60 atm, "
    "under both conventional shock tube operation (with dP5/dt ~ 1-7%/ms) and "
    "near-constant-volume operation (with dP5/dt ~ 0%/ms)."
)
EXPECTED_HEADERS = {
    3: [
        "Post Reflected-Shock Temperature, T5 [K]",
        "Post Reflected-Shock Pressure, P5 [atm]",
        "Ignition Delay Time [milliseconds]",
    ],
    4: [
        "Post Reflected-Shock Temperature, T5 [K]",
        "Post Reflected-Shock Pressure, P5 [atm]",
        "dP5/dt [%/ms]",
        "Ignition Delay Time [milliseconds]",
    ],
}
EXPECTED_CAPTIONS = [
    "Table S1: Post-shock conditions (P5 = 5.5-7.2 atm) and measured ignition data with dP5/dt = 1.2-3.5%/ms.",
    "Table S2: Post-shock conditions (P5 = 5.3-6.8 atm) and measured ignition data with dP5/dt ~ 0%/ms.",
    "Table S3: Post-shock conditions (P5 = 21-29 atm) and measured ignition data with dP5/dt = 5.9-8.0%/ms.",
    "Table S4: Post-shock conditions (P5 = 23-31 atm) and measured ignition data with dP5/dt ~ 0%/ms.",
    "Table S5: Post-shock conditions (P5 = 55-67 atm) and measured ignition data with dP5/dt = 6.2-8.5%/ms.",
    "Table S6: Post-shock conditions (P5 = 54-68 atm) and measured ignition data with dP5/dt ~ 0%/ms.",
]

_INTEGER_RE = re.compile(r"^[1-9][0-9]*$")
_DECIMAL_RE = re.compile(r"^(?:0|[1-9][0-9]*)(?:\.[0-9]+)?$")


class GenerationError(ValueError):
    """A source, Word table, or output contract violation."""


def _clean(value: str) -> str:
    return value.replace("\r\x07", "").replace("\x07", "").strip()


def _read_source(path: Path) -> bytes:
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise GenerationError("source-read-failed:%s" % exc) from exc
    if len(payload) != SOURCE_BYTE_COUNT:
        raise GenerationError("source-byte-count-mismatch")
    if hashlib.sha256(payload).hexdigest() != SOURCE_SHA256:
        raise GenerationError("source-sha256-mismatch")
    return payload


def _validate_number(value: str, integer: bool, label: str) -> None:
    pattern = _INTEGER_RE if integer else _DECIMAL_RE
    if not pattern.fullmatch(value):
        raise GenerationError("table-number-invalid:%s" % label)
    if float(value) <= 0.0:
        raise GenerationError("table-number-nonpositive:%s" % label)


def _table_rows(document: Any) -> List[Dict[str, str]]:
    if int(document.Tables.Count) != 6:
        raise GenerationError("table-count-mismatch:%s" % document.Tables.Count)

    paragraphs = []
    for index in range(1, int(document.Paragraphs.Count) + 1):
        text = _clean(document.Paragraphs(index).Range.Text)
        if text:
            paragraphs.append(text)
    for label, expected in (
        ("title", EXPECTED_TITLE),
        ("experiment-description", EXPECTED_EXPERIMENT_DESCRIPTION),
    ):
        matches = [index for index, text in enumerate(paragraphs) if text == expected]
        if len(matches) != 1:
            raise GenerationError("%s-count-mismatch:%d" % (label, len(matches)))
    positions = []
    for caption in EXPECTED_CAPTIONS:
        matches = [index for index, text in enumerate(paragraphs) if text == caption]
        if len(matches) != 1:
            raise GenerationError("caption-count-mismatch:%s" % (len(positions) + 1))
        positions.append(matches[0])
    if positions != sorted(positions):
        raise GenerationError("caption-order-mismatch")
    title_position = paragraphs.index(EXPECTED_TITLE)
    description_position = paragraphs.index(EXPECTED_EXPERIMENT_DESCRIPTION)
    if not title_position < description_position < positions[0]:
        raise GenerationError("source-prose-order-mismatch")

    records: List[Dict[str, str]] = []
    for table_offset in range(6):
        table_number = table_offset + 1
        table = document.Tables(table_number)
        expected_rows = EXPECTED_ROWS[table_offset]
        expected_columns = EXPECTED_COLUMNS[table_offset]
        if int(table.Rows.Count) != expected_rows + 1:
            raise GenerationError("table-%d-row-count-mismatch" % table_number)
        if int(table.Columns.Count) != expected_columns:
            raise GenerationError("table-%d-column-count-mismatch" % table_number)
        headers = [
            _clean(table.Cell(1, column).Range.Text)
            for column in range(1, expected_columns + 1)
        ]
        if headers != EXPECTED_HEADERS[expected_columns]:
            raise GenerationError("table-%d-header-mismatch" % table_number)

        mode = EXPECTED_MODES[table_offset]
        for row_number in range(2, expected_rows + 2):
            source_values = [
                _clean(table.Cell(row_number, column).Range.Text)
                for column in range(1, expected_columns + 1)
            ]
            t5 = source_values[0]
            p5 = source_values[1]
            _validate_number(t5, True, "S%d-R%d-T5" % (table_number, row_number - 1))
            _validate_number(p5, False, "S%d-R%d-P5" % (table_number, row_number - 1))

            if expected_columns == 4:
                raw_dpdt = source_values[2]
                idt = source_values[3]
                if raw_dpdt == "N.A.":
                    dpdt = ""
                    dpdt_status = "not-available-cell"
                else:
                    _validate_number(
                        raw_dpdt,
                        False,
                        "S%d-R%d-dP5dt" % (table_number, row_number - 1),
                    )
                    dpdt = raw_dpdt
                    dpdt_status = "reported"
            else:
                idt = source_values[2]
                dpdt = ""
                dpdt_status = "not-tabulated-caption-approximately-zero"
            _validate_number(idt, False, "S%d-R%d-IDT" % (table_number, row_number - 1))
            records.append(
                {
                    "table_id": "S%d" % table_number,
                    "row_in_table": str(row_number - 1),
                    "facility_mode": mode,
                    "t5_k": t5,
                    "p5_atm": p5,
                    "dp5_dt_percent_per_ms": dpdt,
                    "dp5_dt_status": dpdt_status,
                    "idt_ms": idt,
                }
            )
    if len(records) != 95:
        raise GenerationError("total-row-count-mismatch:%d" % len(records))
    return records


def _canonical_csv(records: Sequence[Dict[str, str]]) -> bytes:
    import io

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    try:
        return stream.getvalue().encode("ascii")
    except UnicodeEncodeError as exc:
        raise GenerationError("csv-not-ascii") from exc


def extract(source: Path) -> bytes:
    _read_source(source)
    try:
        import pythoncom
        import win32com.client
    except ImportError as exc:
        raise GenerationError("word-com-runtime-unavailable") from exc

    word = None
    document = None
    pythoncom.CoInitialize()
    try:
        word = win32com.client.DispatchEx("Word.Application")
        word.Visible = False
        word.DisplayAlerts = 0
        document = word.Documents.Open(str(source.resolve()), False, True)
        records = _table_rows(document)
        return _canonical_csv(records)
    except GenerationError:
        raise
    except Exception as exc:
        raise GenerationError("word-com-extraction-failed:%s" % type(exc).__name__) from exc
    finally:
        if document is not None:
            try:
                document.Close(False)
            except Exception:
                pass
        if word is not None:
            try:
                word.Quit()
            except Exception:
                pass
        pythoncom.CoUninitialize()


def _write_atomic(path: Path, payload: bytes, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise GenerationError("output-exists-use-overwrite")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=str(path.parent),
            delete=False,
        ) as stream:
            temporary = stream.name
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, str(path))
        temporary = None
    except OSError as exc:
        raise GenerationError("output-write-failed:%s" % exc) from exc
    finally:
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    try:
        args = parser.parse_args(argv)
        payload = extract(args.source)
        _write_atomic(args.output, payload, args.overwrite)
        print(
            "R2_LAM2010_CSV_GENERATED tables=6 rows=95 bytes=%d sha256=%s"
            % (len(payload), hashlib.sha256(payload).hexdigest())
        )
        return 0
    except GenerationError as exc:
        print(
            json.dumps(
                {
                    "schema": "rfam.r2-lam2010-intake-generation/v1",
                    "status": "FAIL",
                    "reason": str(exc),
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2
    except Exception as exc:  # pragma: no cover - fail-closed boundary
        print(
            json.dumps(
                {
                    "schema": "rfam.r2-lam2010-intake-generation/v1",
                    "status": "FAIL",
                    "reason": "unexpected-%s" % type(exc).__name__,
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
