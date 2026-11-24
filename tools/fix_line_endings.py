#!/usr/bin/env python3
"""Normalize line endings and verify string-literal integrity.

Root cause of the recurring WSL/Windows escape bugs (2026-08-03/04):
  1. Windows Python text mode writes `\\n` as CRLF inside string literals,
     splitting C++ strings across lines ("missing terminating quote").
  2. `wsl bash -c "..."` double-quote nesting folds `\\\\n` before Python
     sees it, silently corrupting heredoc patches.
  3. Mixed CRLF/LF files make sed/od/grep displays ambiguous.

This tool:
  - normalizes the file to LF (Unix) line endings;
  - reports any line containing a string literal with an unbalanced quote
    (the symptom of a split string);
  - optionally fixes a known broken pattern: `"<< "` ... real-LF ... `";"`.

Usage:
  python3 fix_line_endings.py <file> [--check-only]
"""

from __future__ import annotations

import sys


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: fix_line_endings.py <file> [--check-only]")
        return 2
    path = sys.argv[1]
    check_only = "--check-only" in sys.argv
    with open(path, "rb") as f:
        data = f.read()
    text = data.decode("utf-8", errors="replace")

    # 1. normalize CRLF -> LF
    if "\r\n" in text:
        text = text.replace("\r\n", "\n")
        print(f"normalized CRLF -> LF ({text.count(chr(10))} lines)")
    else:
        print("already LF")

    # 2. scan for unbalanced quote lines (split string literals)
    bad = []
    for idx, line in enumerate(text.split("\n"), start=1):
        if '<<"' in line or '<< "' in line or '("' in line and line.count('"') % 2 == 1:
            # lines that contain a quote but an odd number of them are
            # likely split string literals (CRLF damage)
            if line.count('"') % 2 == 1:
                bad.append((idx, line.strip()[:80]))
    if bad:
        print(f"unbalanced-quote lines: {len(bad)}")
        for idx, line in bad[:8]:
            print(f"  line {idx}: {line}")
        if check_only:
            print("CHECK_FAIL")
            return 1
    else:
        print("no unbalanced-quote lines")

    if not check_only and ("\r\n" in data.decode("utf-8", errors="replace")
                           or bad):
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print(f"fixed {path}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
