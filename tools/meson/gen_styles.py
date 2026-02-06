#!/usr/bin/env python3
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path


def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return f.read()


def write_style(out_path: Path, token: str, header_glob: str) -> None:
    headers = sorted(glob.glob(header_glob))
    picked: list[str] = []
    for h in headers:
        txt = read_text(h)
        if token in txt:
            picked.append(os.path.basename(h))

    out_lines = [f'#include "{h}"\n' for h in picked]
    out_path.write_text("".join(out_lines), encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: gen_styles.py <srcdir> <outdir>", file=sys.stderr)
        return 2

    srcdir = Path(sys.argv[1]).resolve()
    outdir = Path(sys.argv[2]).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Mirrors src/Make.sh style() patterns
    rules = [
        ("style_app.h",     "APP_CLASS",     "app_*.h"),
        ("style_command.h", "COMMAND_CLASS", "*.h"),
        ("style_diag.h",    "DIAG_CLASS",    "diag_*.h"),
        ("style_dump.h",    "DUMP_CLASS",    "dump_*.h"),
        ("style_pair.h",    "PAIR_CLASS",    "pair_*.h"),
        ("style_region.h",  "REGION_CLASS",  "region_*.h"),
        ("style_solve.h",   "SOLVE_CLASS",   "solve_*.h"),
    ]

    for out_name, token, pat in rules:
        write_style(outdir / out_name, token, str(srcdir / pat))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
