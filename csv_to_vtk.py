#!/usr/bin/env python3
# TinyDEM CSV -> VTK (legacy ASCII) converter
# Usage: csv_to_vtk.py <input.csv|directory> [output.vtk|output_dir]

from __future__ import annotations

import csv
import sys
from pathlib import Path


def _read_csv(path: Path):
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV header missing")
        rows = list(reader)
        return reader.fieldnames, rows


def _write_vtk(out_path: Path, fieldnames, rows):
    n = len(rows)

    # Extract points
    points = []
    for r in rows:
        try:
            points.append((float(r["x"]), float(r["y"]), float(r["z"])))
        except KeyError as e:
            raise ValueError(f"Missing column: {e}") from e

    # Build vector arrays if possible
    vector_defs = {}
    if {"vx", "vy", "vz"}.issubset(fieldnames):
        vector_defs["v"] = ("vx", "vy", "vz")
    if {"wx", "wy", "wz"}.issubset(fieldnames):
        vector_defs["w"] = ("wx", "wy", "wz")

    excluded = {"x", "y", "z"}
    for cols in vector_defs.values():
        excluded.update(cols)

    scalar_fields = [name for name in fieldnames if name not in excluded]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="\n") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("TinyDEM particles\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        # Use double precision to avoid underflow issues with very small values.
        f.write(f"POINTS {n} double\n")
        for x, y, z in points:
            f.write(f"{x} {y} {z}\n")

        f.write(f"VERTICES {n} {n * 2}\n")
        for i in range(n):
            f.write(f"1 {i}\n")

        f.write(f"POINT_DATA {n}\n")

        # VECTORS
        for vec_name, cols in vector_defs.items():
            f.write(f"VECTORS {vec_name} double\n")
            for r in rows:
                f.write(f"{float(r[cols[0]])} {float(r[cols[1]])} {float(r[cols[2]])}\n")

        # SCALARS
        for name in scalar_fields:
            f.write(f"SCALARS {name} double 1\n")
            f.write("LOOKUP_TABLE default\n")
            for r in rows:
                f.write(f"{float(r[name])}\n")


def _convert_file(in_path: Path, out_path: Path):
    fieldnames, rows = _read_csv(in_path)
    _write_vtk(out_path, fieldnames, rows)


def _convert_dir(in_dir: Path, out_dir: Path | None):
    csv_files = sorted(in_dir.glob("particles*.csv"))
    if not csv_files:
        raise FileNotFoundError("No particles*.csv found")
    for csv_path in csv_files:
        if out_dir is None:
            out_path = csv_path.with_suffix(".vtk")
        else:
            out_path = out_dir / (csv_path.stem + ".vtk")
        _convert_file(csv_path, out_path)


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print("Usage: csv_to_vtk.py <input.csv|directory> [output.vtk|output_dir]")
        return 0

    in_path = Path(argv[1]).resolve()
    if not in_path.exists():
        print(f"Input not found: {in_path}", file=sys.stderr)
        return 1

    if in_path.is_dir():
        out_dir = Path(argv[2]).resolve() if len(argv) >= 3 else None
        if out_dir is not None and out_dir.exists() and out_dir.is_file():
            print("Output must be a directory when input is a directory", file=sys.stderr)
            return 1
        _convert_dir(in_path, out_dir)
        return 0

    # input is file
    if len(argv) >= 3:
        out_path = Path(argv[2]).resolve()
    else:
        out_path = in_path.with_suffix(".vtk")
    _convert_file(in_path, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
