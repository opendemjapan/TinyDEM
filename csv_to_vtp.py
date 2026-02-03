#!/usr/bin/env python3
# TinyDEM CSV -> VTK XML PolyData (.vtp) converter
# Usage: csv_to_vtp.py <input.csv|directory> [output.vtp|output_dir]

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


def _write_array(fh, values, comps, indent="        ", per_line=6):
    # values: flat list of numbers (as strings)
    count = 0
    fh.write(indent)
    for v in values:
        fh.write(v)
        fh.write(" ")
        count += 1
        if count % (per_line * comps) == 0:
            fh.write("\n")
            fh.write(indent)
    if not fh.tell() == 0:
        fh.write("\n")


def _write_vtp(out_path: Path, fieldnames, rows, *, no_vectors=False, no_quat=False, only_r=False):
    n = len(rows)
    float_type = "Float64"

    # Points
    points = []
    for r in rows:
        points.extend([r["x"], r["y"], r["z"]])

    # Vector arrays
    vector_defs = {}
    if not no_vectors:
        if {"vx", "vy", "vz"}.issubset(fieldnames):
            vector_defs["v"] = ("vx", "vy", "vz")
        if {"wx", "wy", "wz"}.issubset(fieldnames):
            vector_defs["w"] = ("wx", "wy", "wz")

    # Quaternion as 4-component array if present
    has_q = (not no_quat) and {"q0", "q1", "q2", "q3"}.issubset(fieldnames)

    excluded = {"x", "y", "z"}
    if no_vectors:
        excluded.update({"vx", "vy", "vz", "wx", "wy", "wz"})
    if no_quat:
        excluded.update({"q0", "q1", "q2", "q3"})
    for cols in vector_defs.values():
        excluded.update(cols)
    if has_q:
        excluded.update({"q0", "q1", "q2", "q3"})

    scalar_fields = [name for name in fieldnames if name not in excluded]
    if only_r:
        scalar_fields = ["R"] if "R" in fieldnames else []

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="\n") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <PolyData>\n")
        f.write(f'    <Piece NumberOfPoints="{n}" NumberOfVerts="{n}" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n')

        # PointData (ParaView expects this before geometry in some builds)
        f.write("      <PointData>\n")

        # Vectors
        for name, cols in vector_defs.items():
            f.write(f'        <DataArray type="{float_type}" Name="{name}" NumberOfComponents="3" format="ascii">\n')
            if n > 0:
                vals = []
                for r in rows:
                    vals.extend([r[cols[0]], r[cols[1]], r[cols[2]]])
                _write_array(f, vals, 3)
            f.write("        </DataArray>\n")

        # Quaternion (optional)
        if has_q:
            f.write(f'        <DataArray type="{float_type}" Name="q" NumberOfComponents="4" format="ascii">\n')
            if n > 0:
                vals = []
                for r in rows:
                    vals.extend([r["q0"], r["q1"], r["q2"], r["q3"]])
                _write_array(f, vals, 4)
            f.write("        </DataArray>\n")

        # Scalars
        for name in scalar_fields:
            f.write(f'        <DataArray type="{float_type}" Name="{name}" NumberOfComponents="1" format="ascii">\n')
            if n > 0:
                vals = [r[name] for r in rows]
                _write_array(f, vals, 1)
            f.write("        </DataArray>\n")

        f.write("      </PointData>\n")

        # CellData (empty)
        f.write("      <CellData>\n")
        f.write("      </CellData>\n")

        # Points
        f.write("      <Points>\n")
        f.write(f'        <DataArray type="{float_type}" Name="Points" NumberOfComponents="3" format="ascii">\n')
        if n > 0:
            _write_array(f, points, 3)
        f.write("        </DataArray>\n")
        f.write("      </Points>\n")

        # Verts (one vertex per point)
        f.write("      <Verts>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        if n > 0:
            _write_array(f, [str(i) for i in range(n)], 1)
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        if n > 0:
            _write_array(f, [str(i + 1) for i in range(n)], 1)
        f.write("        </DataArray>\n")
        f.write("      </Verts>\n")

        # Lines (empty)
        f.write("      <Lines>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write("      </Lines>\n")

        # Strips (empty)
        f.write("      <Strips>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write("      </Strips>\n")

        # Polys (empty)
        f.write("      <Polys>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        f.write("        </DataArray>\n")
        f.write("      </Polys>\n")

        f.write("    </Piece>\n")
        f.write("  </PolyData>\n")
        f.write("</VTKFile>\n")


def _convert_file(in_path: Path, out_path: Path, *, no_vectors=False, no_quat=False, only_r=False):
    fieldnames, rows = _read_csv(in_path)
    _write_vtp(out_path, fieldnames, rows, no_vectors=no_vectors, no_quat=no_quat, only_r=only_r)


def _convert_dir(in_dir: Path, out_dir: Path | None, *, no_vectors=False, no_quat=False, only_r=False):
    csv_files = sorted(in_dir.glob("particles*.csv"))
    if not csv_files:
        raise FileNotFoundError("No particles*.csv found")
    for csv_path in csv_files:
        if out_dir is None:
            out_path = csv_path.with_suffix(".vtp")
        else:
            out_path = out_dir / (csv_path.stem + ".vtp")
        _convert_file(csv_path, out_path, no_vectors=no_vectors, no_quat=no_quat, only_r=only_r)


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print("Usage: csv_to_vtp.py [--no-vectors] [--no-quat] [--only-r] <input.csv|directory> [output.vtp|output_dir]")
        return 0
    flags = set()
    args = []
    for a in argv[1:]:
        if a.startswith("--"):
            flags.add(a)
        else:
            args.append(a)
    if not args:
        print("Usage: csv_to_vtp.py [--no-vectors] [--no-quat] [--only-r] <input.csv|directory> [output.vtp|output_dir]")
        return 0

    no_vectors = "--no-vectors" in flags
    no_quat = "--no-quat" in flags
    only_r = "--only-r" in flags

    in_path = Path(args[0]).resolve()
    if not in_path.exists():
        print(f"Input not found: {in_path}", file=sys.stderr)
        return 1

    if in_path.is_dir():
        out_dir = Path(args[1]).resolve() if len(args) >= 2 else None
        if out_dir is not None and out_dir.exists() and out_dir.is_file():
            print("Output must be a directory when input is a directory", file=sys.stderr)
            return 1
        _convert_dir(in_path, out_dir, no_vectors=no_vectors, no_quat=no_quat, only_r=only_r)
        return 0

    # input is file
    if len(args) >= 2:
        out_path = Path(args[1]).resolve()
    else:
        out_path = in_path.with_suffix(".vtp")
    _convert_file(in_path, out_path, no_vectors=no_vectors, no_quat=no_quat, only_r=only_r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
