#!/usr/bin/env python3
# Convert TinyDEM mesh.off (with 1/2/3/4-vertex faces) to ASCII STL.
# Handles special 1-vertex (sphere) and 2-vertex (capsule) faces using Rmesh from tinydem.cpp.

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


def _read_rmesh(tinydem_cpp: Path):
    try:
        text = tinydem_cpp.read_text(encoding="utf-8")
    except OSError:
        return None
    m = re.search(r"Rmesh\s*\[\]\s*=\s*\{([^}]+)\}", text)
    if not m:
        return None
    nums = []
    for tok in m.group(1).split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            nums.append(float(tok))
        except ValueError:
            pass
    return nums if nums else None


def _vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_mul(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a):
    return math.sqrt(_dot(a, a))


def _normalize(a):
    n = _norm(a)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def _facet_normal(a, b, c):
    return _normalize(_cross(_vec_sub(b, a), _vec_sub(c, a)))


def _orthonormal_basis(w):
    # w must be normalized
    if abs(w[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)
    u = _normalize(_cross(w, ref))
    v = _cross(w, u)
    return u, v


def _add_triangle(tris, a, b, c):
    tris.append((a, b, c))


def _add_quad(tris, a, b, c, d):
    _add_triangle(tris, a, b, c)
    _add_triangle(tris, a, c, d)


def _add_sphere(tris, center, r, nseg=24, nstack=12):
    # Full sphere for simplicity
    for i in range(nstack):
        phi0 = math.pi * i / nstack
        phi1 = math.pi * (i + 1) / nstack
        for j in range(nseg):
            th0 = 2 * math.pi * j / nseg
            th1 = 2 * math.pi * (j + 1) / nseg
            p00 = (
                center[0] + r * math.sin(phi0) * math.cos(th0),
                center[1] + r * math.sin(phi0) * math.sin(th0),
                center[2] + r * math.cos(phi0),
            )
            p01 = (
                center[0] + r * math.sin(phi0) * math.cos(th1),
                center[1] + r * math.sin(phi0) * math.sin(th1),
                center[2] + r * math.cos(phi0),
            )
            p10 = (
                center[0] + r * math.sin(phi1) * math.cos(th0),
                center[1] + r * math.sin(phi1) * math.sin(th0),
                center[2] + r * math.cos(phi1),
            )
            p11 = (
                center[0] + r * math.sin(phi1) * math.cos(th1),
                center[1] + r * math.sin(phi1) * math.sin(th1),
                center[2] + r * math.cos(phi1),
            )
            _add_quad(tris, p00, p10, p11, p01)


def _add_capsule(tris, p0, p1, r, nseg=24, nstack=12):
    axis = _vec_sub(p1, p0)
    L = _norm(axis)
    if L == 0:
        _add_sphere(tris, p0, r, nseg=nseg, nstack=nstack)
        return
    w = _vec_mul(axis, 1.0 / L)
    u, v = _orthonormal_basis(w)

    # Cylinder
    ring0 = []
    ring1 = []
    for j in range(nseg):
        th = 2 * math.pi * j / nseg
        radial = _vec_add(_vec_mul(u, math.cos(th)), _vec_mul(v, math.sin(th)))
        ring0.append(_vec_add(p0, _vec_mul(radial, r)))
        ring1.append(_vec_add(p1, _vec_mul(radial, r)))
    for j in range(nseg):
        a = ring0[j]
        b = ring1[j]
        c = ring1[(j + 1) % nseg]
        d = ring0[(j + 1) % nseg]
        _add_quad(tris, a, b, c, d)

    # Hemispheres: approximate by full spheres (simpler & robust)
    _add_sphere(tris, p0, r, nseg=nseg, nstack=nstack)
    _add_sphere(tris, p1, r, nseg=nseg, nstack=nstack)


def _read_off(path: Path):
    raw = path.read_text(encoding="utf-8").splitlines()
    lines = []
    for ln in raw:
        ln = ln.split("#", 1)[0].strip()
        if ln:
            lines.append(ln)
    if not lines or lines[0] != "OFF":
        raise ValueError("Not an OFF file")
    counts = lines[1].split()
    npv, nfaces = int(counts[0]), int(counts[1])
    verts = []
    idx = 2
    for _ in range(npv):
        x, y, z = map(float, lines[idx].split())
        verts.append((x, y, z))
        idx += 1
    faces = []
    for _ in range(nfaces):
        parts = lines[idx].split()
        idx += 1
        n = int(parts[0])
        vidx = list(map(int, parts[1:1 + n]))
        face_index = int(parts[1 + n]) if len(parts) > 1 + n else 0
        faces.append((n, vidx, face_index))
    return verts, faces


def _write_stl(path: Path, tris):
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("solid tinydem\n")
        for a, b, c in tris:
            n = _facet_normal(a, b, c)
            f.write(f"  facet normal {n[0]} {n[1]} {n[2]}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {a[0]} {a[1]} {a[2]}\n")
            f.write(f"      vertex {b[0]} {b[1]} {b[2]}\n")
            f.write(f"      vertex {c[0]} {c[1]} {c[2]}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write("endsolid tinydem\n")


def main(argv):
    base = Path(__file__).resolve().parent
    off_path = base / "mesh.off"
    out_path = base / "mesh.stl"

    if len(argv) >= 2:
        off_path = Path(argv[1]).resolve()
    if len(argv) >= 3:
        out_path = Path(argv[2]).resolve()

    if not off_path.exists():
        print(f"OFF not found: {off_path}", file=sys.stderr)
        return 1

    rmesh = _read_rmesh(base / "tinydem.cpp")
    if rmesh is None:
        rmesh = [0.002, 0.02, 0.01]

    verts, faces = _read_off(off_path)
    tris = []

    for n, vidx, fidx in faces:
        r = rmesh[fidx] if 0 <= fidx < len(rmesh) else rmesh[0]
        if n == 1:
            _add_sphere(tris, verts[vidx[0]], r)
        elif n == 2:
            _add_capsule(tris, verts[vidx[0]], verts[vidx[1]], r)
        elif n == 3:
            a, b, c = (verts[i] for i in vidx)
            _add_triangle(tris, a, b, c)
        elif n == 4:
            a, b, c, d = (verts[i] for i in vidx)
            _add_quad(tris, a, b, c, d)
        else:
            # simple fan triangulation for n>4
            a = verts[vidx[0]]
            for i in range(1, n - 1):
                b = verts[vidx[i]]
                c = verts[vidx[i + 1]]
                _add_triangle(tris, a, b, c)

    _write_stl(out_path, tris)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
