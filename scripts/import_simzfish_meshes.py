"""Import anatomical mesh geometry for the NeuroMechFish body.

Reads the two biologically-derived mesh assets of the simZFish larval
zebrafish (Liu et al., Science Robotics 2025, EPFL BioRob, Apache-2.0)
and writes them into ``meshes/`` re-oriented for MuJoCo:

    simZFish proto axes:  x = lateral, y = up, z = rostro-caudal (head first)
    NeuroMechFish axes:   x = rostro-caudal, y = lateral, z = up

so the transform is a cyclic permutation (x, y, z) -> (z, x, y).
Vertices are re-centred on the geometric bbox centre and rescaled so the
head fits the model's ``L_HEAD`` and the fin matches the trunk scale.

Run:  .venv/bin/python scripts/import_simzfish_meshes.py
"""

from __future__ import annotations

import re
import struct
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
# vendored copies of the two source meshes (see meshes/README.md)
_SRC = _ROOT / "_thirdparty/simzfish"
_DST = _ROOT / "meshes"

# from protos/Zebrafish_4mm_demo.proto body chain (radial head + 6 segments
# of equal length): the built model is 10x a real 4-mm larva
L_HEAD = 1.4e-3          # [m]  == neuromechfish.body.L_HEAD
FIN_WIDTH = 5.2e-4       # [m]  dorsoventral spread of the caudal fin plate

MESHES = {
    # map the raw along-body span onto the model's head length
    "Zebrafish_head.stl": ("head.stl", L_HEAD, 2),      # axis 2 (z) -> L_HEAD
    "Tail fin.stl": ("fin.stl", FIN_WIDTH, 1),          # axis 1 (y) -> FIN_WIDTH
}


def parse_ascii_stl(path: Path) -> np.ndarray:
    txt = path.read_text()
    verts = re.findall(
        r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", txt)
    return np.array([tuple(map(float, v)) for v in verts], np.float64).reshape(-1, 3)


def write_binary_stl(path: Path, tris: np.ndarray) -> None:
    n = len(tris)
    with open(path, "wb") as fh:
        fh.write(b"NeuroMechFish anatomical mesh (simZFish, Apache-2.0)"[:80].ljust(80, b"\0"))
        fh.write(struct.pack("<I", n))
        for t in tris:
            a, b, c = t
            normal = np.cross(b - a, c - a)
            nl = np.linalg.norm(normal)
            normal = normal / nl if nl > 0 else np.zeros(3)
            fh.write(struct.pack("<12fH",
                                 *normal, *a, *b, *c, 0))


def transform(verts: np.ndarray, scale: float) -> np.ndarray:
    v = verts[:, [2, 0, 1]]              # cyclic perm (head-first check below)
    v = v * scale
    center = (v.min(0) + v.max(0)) / 2
    return v - center                    # recentre on bbox midpoint


def main() -> None:
    _DST.mkdir(exist_ok=True)
    for src_name, (dst_name, target_size, axis) in MESHES.items():
        src = _SRC / src_name
        if not src.exists():
            print(f"source mesh missing: {src}", file=sys.stderr)
            sys.exit(1)
        verts = parse_ascii_stl(src)
        raw_span = (verts.max(0) - verts.min(0))[axis]
        scale = target_size / raw_span
        tris = transform(verts, scale).reshape(-1, 3, 3)
        write_binary_stl(_DST / dst_name, tris)
        span = tris.reshape(-1, 3)
        p = span.max(0) - span.min(0)
        print(f"{dst_name}: {len(tris)} tris  "
              f"span x {p[0]*1e3:.3f} mm, y {p[1]*1e3:.3f} mm, z {p[2]*1e3:.3f} mm")


if __name__ == "__main__":
    main()
