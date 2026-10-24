"""Reader for the Fish1 VAST `.vsvi` container descriptors.

The three Petkova_ZFish1_*.vsvi files are JSON descriptors written by VAST.
They do not contain voxel data themselves; they point at sharded,
neuroglancer-precomputed volumes on the public bucket
``storage.googleapis.com``:

  * Petkova_ZFish1_EM_Clahe23.vsvi        -> fish1-public/clahe_231218        (EM, uint8)
  * Petkova_ZFish1_Base_Segmentation.vsvi -> fish1-public/seg_241003          (supervoxels, uint64)
  * Petkova_ZFish1_Agglo_Segmentation.vsvi -> fish1-public/seg_241003_agg241003 (agglomerated, uint64)

This module parses a descriptor, derives the bucket prefix from
``SourceFileNameTemplate`` and streams individual chunks through
cloud-volume (``neuroglancer_uint64_sharded_v1`` shards, HTTP range
requests), so only requested ROIs are ever downloaded.

Coordinates: CloudVolume slicing is ``vol[x0:x1, y0:y1, z0:z1]``; returned
arrays come back in that same (x, y, z[, c]) order from the 12.x frontend.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]

EM_VSVI = PROJECT_ROOT / "Petkova_ZFish1_EM_Clahe23.vsvi"
BASE_SEG_VSVI = PROJECT_ROOT / "Petkova_ZFish1_Base_Segmentation.vsvi"
AGGLO_SEG_VSVI = PROJECT_ROOT / "Petkova_ZFish1_Agglo_Segmentation.vsvi"

# Documented example cell (Fish1 programmatic-access notebook), 8x8x30 nm voxels.
EXAMPLE_CELL_XYZ = (75765, 38366, 6682)


def load_descriptor(path: str | Path) -> dict:
    """Parse the header fields of a .vsvi descriptor.

    VAST .vsvi files are JSON-like but not strict JSON (the embedded
    ``neuroglancer_sharded`` block uses unquoted keys / colon syntax), so
    the handful of fields needed to locate the data are regex-extracted.
    """
    text = Path(path).read_text(errors="replace")

    def grab_str(key: str) -> str:
        m = re.search(rf'"{key}"\s*:\s*"([^"]*)"', text)
        if not m:
            raise KeyError(f"{key} not found in {path}")
        return m.group(1)

    out = {
        "ServerName": grab_str("ServerName"),
        "SourceFileNameTemplate": grab_str("SourceFileNameTemplate"),
        "SourceDataEncoding": grab_str("SourceDataEncoding"),
    }
    m = re.search(r'"ServerPort"\s*:\s*(\d+)', text)
    out["ServerPort"] = int(m.group(1)) if m else 443
    return out


def bucket_prefix(descriptor: dict) -> str:
    """``/fish1-public/clahe_231218/%s/`` -> ``/fish1-public/clahe_231218``."""
    tmpl = descriptor["SourceFileNameTemplate"]
    return tmpl.replace("%s", "").rstrip("/")


def http_url(descriptor: dict) -> str:
    host = descriptor["ServerName"]
    port = int(descriptor.get("ServerPort", 443))
    scheme = "https" if port == 443 else "http"
    return f"{scheme}://{host}{bucket_prefix(descriptor)}"


def volume(descriptor: dict | str | Path, mip: int = 0):
    """Open a CloudVolume for the volume referenced by a .vsvi descriptor."""
    from cloudvolume import CloudVolume

    if not isinstance(descriptor, dict):
        descriptor = load_descriptor(descriptor)
    return CloudVolume(f"precomputed://{http_url(descriptor)}", mip=mip,
                       progress=False, bounded=True)


def fetch(descriptor: dict | str | Path, mip: int,
          x: int, y: int, z: int, size: int | tuple[int, int, int]) -> np.ndarray:
    """Fetch an ROI centred on (x, y, z) (voxel coords at that mip).

    ``size`` is (sx, sy, sz) voxels; a scalar sizes all three axes.
    Returns an (sx, sy, sz) array.
    """
    if isinstance(size, int):
        size = (size, size, size)
    sx, sy, sz = size
    vol = volume(descriptor, mip)
    x0, x1 = x - sx // 2, x - sx // 2 + sx
    y0, y1 = y - sy // 2, y - sy // 2 + sy
    z0, z1 = z - sz // 2, z - sz // 2 + sz
    shape = tuple(int(s) for s in vol.shape[:3])
    for a, b, ext in ((x0, x1, shape[0]), (y0, y1, shape[1]), (z0, z1, shape[2])):
        if a < 0 or b > ext:
            raise ValueError(f"ROI [{a}:{b}] outside volume extent 0:{ext}")
    block = np.asarray(vol[x0:x1, y0:y1, z0:z1, 0])
    return block.reshape(sx, sy, sz)


def at_mip(xyz: tuple[int, int, int], from_res: tuple[int, int, int],
           to_res: tuple[int, int, int]) -> tuple[int, int, int]:
    """Rescale voxel coordinates between resolutions (z usually unchanged)."""
    return tuple(int(round(c * f / t)) for c, f, t in zip(xyz, from_res, to_res))
