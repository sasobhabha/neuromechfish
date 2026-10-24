"""Fetch an EM + segmentation sample from the Fish1 release using the
three Petkova_ZFish1_*.vsvi descriptors and export figures.

Outputs (outputs/fish1_sample/):
  * em_slice.png            - raw EM slice (16 nm)
  * seg_overlay.png         - EM with agglomerated-segment outlines
  * seg_panels.png          - EM | supervoxels | agglomerated segments
  * sample.npz              - raw blocks (EM, base seg, agglo seg)
  * summary.txt             - label counts etc.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neuromechfish import vsvi  # noqa: E402

OUT = vsvi.PROJECT_ROOT / "outputs" / "fish1_sample"
OUT.mkdir(parents=True, exist_ok=True)

# The documented example cell is in 8 nm voxels; seg volumes top out at 16 nm.
X, Y, Z = vsvi.at_mip(vsvi.EXAMPLE_CELL_XYZ, (8, 8, 30), (16, 16, 30))
SIZE = (192, 160, 16)  # x, y, z voxels @16 nm (asymmetric: asserts axis order)

print(f"Fetching ROI @16nm centre=({X},{Y},{Z}) size={SIZE} ...")
em = vsvi.fetch(vsvi.EM_VSVI, mip=1, x=X, y=Y, z=Z, size=SIZE)
base = vsvi.fetch(vsvi.BASE_SEG_VSVI, mip=0, x=X, y=Y, z=Z, size=SIZE)
agg = vsvi.fetch(vsvi.AGGLO_SEG_VSVI, mip=0, x=X, y=Y, z=Z, size=SIZE)

# Confirm array orientation: extents must map to (x, y, z)
assert em.shape == SIZE, f"unexpected EM block shape {em.shape} != {SIZE}"
print(f"EM {em.shape} range [{em.min()},{em.max()}] mean {em.mean():.1f}")
print(f"base seg: {len(np.unique(base))} labels, agglo seg: {len(np.unique(agg))} labels")

sl = em.shape[2] // 2  # middle z slice
em2d = em[:, :, sl]
base2d = base[:, :, sl]
agg2d = agg[:, :, sl]

# Pick slice with most labels for display if this one is empty
if len(np.unique(agg2d)) < 3:
    counts = [len(np.unique(agg[:, :, k])) for k in range(agg.shape[2])]
    sl = int(np.argmax(counts))
    em2d, base2d, agg2d = em[:, :, sl], base[:, :, sl], agg[:, :, sl]

fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(em2d.T, cmap="gray", origin="lower")
ax.set_title(f"Fish1 EM — z-slice {sl} @16 nm")
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT / "em_slice.png", dpi=150)
plt.close(fig)

# outlines of agglomerated segments over EM
def outlines(lab):
    edge = np.zeros(lab.shape, bool)
    edge[:, :-1] |= lab[:, :-1] != lab[:, 1:]
    edge[:, 1:] |= lab[:, :-1] != lab[:, 1:]
    edge[:-1, :] |= lab[:-1, :] != lab[1:, :]
    edge[1:, :] |= lab[:-1, :] != lab[1:, :]
    return edge

fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(em2d.T, cmap="gray", origin="lower")
edge = outlines(agg2d)
overlay = np.zeros((*edge.shape, 4))
overlay[edge] = (0.2, 1.0, 0.3, 0.9)
ax.imshow(overlay.transpose(1, 0, 2), origin="lower")
ax.set_title(f"EM + agglomerated segment outlines ({len(np.unique(agg2d))} segs)")
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT / "seg_overlay.png", dpi=150)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
axes[0].imshow(em2d.T, cmap="gray", origin="lower")
axes[0].set_title("EM (clahe_231218)")
for ax, data, title in (
    (axes[1], base2d, f"Supervoxels ({len(np.unique(base2d))})"),
    (axes[2], agg2d, f"Agglomerated ({len(np.unique(agg2d))})"),
):
    ids = np.unique(data)
    rng = np.random.default_rng(0)
    color = np.zeros((*data.shape, 4))
    for i, sid in enumerate(ids):
        if sid == 0:
            continue
        color[data == sid] = (*rng.random(3), 0.85)
    ax.imshow(color.transpose(1, 0, 2), origin="lower")
    ax.set_title(title)
for ax in axes:
    ax.axis("off")
fig.suptitle("Petkova_ZFish1_*.vsvi volumes — same ROI, three channels")
fig.tight_layout()
fig.savefig(OUT / "seg_panels.png", dpi=150)
plt.close(fig)

np.savez_compressed(OUT / "sample.npz", em=em, base_seg=base, agglo_seg=agg,
                    xyz_nm16=(X, Y, Z), size=SIZE)

summary = f"""Fish1 sample export
descriptor prefix EM    : {vsvi.bucket_prefix(vsvi.load_descriptor(vsvi.EM_VSVI))}
descriptor prefix base  : {vsvi.bucket_prefix(vsvi.load_descriptor(vsvi.BASE_SEG_VSVI))}
descriptor prefix agglo : {vsvi.bucket_prefix(vsvi.load_descriptor(vsvi.AGGLO_SEG_VSVI))}
ROI centre (16 nm vox): ({X}, {Y}, {Z})
ROI size (x,y,z)     : {SIZE}
EM block             : dtype={em.dtype} range=[{em.min()},{em.max()}]
supervoxel labels    : {len(np.unique(base))}
agglomerated labels  : {len(np.unique(agg))}
displayed z-slice    : {sl}
"""
(OUT / "summary.txt").write_text(summary)
print(summary)
print(f"Wrote {OUT}")
