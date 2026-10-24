#!/usr/bin/env python
"""Build the zebrafish sensory-motor-pathway connectome.

Inputs
------
_work/Network_Creation/SMP_generated.dot   (produced by scripts/make_smp_dot.m)
_work/Network_Creation/*.mat               (rsMap, dendriteMap, SpineList, ...)

Outputs (data/processed/)
--------------------------
connectome.npz   CSR adjacency (binary), names, per-node metadata arrays
nodes.csv        idx, name, category, side, region, segment, cell_type
edges.csv        src_idx, dst_idx, src_name, dst_name (unique directed pairs)

Also ports ``updateMatrixDed.m`` (the post-processing step of the original
makeSMP.m): tops up sensory -> RS wiring to 25 outgoing connections per
sensory cell, rewires the neuromast layer, and trims/grows the RS-internal
block to exactly 70,000 connections.
"""

import re
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SMP_DIR = ROOT / "_work" / "Network_Creation"
DOT = SMP_DIR / "SMP_generated.dot"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

# Quoted node declarations look like `"Name";` but `label = "..."` must not
# be mistaken for a node declaration.
TOKEN_RE = re.compile(
    r'"([^"]+)"\s*->\s*"([^"]+)"'          # edge
    r'|(?<!label = )"([^"]+)"\s*;'          # node declaration
)
SPINAL_RE = re.compile(r"^(\d+)(In|Pri|Sec)_([LR])_(\d+)$")
BRAIN_RE = re.compile(r"^(Mid|Ro1|Ro2|Ro3|Mi1|Mi2|Mi3|Ca|T)_([ML])([LR])_(\d+)$")


def parse_dot(path: Path):
    """Single ordered scan: names in first-appearance order (LL, RS, spine)."""
    text = path.read_text()
    names: list[str] = []
    index: dict[str, int] = {}
    edges: list[tuple[int, int]] = []

    def node(n: str) -> int:
        i = index.get(n)
        if i is None:
            i = index[n] = len(names)
            names.append(n)
        return i

    for m in TOKEN_RE.finditer(text):
        if m.group(1) is not None:
            edges.append((node(m.group(1)), node(m.group(2))))
        else:
            node(m.group(3))
    return names, edges


def classify(name: str) -> dict:
    if name.startswith("Neuromast_"):
        return dict(category="neuromast", side=name.split("_")[1][0],
                    region="", segment=-1, cell_type="neuromast")
    if name.startswith("Sens_"):
        return dict(category="sensory", side=name.split("_")[1][0],
                    region="", segment=-1, cell_type="sensory")
    if name.startswith("Des_"):
        parts = name.split("_")            # Des_<type>_<side>_<n>
        return dict(category="descending", side=parts[2],
                    region="", segment=-1, cell_type=parts[1])
    m = SPINAL_RE.match(name)
    if m:
        seg, kind, side, idx = int(m.group(1)), m.group(2), m.group(3), int(m.group(4))
        cat = {"In": "spinal_in", "Pri": "motorn_primary", "Sec": "motorn_secondary"}[kind]
        return dict(category=cat, side=side, region="spine", segment=seg,
                    cell_type=f"{kind}{idx}")
    m = BRAIN_RE.match(name)
    if m:
        region, div, side, idx = m.groups()
        return dict(category="brain_in", side=side, region=region,
                    segment=-1, cell_type=f"{div}{idx}")
    return dict(category="unknown", side="", region="", segment=-1, cell_type="")


# --------------------------------------------------------------------------
# Port of updateMatrixDed.m
# --------------------------------------------------------------------------
SENS_TO_RS = np.array([
    [1, 2, 3, 4, 0, 0, 0],
    [1, 2, 3, 4, 5, 6, 0],
    [3, 4, 5, 6, 7, 8, 0],
    [5, 6, 7, 8, 9, 10, 11],
    [7, 8, 9, 10, 11, 12, 13],
    [9, 10, 11, 12, 13, 14, 15],
    [12, 13, 14, 15, 16, 17, 0],
    [14, 15, 16, 17, 18, 19, 20],
    [16, 17, 18, 19, 20, 0, 0],
])
RS_SEC = [("Mid_M", "Mid_L"), ("Ro1_M", "Ro1_L"), ("Ro2_M", "Ro2_L"),
          ("Ro3_M", "Ro3_L"), ("Mi1_M", "Mi1_L"), ("Mi2_M", "Mi2_L"),
          ("Mi3_M", "Mi3_L"), ("Ca_M", "Ca_L"), ("T_M", "T_L")]
NUM_CONS_PER_CELL = np.array([16, 16, 19, 19, 10, 10, 7, 7, 8, 8,
                              8, 10, 10, 6, 6, 12, 12, 11, 11, 11])
NUM_NEURO_CONN = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5])
NEW_MAX_CONNEC = 25
RS_TARGET_NNZ = 70_000


def update_matrix_ded(adj: sp.csr_matrix, names: list[str],
                      rng: np.random.Generator) -> sp.csr_matrix:
    """Port of updateMatrixDed.m operating on an edge set (binary graph)."""
    n = len(names)
    index = {nm: i for i, nm in enumerate(names)}
    is_ll = np.array([nm.startswith(("Neuromast_", "Sens_")) for nm in names])
    is_spine = np.array([SPINAL_RE.match(nm) is not None for nm in names])
    is_rs = ~(is_ll | is_spine)
    rs_idx = np.flatnonzero(is_rs)
    t_idx = np.flatnonzero(np.array([nm.startswith("T_") for nm in names]))

    present: set[tuple[int, int]] = set(zip(*map(tuple, adj.nonzero())))
    prefix_pool: dict[str, np.ndarray] = {}

    def section_pool(prefix: str) -> np.ndarray:
        """Brain interneurons of one RS section prefix (e.g. 'Mid_ML')."""
        pool = prefix_pool.get(prefix)
        if pool is None:
            pool = np.array([i for i, nm in enumerate(names)
                             if nm.startswith(prefix + "_") and is_rs[i]])
            prefix_pool[prefix] = pool
        return pool

    # --- 1. top up sensory -> RS interneuron connections to 25 per cell -----
    for side in ("L", "R"):
        for sens in range(1, 21):
            rows, _ = np.where(SENS_TO_RS == sens)
            targets = []
            for r in rows:
                targets.append(section_pool(RS_SEC[r][0] + side))
                targets.append(section_pool(RS_SEC[r][1] + side))
            rs_int = np.concatenate(targets) if targets else np.array([], int)
            if rs_int.size:
                rs_int = rs_int[rng.permutation(rs_int.size)]
            tmp_sens = index[f"Sens_{side}{sens}"]
            for j in range(NEW_MAX_CONNEC - NUM_CONS_PER_CELL[sens - 1]):
                present.add((tmp_sens, int(rs_int[j])))

    # --- 2. zero the neuromast -> sensory block, then rewire ----------------
    present = {e for e in present if not (is_ll[e[0]] and is_ll[e[1]])}

    for side in ("L", "R"):
        remain = NUM_NEURO_CONN[rng.permutation(15)].copy()
        for sens in range(1, 21):
            sens_i = index[f"Sens_{side}{sens}"]
            nz = np.flatnonzero(remain > 0)
            if nz.size == 0:
                # MATLAB isempty branch: force-connect Neuromast_15
                nm_i = index[f"Neuromast_{side}15"]
                present.add((nm_i, sens_i))
                prev = index.get(f"Sens_{side}{sens - 1}")
                if prev is not None:
                    present.discard((nm_i, prev))
                continue
            start_p = int(nz[0])
            take_two = not (start_p == 14 or
                            (start_p + 1 < 15 and remain[start_p + 1] == 0))
            picks = [start_p] + ([start_p + 1] if take_two else [])
            for p in picks:
                if p < 15 and remain[p] > 0:
                    edge = (index[f"Neuromast_{side}{p + 1}"], sens_i)
                    if edge not in present:      # MATLAB: prints ERROR on dup
                        present.add(edge)
                        remain[p] -= 1

        # fill leftover capacity (MATLAB "Did NOT connect everything" path)
        while remain.sum() > 0:
            nm_p = int(np.flatnonzero(remain > 0)[0])
            nm_i = index[f"Neuromast_{side}{nm_p + 1}"]
            placed = 0
            for s in range(1, 21):
                if placed >= remain[nm_p] or remain[nm_p] == 0:
                    break
                edge = (nm_i, index[f"Sens_{side}{s}"])
                if edge not in present:
                    present.add(edge)
                    remain[nm_p] -= 1
                    placed += 1
            if placed == 0:
                remain[nm_p] = 0              # avoid infinite loops

    # --- 3. resize the RS-internal block to exactly 70,000 ------------------
    rs_internal = sum(1 for i, j in present if is_rs[i] and is_rs[j])
    conn_rmv = rs_internal - RS_TARGET_NNZ
    in_t = set(t_idx.tolist())
    t_present = [(i, j) for i, j in present if i in in_t and j in in_t]
    if conn_rmv > 0:
        order = rng.permutation(len(t_present))
        for k in order[:conn_rmv]:
            present.discard(t_present[int(k)])
    elif conn_rmv < 0:
        need = -conn_rmv
        missing = [(a, b) for a in t_idx for b in t_idx
                   if a != b and (a, b) not in present]
        order = rng.permutation(len(missing))
        for k in order[:need]:
            present.add(missing[int(k)])

    rows = np.fromiter((e[0] for e in present), np.int64, len(present))
    cols = np.fromiter((e[1] for e in present), np.int64, len(present))
    return sp.coo_matrix((np.ones(len(present), np.uint8), (rows, cols)),
                         shape=(n, n)).tocsr()


# --------------------------------------------------------------------------
def main():
    if not DOT.exists():
        sys.exit(f"{DOT} missing - run scripts/make_smp_dot.m first")

    names, edges = parse_dot(DOT)
    print(f"parsed {len(names)} nodes, {len(edges)} edge statements from {DOT.name}")

    cats = [classify(nm) for nm in names]
    uniq, counts = np.unique([c["category"] for c in cats], return_counts=True)
    print("node categories:", dict(zip(uniq.tolist(), counts.tolist())))
    if len(names) != 2616:
        print(f"WARNING: expected 2616 nodes (70 LL + 1196 RS + 1350 spine), got {len(names)}")

    rows = np.array([e[0] for e in edges])
    cols = np.array([e[1] for e in edges])
    data = np.ones(len(edges), np.uint8)
    adj = sp.coo_matrix((data, (rows, cols)), shape=(len(names),) * 2).tocsr()
    adj.data[:] = 1                                 # binary, like drawDot output
    adj.eliminate_zeros()
    print(f"unique directed edges after parse: {adj.nnz}")

    rng = np.random.default_rng(42)
    adj = update_matrix_ded(adj, names, rng)
    adj.eliminate_zeros()
    print(f"edges after updateMatrixDed: {adj.nnz}")

    # ---- category-level stats ------------------------------------------
    def block(pred_s, pred_t):
        s = np.array([bool(pred_s(nm)) for nm in names])
        t = np.array([bool(pred_t(nm)) for nm in names])
        return adj[np.ix_(np.flatnonzero(s), np.flatnonzero(t))].nnz

    stats = {
        "neuromast->sensory": block(lambda n: n.startswith("Neuromast_"),
                                    lambda n: n.startswith("Sens_")),
        "sensory->RS": block(lambda n: n.startswith("Sens_"),
                             lambda n: not n.startswith(("Sens_", "Neuromast_"))
                             and not SPINAL_RE.match(n)),
        "RS internal": block(lambda n: not n.startswith(("Sens_", "Neuromast_"))
                             and not SPINAL_RE.match(n),
                             lambda n: not n.startswith(("Sens_", "Neuromast_"))
                             and not SPINAL_RE.match(n)),
        "RS->spine": block(lambda n: not n.startswith(("Sens_", "Neuromast_"))
                           and not SPINAL_RE.match(n),
                           lambda n: SPINAL_RE.match(n) is not None),
        "spine internal": block(lambda n: SPINAL_RE.match(n) is not None,
                                lambda n: SPINAL_RE.match(n) is not None),
        "total": adj.nnz,
    }
    readme = {
        "neuromast->sensory": 70, "sensory->RS": 1000, "RS internal": 70000,
        "RS->spine": 34800, "spine internal": 61244, "total": 167114,
    }
    print(f"{'block':<20}{'built':>10}{'README':>10}")
    for k in stats:
        print(f"{k:<20}{stats[k]:>10}{readme[k]:>10}")

    # ---- save ------------------------------------------------------------
    np.savez_compressed(
        OUT / "connectome.npz",
        indptr=adj.indptr.astype(np.int32),
        indices=adj.indices.astype(np.int32),
        data=np.ones(adj.nnz, np.uint8),
        shape=np.array(adj.shape),
        names=np.array(names, dtype=object),
        category=np.array([c["category"] for c in cats], dtype=object),
        side=np.array([c["side"] for c in cats], dtype=object),
        region=np.array([c["region"] for c in cats], dtype=object),
        segment=np.array([c["segment"] for c in cats], dtype=np.int16),
        cell_type=np.array([c["cell_type"] for c in cats], dtype=object),
    )
    with open(OUT / "nodes.csv", "w") as fh:
        fh.write("idx,name,category,side,region,segment,cell_type\n")
        for i, nm in enumerate(names):
            c = cats[i]
            fh.write(f"{i},{nm},{c['category']},{c['side']},{c['region']},"
                     f"{c['segment']},{c['cell_type']}\n")
    src, dst = adj.nonzero()
    with open(OUT / "edges.csv", "w") as fh:
        fh.write("src_idx,dst_idx,src_name,dst_name\n")
        for s, d in zip(src.tolist(), dst.tolist()):
            fh.write(f"{s},{d},{names[s]},{names[d]}\n")

    print(f"\nwrote {OUT}/connectome.npz, nodes.csv ({len(names)}), edges.csv ({adj.nnz})")


if __name__ == "__main__":
    main()
