"""Connectome-driven neural controller with the PUBLISHED spinal swim CPG.

Architecture (two layers, both documented):

1. CONNECTOME LAYER — the 2,616-node sensory-motor-pathway connectome
   (data/processed/connectome.npz) as a leaky firing-rate network with
   spike-frequency adaptation.  It receives (a) the tonic reticulospinal
   locomotor command (G_CMD onto descending units), (b) the lateral-line
   reflex (G_SENS, neuromast/afferent units, side-signed by segment), and
   (c) the Mauthner-like startle burst.  Its DESCENDING-pool output is sent
   to layer 2 as a slow gate envelope (TAU_GATE): the documented
   reticulospinal function of setting swim vigor.  A faster sensory-driven
   component differentially biases L vs R drive (heading stabilization; the
   lateral line of the larva projects via the posterior lateral line nerve
   to hindbrain command centers — Metcalfe et al. 2008).

2. SPINAL CPG LAYER — the published 8-population larval-zebrafish swim
   circuit, ported verbatim into neuromechfish/spine.py (Wandler et al.
   2025 eLife; cell-type-specific connectivity rules, fast/slow modules,
   distance-scaled axonal delays, tonic drive on E and I populations).  The
   connectome gate multiplies the published tonic drive per side
   (fastdrivemod/slowdrivemod = the example_run.py values 2.0 / 0.5).

Motor output: the ipsilateral E (V2a) population rate of each circuit
segment drives that segment's muscle (Roussel et al. 2021: MN firing →
muscle → body angle), low-pass filtered (TAU_ACT), resampled from the 30
biological somites of the circuit to the 25 mechanical hemisegments of
neuromechfish/body.py.

Nothing in the gait is scripted: beat rhythm, left/right alternation, and
the rostrocaudal traveling wave all emerge from the published circuit, and
vigor/steering from the connectome command and lateral-line feedback.
"""
# The previous version grafted a Matsuoka half-center pair per body segment
# (W_PAIR/G_PAIR/CHAIN/GRAD, a caudal amplitude ramp, a tanh differential
# readout).  All of it is removed: the published circuit supplies rhythm,
# alternation and the wave from documented cell-type-resolved wiring.

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from neuromechfish import spine as spine_mod
from neuromechfish.spine import TSTEP_S

TAU_M = 0.015          # membrane / rate time constant [s]
TAU_E = 0.006          # excitatory synaptic filter [s]
TAU_I = 0.060          # inhibitory synaptic filter [s]
G = 0.35               # global synaptic gain
E_I_BALANCE = 0.45     # relative strength of inhibitory synapses (gain knob)

G_CMD = 1.4            # tonic command drive onto descending units [a.u.]
G_ADAPT = 2.0          # spike-frequency adaptation strength [a.u.]
TAU_ADAPT = 0.15       # adaptation time constant [s] -> sets beat period
G_SENS = 0.25          # lateral-line feedback gain (tanh-scaled): at this
                       # level the reflex STABILIZES heading; at 0.5 the loop
                       # amplifies yaw (the fish wobbles), at 0 it swims badly
V_SCALE = 0.05         # body lateral velocity normalization [m/s]

# descending command -> published tonic drive gate (DRIVE coupling in
# spine.py; normalized so the nominal G_CMD reproduces the published
# example drive of fastdrivemod=2.0 / slowdrivemod=0.5 at gate ≈ 1)
G_GATE = 1.65           # envelope -> drive-scale gain (normalized: gate env
                        # ~0.59 at nominal G_CMD drives the published circuit
                        # at its example drive, scale ~= 1.0 — measured)
TAU_GATE = 0.05        # reticulospinal gate envelope time constant [s]
TAU_SENSE = 0.012      # lateral-line component of the gate [s]


def saturate(x: np.ndarray) -> np.ndarray:
    p = np.maximum(x, 0.0)
    return p / (1.0 + 0.5 * p)


class SMPNetwork:
    """Recurrent rate network over the zebrafish sensory-motor pathway."""

    def __init__(self, npz_path: str, seed: int = 0):
        d = np.load(npz_path, allow_pickle=True)
        self.names = d["names"].astype(str)
        self.W = sp.csr_matrix(
            (d["data"].astype(np.float32), d["indices"], d["indptr"]),
            shape=tuple(d["shape"]),
        )
        self.n = self.W.shape[0]
        self.cats = d["category"].astype(str)
        self.segments = d["segment"].astype(int)
        self.rng = np.random.default_rng(seed)

        outdeg = np.asarray(self.W.sum(axis=1)).ravel()
        scale = G / np.sqrt(np.maximum(outdeg, 1.0))
        signs = self._edge_signs()

        # build signed weight matrices explicitly from COO triplets (scipy's
        # chained sparse.multiply with dense columns silently zeroes data)
        coo = self.W.tocoo()
        src, dst = coo.row, coo.col
        w_edge = scale[src] * signs[src]
        is_inhib = signs[src] < 0
        self.W_E = sp.coo_matrix(
            (np.where(is_inhib, 0.0, w_edge), (src, dst)), shape=self.W.shape
        ).tocsr()
        self.W_I = sp.coo_matrix(
            (np.where(is_inhib, E_I_BALANCE * w_edge, 0.0), (src, dst)),
            shape=self.W.shape).tocsr()
        self.W_E.eliminate_zeros()
        self.W_I.eliminate_zeros()
        # store transposed: postsynaptic propagation is over target columns
        self.W_E = self.W_E.T.tocsr()
        self.W_I = self.W_I.T.tocsr()
        # inhibition/driver balance: connectome has ~2.3x more inhibitory than
        # excitatory synapses onto motor pools; scale inhibitory strength so
        # driven equilibrium sits mid-range (E_I_BALANCE, a gain knob).

        # index helpers by category
        self.idx = {c: np.where(self.cats == c)[0] for c in np.unique(self.cats)}
        self.des = self.idx["descending"]
        self.sens_units = np.concatenate([self.idx["neuromast"], self.idx["sensory"]])
        self.mn = np.concatenate([self.idx["motorn_primary"],
                                  self.idx["motorn_secondary"]])

        # per-side list of motor-neuron indices per segment 0..24
        self.mn_lists = [[[] for _ in range(25)] for _ in range(2)]
        for i in self.mn:
            side = 0 if "_L_" in self.names[i] else 1
            self.mn_lists[side][self.segments[i] - 1].append(i)

        # lateral-line units: every neuromast/afferent gets a body segment and
        # side sign (L +1 / R -1); count = 15 neuromasts, 20 afferents per side
        self.sens_seg = np.zeros(len(self.names), np.int32)
        self.sens_sign = np.zeros(len(self.names), np.float32)
        counts = {"neuromast": 15, "sensory": 20}
        for c, cnt in counts.items():
            for i in self.idx[c]:
                name = self.names[i]
                # names like Neuromast_L7 / Sens_R12
                m = name.rsplit("_", 1)[-1]
                side_sign = 1.0 if m[0] == "L" else -1.0
                num = int(m[1:])
                self.sens_seg[i] = 1 + int(round((num - 1) * 22.0 / (cnt - 1)))  # 2..24
                self.sens_sign[i] = side_sign

        self.x = 0.02 * self.rng.standard_normal(self.n)
        self.Ie = np.zeros(self.n)
        self.Ii = np.zeros(self.n)
        self.fat = np.zeros(self.n)          # spike-frequency adaptation state
        self.actL = np.zeros(25)
        self.actR = np.zeros(25)
        self.t = 0.0

        # ---------- descending command split per side + gate envelopes ----------
        self.des_L = [i for i in self.des if "_L_" in self.names[i]]
        self.des_R = [i for i in self.des if not ("_L_" in self.names[i])]
        self._gate = np.zeros(2)             # slow reticulospinal envelope
        self._sense = np.zeros(2)            # fast lateral-line component

        # ---------- the published spinal swim circuit (CPG layer) ----------
        self.spine = spine_mod.ZebrafishSpine()

    def _edge_signs(self) -> np.ndarray:
        """Per-**edge** sign array over the COO representation of ``self.W``.

        Family sign per presynaptic cell, with the commissural-inhibition
        override: a spinal interneuron that targets the *opposite* side acts
        inhibitory (mutual half-center inhibition); its ipsilateral targets
        keep the family sign.
        """
        base = assign_signs(self.names, self.cats)
        coo = self.W.tocoo()
        src, dst = coo.row, coo.col
        s = base[src].copy()
        cross = np.array(["_L_" in self.names[i] for i in src]) != \
            np.array(["_L_" in self.names[i] for i in dst])
        mask = cross & (self.cats[src] == "spinal_in")
        s[mask] = -1.0
        return s

    # ------------------------------------------------------------- dynamics
    def step(self, dt: float, I_ext: np.ndarray):
        y = saturate(self.x)                        # synaptic output, unfatigued
        self.Ie += (-self.Ie + self.W_E @ y) * dt / TAU_E
        self.Ii += (-self.Ii + self.W_I @ y) * dt / TAU_I
        # Matsuoka adaptation: slow self-inhibition on the DRIVE — each cell
        # fatigues its own excitability on the scale of TAU_ADAPT, which turns
        # the symmetric mutual-inhibition equilibrium into a self-sustaining
        # left/right alternating oscillator.
        self.fat += (-self.fat + G_ADAPT * y) * dt / TAU_ADAPT
        self.x = (self.x + (-self.x + self.Ie - self.Ii - self.fat + I_ext) * dt / TAU_M)
        self.t += dt

        # ----- connectome -> published tonic drive (spine.py DRIVE coupling)
        env = (y[self.des_L].mean(), y[self.des_R].mean())
        self._gate += (np.array(env) - self._gate) * dt / TAU_GATE
        # lateral-line fast component: differential L vs R (steering/stabilize)
        sense_diff = y[self.des_L].mean() - y[self.des_R].mean()
        self._sense += (sense_diff - self._sense) * dt / TAU_SENSE
        # gate: nominal G_CMD reproduces the published example drive (=1.0)
        scale = np.clip(G_GATE * self._gate
                        + 0.1 * G_SENS * self._sense / (V_SCALE * env[0]),
                        0.15, 1.6)
        # the published circuit: run its fixed 0.1 ms substep grid ON the real
        # time.  A physics dt finer than 0.1 ms needs a fractional accumulator
        # (n_sub = round(dt/1e-4) would floor to 1 and dilate the CPG 2x,
        # halving the beat frequency); coarser dt still fans out correctly.
        self._acc = getattr(self, "_acc", 0.0) + dt
        while self._acc >= TSTEP_S - 1e-12:
            self._acc -= TSTEP_S
            self.spine.substep(scale)

    def readout(self, dt: float) -> np.ndarray:
        """Published circuit's E (V2a) pools -> filtered bilateral muscles.

        See spine.py for the documented MUSCLE coupling (Roussel et al.
        2021: MN firing → muscle output → body angle).
        """
        act = self.spine.muscles(dt)               # (2, 25)
        self.actL, self.actR = act[0], act[1]
        return act

    # --------------------------------------------------------------- inputs
    def command_input(self) -> np.ndarray:
        """Descending locomotor command: tonic reticulospinal drive only.

        The beat, its L/R alternation and the head->tail traveling wave all
        emerge from network dynamics — the connectome's global oscillation
        paces the grafted segmental half-centers; nothing is time-scripted.
        """
        I = np.zeros(self.n)
        I[self.des] = G_CMD            # tonic reticulospinal drive
        return I

    def sensory_input(self, v_y_seg: np.ndarray) -> np.ndarray:
        """Lateral-line drive from per-segment lateral velocity (25,)."""
        I = np.zeros(self.n)
        drive = G_SENS * np.tanh(v_y_seg / V_SCALE)
        r = self.sens_units
        I[r] = self.sens_sign[r] * drive[self.sens_seg[r]]
        return I

    def startle_input(self, t: float, t0: float, dur: float = 0.02) -> np.ndarray:
        """Brief Mauthner-like escape burst onto all descending units."""
        I = np.zeros(self.n)
        if t0 <= t < t0 + dur:
            I[self.des] = 4.0
        return I


def assign_signs(names: np.ndarray, cats: np.ndarray) -> np.ndarray:
    """Presynaptic E/I assignment by cell family (documented modeling choice)."""
    signs = np.ones(len(names), np.float32)
    for i, (name, cat) in enumerate(zip(names, cats)):
        if cat == "spinal_in":
            # segment-family cells: numbers 1..4 excitatory, 5..14 inhibitory
            signs[i] = -1.0 if int(name.rsplit("_", 1)[-1]) >= 5 else 1.0
    return signs
