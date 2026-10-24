"""Prototype: find sign rules / gains where the connectome self-oscillates.

Variants of the E/I assignment are applied on top of the same graph and
driven only by tonic descending drive (no scripted CPG).  We measure
motor-neuron (L/R) alternation, oscillation amplitude and frequency.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import scipy.sparse as sp

from neuromechfish import network as net


def rebuild(nm, signs_src):
    """Rebuild W_E / W_I with a per-edge sign array (signs is fn(src,dst)->array)."""
    coo = nm.W.tocoo()
    src, dst = coo.row, coo.col
    outdeg = np.asarray(nm.W.sum(axis=1)).ravel()
    scale = net.G / np.sqrt(np.maximum(outdeg, 1.0))
    s = signs_src(src, dst)
    w_edge = scale[src] * s
    is_inhib = s < 0
    WE = sp.coo_matrix((np.where(is_inhib, 0.0, w_edge), (src, dst)),
                       shape=nm.W.shape).tocsr().T.tocsr()
    WI = sp.coo_matrix((np.where(is_inhib, net.E_I_BALANCE * w_edge, (src, dst)),
                        ) if False else (np.where(is_inhib, net.E_I_BALANCE * w_edge, 0.0), (src, dst)),
                       shape=nm.W.shape).tocsr().T.tocsr()
    WE.eliminate_zeros(); WI.eliminate_zeros()
    nm.W_E, nm.W_I = WE, WI
    return nm


def rule_commissural_inhib(src, dst, nm):
    """Spinal interneurons inhibit the *opposite* side (CiD motif); ipsilateral
    projections keep the cell-family sign."""
    base = net.assign_signs(nm.names, nm.cats)
    s = base[src].copy()
    cross = ("_L_" in nm.names[src]) != ("_L_" in nm.names[dst])
    s[cross & (nm.cats[src] == "spinal_in")] = -1.0
    return s


def rule_all_commissural_inhib(src, dst, nm):
    """Every crossing edge is inhibitory, regardless of family."""
    base = net.assign_signs(nm.names, nm.cats)
    s = base[src].copy()
    cross = ("_L_" in nm.names[src]) != ("_L_" in nm.names[dst])
    s[cross] = -1.0
    return s


def run_with_adaptation(label, g_cmd=1.1, e_i=0.45, g=0.12, g_adapt=0.8,
                        tau_a=0.15, watch=4.0):
    """Commissural half-center + Matsuoka-style adaptation on all units."""
    net.G_CMD = g_cmd
    net.E_I_BALANCE = e_i
    net.G = g
    nm = net.SMPNetwork("data/processed/connectome.npz", seed=0)
    coo = nm.W.tocoo()
    s_arr = rule_commissural_inhib(coo.row, coo.col, nm)
    rebuild(nm, lambda s_, d_: s_arr)

    a_fat = np.zeros(nm.n)

    def step(dt, I_ext):
        nonlocal a_fat
        s_raw = net.saturate(nm.x)
        s_eff = np.maximum(0.0, s_raw - a_fat)          # fatigued firing rate
        nm.Ie += (-nm.Ie + nm.W_E @ s_eff) * dt / net.TAU_E
        nm.Ii += (-nm.Ii + nm.W_I @ s_eff) * dt / net.TAU_I
        nm.x += (-nm.x + nm.Ie - nm.Ii + I_ext) * dt / net.TAU_M
        a_fat += (-a_fat + g_adapt * np.maximum(0.0, s_raw - a_fat)) \
                 * dt / tau_a
        nm.t += dt
        return s_eff

    def readout(dt):
        s_eff = np.maximum(0.0, net.saturate(nm.x) - a_fat)
        rates = s_eff
        mn_act = np.zeros((2, 25))
        for side in (0, 1):
            for k in range(25):
                ids = nm.mn_lists[side][k]
                if ids:
                    mn_act[side, k] = rates[ids].mean()
        if not hasattr(nm, "_baseline"):
            nm._baseline = mn_act.copy()
        else:
            nm._baseline += (mn_act - nm._baseline) * dt / 1.0
        mn_act = (mn_act - nm._baseline) / net.MN_NORM + 0.5
        nm.actL += (np.clip(mn_act[0], 0, 1) - nm.actL) * dt / net.TAU_ACT
        nm.actR += (np.clip(mn_act[1], 0, 1) - nm.actR) * dt / net.TAU_ACT
        return np.stack([nm.actL, nm.actR])

    I_cmd = nm.command_input()
    for _ in range(2000):
        step(1e-3, I_cmd); nm.readout(1e-3)
    n_steps = int(watch / 1e-3)
    tr = np.zeros((2, 25, n_steps))
    for k in range(n_steps):
        step(1e-3, nm.command_input())
        tr[:, :, k] = readout(1e-3)

    aL, aR = tr[0], tr[1]
    aL5 = aL[4] - aL[4].mean(); aR5 = aR[4] - aR[4].mean()
    corr = np.corrcoef(aL5, aR5)[0, 1]
    f = np.fft.rfftfreq(n_steps, 1e-3)
    P = np.abs(np.fft.rfft(aL5))
    print(f"{label:32s} seg5 ptp {np.ptp(aL5):.3f}  L/R corr {corr:+.2f}  "
          f"peak {f[P.argmax()]:.2f} Hz  |act| {np.abs(aL).mean():.3f}")
    return tr


def run_variant(label, g_cmd=1.1, e_i=0.45, g=0.12, rule="default", watch=4.0):
    net.G_CMD = g_cmd
    net.E_I_BALANCE = e_i
    net.G = g
    nm = net.SMPNetwork("data/processed/connectome.npz", seed=0)
    if rule != "default":
        coo = nm.W.tocoo()
        s_arr = {"commissural": rule_commissural_inhib,
                 "all_cross": rule_all_commissural_inhib}[rule](
            coo.row, coo.col, nm)
        nm = rebuild(nm, lambda s_, d_: s_arr)

    def tonic(self):
        I = np.zeros(self.n)
        I[self.des] = net.G_CMD
        return I
    orig_cmd = net.SMPNetwork.command_input
    net.SMPNetwork.command_input = tonic
    try:
        I_cmd = nm.command_input()
        for _ in range(2000):
            nm.step(1e-3, I_cmd); nm.readout(1e-3)
        n_steps = int(watch / 1e-3)
        tr = np.zeros((2, 25, n_steps))
        for k in range(n_steps):
            nm.step(1e-3, nm.command_input())
            tr[:, :, k] = nm.readout(1e-3)
    finally:
        net.SMPNetwork.command_input = orig_cmd

    aL, aR = tr[0], tr[1]
    aL5 = aL[4] - aL[4].mean(); aR5 = aR[4] - aR[4].mean()
    corr = np.corrcoef(aL5, aR5)[0, 1]
    f = np.fft.rfftfreq(n_steps, 1e-3)
    P = np.abs(np.fft.rfft(aL5))
    anti = max(corr, -corr)
    print(f"{label:32s} seg5 ptp {np.ptp(aL5):.3f}  L/R corr {corr:+.2f}  "
          f"peak {f[P.argmax()]:.2f} Hz  |act| {np.abs(aL).mean():.3f}")


if __name__ == "__main__":
    for ga in (0.6, 0.8, 1.0, 1.3):
        run_with_adaptation(f"adapt {ga}", g_adapt=ga)
    run_with_adaptation("adapt 0.8 G .2", g_adapt=0.8, g=0.2)
    run_with_adaptation("adapt 1.0 E/I .75", g_adapt=1.0, e_i=0.75)
    run_with_adaptation("adapt 0.8 tau .25", g_adapt=0.8, tau_a=0.25)
    run_with_adaptation("adapt 1.0 cmd 1.6", g_adapt=1.0, g_cmd=1.6)
