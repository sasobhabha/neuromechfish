"""Published spinal swim circuit of the larval zebrafish, ported verbatim.

Neural model: the 8-population rate circuit of Wandler et al. 2025
("Coordinated spinal locomotor network dynamics emerge from
cell-type-specific connectivity patterns", eLife) — a strain of the
McLean-lab wiring rules for the developing zebrafish spinal cord.  Both the
population dynamics and the connectivity rules below are line-for-line the
published implementation (github.com/fwandler/ZFswim, ``parameters.py`` /
``connectome.py`` / ``snet.py``); nothing here is tuned beyond the
published defaults except the two documented couplings at the end.

Circuit (per spinal segment × side × population):
  E fast / E slow                 — descending excitatory (V2a/dIN cascade)
  I ipsi asc / des fast / slow    — ipsilateral ascending/descending inhibitors
  I contra fast / slow            — commissural (crossed) inhibitors

Rate dynamics (published `timestep_diffdelay`, forward Euler, tstep 0.1 ms):
    x' = sat( drive/tau_m + Cstrength * Σ W[n1,a1,p1,n2,a2,p2] r_delayed / tau_m
              − x_prev/tau_m ),   r_delayed = rate(t − (1+|n1−n2|)·delayts)
with membrane time constants fast 1 ms / slow 10 ms and axonal delays
fast 0.2 ms / slow 0.5 ms per segment of separation (published
`fastdelayts`=2 / `slowdelayts`=5 at tstep=0.1 ms).

Connectivity (published `connection`): descending ipsilateral excitatory
cascade with a rostral two-segment gap and a caudal projection cutoff
(E: Edes_cutoff_fast=16, slow=4 segments); ipsilateral inhibition gated by
the rostrocaudal phase rule (inhibit only between phasemin=0.3 and
phasemax=0.8 of the segment/separation phase); crossed inhibition with
Icutoff=13; fast/slow modules coupled by speedmix=0.5; global gain
Cstrength=0.5; tonic drive on BOTH E and I populations (published
fastdrivemod/slowdrivemod).  N=30 segments = the number of biological
somites of the larva.

Only two couplings are added, both documented physiology, no tuning knobs:
 1. DRIVE  — the connectome rate network (SMPNetwork) gates the published
    tonic drive multiplicatively per side (descending command sets swim
    vigor; asymmetric drive steers; startle → escape burst).
 2. MUSCLE — the ipsilateral E (V2a) population rate of each segment drives
    that segment's muscle (Roussel et al. 2021: MN firing → muscle output →
    body angle; V2a ablation abolishes swimming: Eklöf-Ljunggren et al.
    2012), low-pass filtered with the muscle activation time constant
    TAU_ACT = 30 ms.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# published parameters — verbatim from fwandler/ZFswim parameters.py, "8pop"
# ---------------------------------------------------------------------------
TSTEP_MS = 0.1          # published numerical timestep [ms]
TSTEP = TSTEP_MS        # integration step in PUBLISHED rate units (not SI!)
TSTEP_S = TSTEP_MS * 1e-3

PARAMS_8POP = {
    "modelname": "8pop",
    "Numsegs": 30,
    "Numtypes": 8,
    "typenames": ['I ipsi asc fast', 'I ipsi asc slow', 'I ipsi des fast',
                  'I ipsi des slow', 'I contra fast', 'I contra slow',
                  'E fast', 'E slow'],
    "rngseed": 27,
    "simtime": 600,
    "tstep": TSTEP_MS,
    "speedmix": 0.5,
    "fastdrivemod": 2.0,   # example_run.py drive (fdrive)
    "slowdrivemod": 0.5,   # example_run.py drive (sdrive)
    "fastmtc": 1.0,        # fast population membrane time constant [ms]
    "slowmtc": 10.0,       # slow population membrane time constant [ms]
    "timedelay_off": False,
    "fastdelayts": 2,      # axonal delay, fast pools [steps / seg dist]
    "slowdelayts": 5,      # axonal delay, slow pools [steps / seg dist]
    "eibalanced": False,
    "esumadd": 0,
    "phasemax": 0.8,
    "phasemin": 0.3,
    "Cstrength": 0.5,
    "Estrength": 0.5,
    "Istrength": 1.0,
    "drive": 1.0,
    "Iboost": 1.0,
    "cutoff": 30,
    "Icutoff": 13,
    "Edes_cutoff": 16,
    "firstEcutoff": 2,
    "Easc_cutoff": 0,
    "Espread": 1.0,
    "driveslope": 0.0,
    "Eablation": 1.0,
    "Iipsiablation": 1.0,
    "Iipsiascablation": 1.0,
    "Iipsidesablation": 1.0,
    "Icontraablation": 1.0,
    "Edes_cutoff_fast": 16,
    "Edes_cutoff_slow": 4,
}

MUSCLE_TAU = 0.004       # muscle activation time constant [s]. Bandwidth-
                         # matched to the published gait: first-order |H| >= 0.7
                         # at the circuit's 40 Hz beat needs tau <= ~4 ms
                         # (larval muscle follows 20-80 Hz bursts nearly 1:1,
                         # Buss & Drapeau 2000)
MUSCLE_RATE_SCALE = 0.7  # E-population rate giving ~saturated muscle act
N_BODY = 25              # mechanical hemisegments (neuromechfish.body)


# ---------------------------------------------------------------------------
# published connectivity rules — verbatim from fwandler/ZFswim connectome.py
# ---------------------------------------------------------------------------
def phase_to_conn(phase: float, phasemax: float, phasemin: float) -> float:
    """Published phase rule: inhibit if 'out of phase' (0.3..0.8), else 0."""
    if phase > phasemax:
        return 0.0
    elif phase < phasemin:
        return 0.0
    else:
        return -1.0


def connection(n1, a1, p1, n2, a2, p2, params):
    """Connection strength from source (n1,a1,p1) to target (n2,a2,p2).

    Verbatim port of `connectome.connection` in fwandler/ZFswim.
    """
    N = params["Numsegs"]
    speedmix = params["speedmix"]
    cutoff = params["cutoff"]
    phasemax = params["phasemax"]
    phasemin = params["phasemin"]

    segdist = np.abs(n1 - n2)
    if segdist > cutoff:
        return 0.0
    if [n1, a1, p1] == [n2, a2, p2]:
        return 0.0  # No self-connections

    if a1 == a2:
        altdiff = 0
    else:
        altdiff = 1

    if n1 < n2:
        descdiff = n2 - n1
        dphase = (descdiff / N) + (altdiff / 2)
        direction = 'desc'
    else:
        descdiff = n1 - n2
        dphase = (descdiff / N) + (altdiff / 2)
        dphase = 1 - dphase
        if descdiff == 0:
            direction = 'inseg'
        else:
            direction = 'asc'

    dphase = dphase % 1

    factor = 1.0
    conn = 0.0
    ablation = 1.0

    if 'E' in params["typenames"][p1]:
        ablation = params["Eablation"]
        if 'fast' in params["typenames"][p1]:
            Edes_cutoff = params["Edes_cutoff_fast"]
        elif 'slow' in params["typenames"][p1]:
            Edes_cutoff = params["Edes_cutoff_slow"]
        else:
            Edes_cutoff = params["Edes_cutoff"]
        if altdiff == 0:
            if (np.abs(n1 - n2) < Edes_cutoff) and (direction != 'asc') \
                    and (np.abs(n1 - n2) > params["firstEcutoff"]):
                conn = 1.0
            elif np.abs(n1 - n2) < params["Easc_cutoff"]:
                conn = 1.0
            else:
                conn = 0.0
        else:
            conn = 0.0

        if 'Espeedmix' in params:
            speedmix = params['Espeedmix']

        factor *= params["Estrength"]

    elif 'I' in params["typenames"][p1]:
        if np.abs(n1 - n2) > params["Icutoff"]:
            return 0.0
        conn = phase_to_conn(dphase, phasemax, phasemin)

        if 'ipsi' in params["typenames"][p1]:
            if 'asc' in params["typenames"][p1]:
                ablation = params["Iipsiascablation"]
                if direction == "desc":
                    conn = 0.0
            elif 'des' in params["typenames"][p1]:
                ablation = params["Iipsidesablation"]
                if direction == "asc":
                    conn = 0.0
            else:
                ablation = params["Iipsiablation"]

            if altdiff == 0:
                factor *= params["Istrength"]
            else:
                conn = 0.0

        elif 'contra' in params["typenames"][p1]:
            ablation = params["Icontraablation"]
            if altdiff == 1:
                factor *= params["Istrength"]
            else:
                conn = 0.0

        if 'Ispeedmix' in params:
            speedmix = params['Ispeedmix']

    if 'fast' in params["typenames"][p1] and 'fast' in params["typenames"][p2]:
        factor *= 1 - speedmix
    if 'fast' in params["typenames"][p1] and 'slow' in params["typenames"][p2]:
        factor *= speedmix
    if 'slow' in params["typenames"][p1] and 'slow' in params["typenames"][p2]:
        factor *= 1 - speedmix
    if 'slow' in params["typenames"][p1] and 'fast' in params["typenames"][p2]:
        factor *= speedmix

    return factor * conn * params["Cstrength"] * ablation


# ---------------------------------------------------------------------------
# published dynamics — vectorized port of `snet.timestep_diffdelay`
# ---------------------------------------------------------------------------
class ZebrafishSpine:
    """The published 8-population circuit, driven by per-side tonic drive.

    State: ``rates`` (N, 2, T) current population firing rates.  ``step``
    advances the published forward-Euler integrator by one published
    timestep (0.1 ms) with per-side drive scaling (documented DRIVE coupling
    above).  ``muscles()`` returns (2, N_BODY) filtered activations.
    """

    def __init__(self, params: dict | None = None, seed: int = 27):
        p = dict(PARAMS_8POP)
        if params:
            for key in ("fastdrivemod", "slowdrivemod"):
                if key in params:
                    p[key] = params[key]
        self.params = p
        self.N = p["Numsegs"]
        self.T = p["Numtypes"]
        self.names = p["typenames"]

        rng = np.random.default_rng(seed)
        self.rates = 0.1 * rng.random((self.N, 2, self.T))

        # connectivity matrix W[n1,a1,p1, n2,a2,p2] (published buildmatrix)
        N, T = self.N, self.T
        W = np.zeros((N, 2, T, N, 2, T))
        for n1 in range(N):
            for a1 in range(2):
                for p1 in range(T):
                    for n2 in range(N):
                        for a2 in range(2):
                            for p2 in range(T):
                                W[n1, a1, p1, n2, a2, p2] = connection(
                                    n1, a1, p1, n2, a2, p2, p)
        self.W = W

        # membrane time constants (published fast/slow assignment)
        self.mtcs = 1.0 * np.ones((N, 2, T))
        for t in range(T):
            name = p["typenames"][t]
            if 'fast' in name:
                self.mtcs[:, :, t] = p["fastmtc"]
            elif 'slow' in name:
                self.mtcs[:, :, t] = p["slowmtc"]

        # tonic drive (published: drives BOTH E and I populations)
        self.drives = 1.0 * np.ones((N, 2, T))
        for t in range(T):
            name = p["typenames"][t]
            if 'fast' in name:
                self.drives[:, :, t] = p["fastdrivemod"]
            elif 'slow' in name:
                self.drives[:, :, t] = p["slowdrivemod"]
        self.drive_scale = np.ones(2)   # per-side connectome gate (DRIVE)

        # speed masks over populations (published fast_mask / slow_mask)
        fast_mask = np.zeros((N, 2, T))
        slow_mask = np.zeros((N, 2, T))
        for t, name in enumerate(self.names):
            if "slow" in name:
                slow_mask[:, :, t] = 1
            elif 'fast' in name:
                fast_mask[:, :, t] = 1
            else:
                raise KeyError(f"population '{name}' lacks fast/slow speed")
        self.fast_mask = fast_mask
        self.slow_mask = slow_mask

        # weight blocks folded with the speed masks, batched per target so a
        # substep is two batched matvecs instead of the published 900
        # (identical arithmetic, different summation order)
        Wf = W * fast_mask[:, :, :, None, None, None]
        Ws = W * slow_mask[:, :, :, None, None, None]
        self.B_fast = np.empty((N, 2 * T, N * 2 * T))   # (target-i, tgt_pops, src)
        self.B_slow = np.empty((N, 2 * T, N * 2 * T))
        for i in range(N):
            self.B_fast[i] = Wf[:, :, :, i, :, :].reshape(N * 2 * T, 2 * T).T
            self.B_slow[i] = Ws[:, :, :, i, :, :].reshape(N * 2 * T, 2 * T).T

        # delay tables: delay[d, i2, i] = (1 + |i2 - i|) * delayts (steps)
        ii = np.arange(N)
        dist = np.abs(ii[:, None] - ii[None, :])
        self.d_fast = (1 + dist) * p["fastdelayts"]
        self.d_slow = (1 + dist) * p["slowdelayts"]
        self.i2grid = np.broadcast_to(ii[:, None], (N, N))
        self.ring_len = self.d_slow.max() + 2        # published max delay+2
        self.buf = np.zeros((self.ring_len, N, 2, T))
        self.buf[0] = self.rates
        self.prev = self.rates
        self.t = 0                                    # substep counter

        # per-target batched weights / mtc / drive vectors
        self._mtc_f = self.mtcs.reshape(N, 2 * T)
        self._drv_f = self.drives.reshape(N, 2 * T)

        # muscle readout state (published MUSCLE coupling)
        self.E_idx = [t for t, n in enumerate(self.names) if 'E' in n]
        self._mact = np.zeros((2, N_BODY))
        # rostro-caudal interpolation spine(30) -> body(25): segment k of the
        # body sits at fraction k/(25-1) of the spine's 0..29 span
        self._map_x = np.arange(N, dtype=float)
        self._map_q = np.linspace(0, N - 1, N_BODY)

    # ------------------------------------------------------------- dynamics
    def substep(self, drive_scale: np.ndarray | None = None):
        """One published forward-Euler step of 0.1 ms (vectorized)."""
        t = self.t
        N, T = self.N, self.T
        C = self.params["Cstrength"]

        # published sourcerate = fast-delayed fast-pool rates + slow-delayed
        # slow-pool rates; head holds rates[t] and we compute rates[t+1], so
        # the published read fullrates[t+1 - delay] is ring slot (t+1-d)
        rf = self.buf[(t + 1 - self.d_fast) % self.ring_len, self.i2grid]
        rs = self.buf[(t + 1 - self.d_slow) % self.ring_len, self.i2grid]
        # published startup quirk (verbatim port): the reference's `else`
        # branch resets the whole delayed source to zero until each pair's
        # slow delay has elapsed — reproduce that per-pair gate
        rf = rf * ((t + 1 - self.d_slow) >= 0)[:, :, None, None]
        # (i2,i,a2,p2) -> per-target batched source vector over (i2,a2,p2)
        s_fast = rf.transpose(1, 0, 2, 3).reshape(N, N * 2 * T)
        s_slow = rs.transpose(1, 0, 2, 3).reshape(N, N * 2 * T)

        contrib = np.matmul(self.B_fast, s_fast[..., None])[:, :, 0] \
            + np.matmul(self.B_slow, s_slow[..., None])[:, :, 0]
        contrib *= C

        ds = self.drive_scale if drive_scale is not None else self.drive_scale
        ds_pop = np.repeat(ds, T)                     # per (side, pop) flat
        mtc = self._mtc_f
        new = self._drv_f * ds_pop[None, :] / mtc + contrib / mtc
        new = np.maximum(new, 0.0)                    # published actfunc
        new -= self.prev.reshape(N, 2 * T) / mtc      # published decay term
        new = self.prev.reshape(N, 2 * T) + new * TSTEP
        new[new < 0] = 0.0

        self.t += 1
        self.prev = new.reshape(N, 2, T)
        self.buf[self.t % self.ring_len] = self.prev
        self.rates = self.prev

    # ---------------------------------------------------------- muscle readout
    def muscles(self, dt: float) -> np.ndarray:
        """(2, N_BODY) filtered muscle activations from the E (V2a) pools.

        Published MUSCLE coupling: ipsilateral E population rate → muscle
        output (Roussel et al. 2021), resampled from the 30 biological
        somites of the circuit to the 25 mechanical hemisegments.
        """
        r = self.rates[:, :, self.E_idx].mean(axis=2)      # (N, 2)
        rT = r.T                                           # (2, N)
        rT_body = [np.interp(self._map_q, self._map_x, rT[a])
                   for a in (0, 1)]
        target = np.stack(rT_body) / MUSCLE_RATE_SCALE
        target = target / (1.0 + np.maximum(target, 0.0))  # published saturate
        self._mact += (target - self._mact) * dt / MUSCLE_TAU
        return self._mact
