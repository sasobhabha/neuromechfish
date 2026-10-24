"""Port validation: neuromechfish.spine vs the published reference implementation.

Compares the vectorized port (ZebrafishSpine) against the verbatim published
code (fwandler/ZFswim, `snet.py`/`parameters.py`), step by step, and then
computes the published analysis metrics (Wandler et al. 2025
`analysis_functions.py`) on both.

Run:
  .venv/bin/python scripts/validate_spine.py          # fresh reference run
  .venv/bin/python scripts/validate_spine.py --reuse  # reuse saved reference
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "_work" / "ZFswim"))

import numpy as np
from scipy.stats import circmean

from analysis_functions import get_amp, get_freq, get_phase  # noqa: E402  (published)

REF_PATH = _ROOT / "_work" / "ref_fullrates_8pop.npy"


def run_reference() -> np.ndarray:
    from parameters import load_default_params  # published
    from snet import SpinalNetwork              # published

    params = load_default_params("8pop")
    Spine = SpinalNetwork(params)
    Spine.params["fastdrivemod"] = 2.0
    Spine.params["slowdrivemod"] = 0.5
    Spine.update_model()
    t0 = time.time()
    Spine.simulation()
    print(f"published reference: {time.time() - t0:.1f} s")
    full = np.asarray(Spine.fullrates, dtype=np.float64)
    np.save(REF_PATH, full.astype(np.float32))
    return full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reuse", action="store_true", help="reuse cached reference run")
    args = ap.parse_args()

    if args.reuse and REF_PATH.exists():
        ref = np.load(REF_PATH).astype(np.float64)
        print("reusing cached published reference run")
    else:
        ref = run_reference()

    from neuromechfish.spine import ZebrafishSpine

    sp = ZebrafishSpine()
    t0 = time.time()
    rates = np.empty_like(ref)
    first_diff_step = None
    for t in range(ref.shape[0]):
        rates[t] = sp.rates
        d = np.max(np.abs(rates[t] - ref[t]))
        if d > 1e-3 and first_diff_step is None:
            first_diff_step = t
        sp.substep()
    print(f"port: {time.time() - t0:.1f} s")
    print("first step with |port - reference| > 1e-3:", first_diff_step)
    print("max |port - published reference| =", float(np.max(np.abs(rates - ref))))

    # published analysis metric set (mirrors analysis_functions.py main)
    def metrics(full):
        freqs = np.zeros(full.shape[1:])
        amps = np.zeros(full.shape[1:])
        phases = np.zeros(full.shape[1:])
        for seg in range(full.shape[1]):
            for side in range(2):
                for pop in range(full.shape[3]):
                    r = full[:, seg, side, pop]
                    f = get_freq(r)
                    freqs[seg, side, pop] = f
                    amps[seg, side, pop] = get_amp(r)
                    phases[seg, side, pop] = get_phase(r, f)
        f_avg = np.mean(np.sum(freqs * amps, axis=(1, 2)) / np.sum(amps, axis=(1, 2)))
        ph_seg = np.sum(np.angle(np.exp(-1j * (phases[1:] - phases[:-1])))
                        * amps[1:], axis=2) / np.sum(amps[1:], axis=2)
        ph_lr = np.sum(np.angle(np.exp(-1j * (phases[:, 0, :] - phases[:, 1, :])))
                       * amps[:, 0, :], axis=1) / np.sum(amps[:, 0, :], axis=1)
        return (f_avg, float(np.mean(amps)),
                float(np.degrees(circmean(ph_lr))),
                float(np.degrees(circmean(ph_seg))))

    for name, full in (("published reference", ref), ("port", rates)):
        f, a, lr, seg = metrics(full)
        print(f"{name}: freq {f:.1f} Hz, amp {a:.3f}, "
              f"LR phase {lr:+.0f} deg, interseg {seg:+.1f} deg/seg")

    e_pops = [t for t, n in enumerate(sp.names) if "E" in n]
    amp_E = np.mean([get_amp(ref[:, s, a, p])
                     for s in range(30) for a in range(2) for p in e_pops])
    print(f"E-population mean peak-to-peak amplitude: {amp_E:.3f}")


if __name__ == "__main__":
    main()
