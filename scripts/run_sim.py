#!/usr/bin/env python
"""NeuroMechFish: closed-loop connectome-driven swimming simulation.

Couples the 2,616-node zebrafish sensory-motor-pathway connectome
(rate network) to the MuJoCo larval-zebrafish body via bilateral muscle
activations, with lateral-line-like sensory feedback.

Usage:
  .venv/bin/python scripts/run_sim.py [--dur 6] [--cmd 1.1] [--startle 3]
  .venv/bin/python scripts/run_sim.py --viewer        # native MuJoCo viewer window
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

import matplotlib
matplotlib.use("Agg")                 # under mjpython, main() runs off-thread
import matplotlib.pyplot as plt
import mujoco
import mujoco.viewer
import numpy as np

from neuromechfish import body as fish_body
from neuromechfish import network as fish_net

OUT = _ROOT / "results"


def _ensure_mjpython():
    """Re-exec the script under ``mjpython`` (required by launch_passive)."""
    if sys.platform != "darwin" or os.environ.get("_NMF_IN_MJPYTHON") == "1":
        return
    mjpython = Path(sys.executable).parent / "mjpython"
    if mjpython.exists():
        os.environ["_NMF_IN_MJPYTHON"] = "1"   # survives os.execv
        os.execv(str(mjpython), [str(mjpython), *sys.argv])


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dur", type=float, default=6.0)
    p.add_argument("--cmd", type=float, default=None, help="override command gain")
    p.add_argument("--g-sens", type=float, default=None)
    p.add_argument("--startle", type=float, default=None, help="startle time [s]")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--viewer", action="store_true",
                   help="live NATIVE MuJoCo viewer window (draggable camera, "
                        "user UI; launched under mjpython automatically)")
    p.add_argument("--rec-limit", type=float, default=6.0,
                   help="seconds to record in viewer mode for metrics/plots")
    p.add_argument("--video", action="store_true",
                   help="headless render of the swim to results/fish_swim.mp4")
    return p.parse_args()


def main():
    args = parse_args()
    OUT.mkdir(exist_ok=True)

    print("building fish body model...", flush=True)
    model = fish_body.build_model()
    data = mujoco.MjData(model)
    cb = fish_body.water_forces(model)
    mujoco.set_mjcb_passive(cb)

    net = fish_net.SMPNetwork(_ROOT / "data/processed/connectome.npz", seed=args.seed)
    if args.cmd is not None:
        fish_net.G_CMD = args.cmd
    if args.g_sens is not None:
        fish_net.G_SENS = args.g_sens

    acts = fish_body.actuator_layout(model)          # {'L': ids, 'R': ids}
    seg_bodies = fish_body.segment_body_ids(model)   # body id per segment
    head_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "head")

    viewer = None
    if args.viewer:
        _ensure_mjpython()                          # relaunch under mjpython
        print("opening native MuJoCo viewer window...", flush=True)
        viewer = mujoco.viewer.launch_passive(model, data)
        # zoom onto the larva (the default frame is dominated by the floor)
        vcam = viewer.cam
        vcam.type = mujoco.mjtCamera.mjCAMERA_TRACKING
        vcam.trackbodyid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "seg_10")
        vcam.azimuth, vcam.elevation, vcam.distance = 90.0, -10.0, 0.0075
        print("viewer running:", viewer.is_running(), flush=True)

    renderer = ffmpeg_proc = cam = None
    if args.video:
        import subprocess
        W, H = 960, 540
        renderer = mujoco.Renderer(model, height=H, width=W)
        cam = mujoco.MjvCamera()
        cam.azimuth, cam.elevation = 90.0, -8.0     # side view
        cam.distance = 0.008                        # fish is ~1.5 mm long
        cmd = ["ffmpeg", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{W}x{H}", "-r", "25", "-i", "-",
               "-c:v", "libx264", "-pix_fmt", "yuv420p",
               str(OUT / "fish_swim.mp4")]
        ffmpeg_proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                       stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL)

    # fixed-length headless run vs. viewer run (window drives until rec_limit)
    eff_dur = args.dur if args.video else (args.rec_limit if viewer is not None else args.dur)
    n_steps = int(eff_dur / model.opt.timestep)

    print(f"simulating {eff_dur:.1f} s ({n_steps} steps)...", flush=True)
    # warm the controller for 1 s (no physics) so muscles start mid-rhythm
    I_cmd = net.command_input()
    for _ in range(1000):
        net.step(model.opt.timestep, I_cmd)
        net.readout(model.opt.timestep)
        I_cmd = net.command_input()      # CPG is time-dependent: refresh!

    # recording buffers
    n_rec_total = n_steps
    n_rec = n_rec_total
    rec_t, rec_x, rec_y = [], [], []
    rec_q = np.zeros((n_rec_total, 25))
    rec_act = np.zeros((n_rec_total, 2, 25))
    rec_speed = np.zeros(n_rec_total)

    for k in range(n_steps):
        t = k * model.opt.timestep

        # ---- sensory: lateral velocity of each segment (body frame) ----
        v_y = fish_body.seg_lateral_vel(model, data, seg_bodies)
        I = net.command_input() + net.sensory_input(v_y)
        if args.startle is not None:
            I = I + net.startle_input(t, args.startle)

        net.step(model.opt.timestep, I)
        act = net.readout(model.opt.timestep)            # (2, 25)

        ctrl = np.zeros(model.nu)
        ctrl[acts["L"]] = act[0]
        ctrl[acts["R"]] = act[1]
        data.ctrl[:] = ctrl
        mujoco.mj_step(model, data)

        if viewer is not None and k % 5 == 0:        # live native window
            if not viewer.is_running():              # user closed the window
                print("viewer closed — continuing headless...")
                viewer = None
            else:
                viewer.sync()

        if renderer is not None and k % 40 == 0:     # 25 fps at dt=1 ms
            tail_b = seg_bodies[-1]
            lookat = 0.5 * (np.array(data.xpos[head_id]) + np.array(data.xpos[tail_b]))
            cam.lookat[:] = lookat
            cam.lookat[2] = 2e-3
            renderer.update_scene(data, camera=cam)
            frame = renderer.render()
            ffmpeg_proc.stdin.write(frame.tobytes())

        rec_t.append(t)
        rec_x.append(data.xpos[head_id, 0]); rec_y.append(data.xpos[head_id, 1])
        rec_q[k] = data.qpos[7: 7 + 25]                  # after freejoint (7)
        rec_act[k] = act
        rec_speed[k] = np.linalg.norm(data.qvel[0:3]) if k > 0 else 0.0

    if renderer is not None:
        renderer.close()
        ffmpeg_proc.stdin.close()
        ffmpeg_proc.wait()
        print(f"wrote {OUT}/fish_swim.mp4")

    rec_t = np.array(rec_t); rec_x = np.array(rec_x); rec_y = np.array(rec_y)
    rec_speed = np.array(rec_speed)

    # ---------------------------------------------------------------- metrics
    width = 1.0            # skip warm-up window
    m = rec_t > width
    if m.sum() <= 2:
        print("run too short for metrics (closed window early?) — no plots written")
        return
    displacement = np.hypot(rec_x[m][-1] - rec_x[m][0], rec_y[m][-1] - rec_y[m][0])
    bl = fish_body.L_HEAD + fish_body.N_SEG * fish_body.L_SEG
    mean_speed = displacement / (rec_t[m][-1] - rec_t[m][0])
    tail_q = rec_q[m, -5]                                   # segment 20
    tail_q -= tail_q.mean()
    crossings = np.sum(np.abs(np.diff(np.sign(tail_q))) > 0)
    beat_freq = crossings / 2 / (rec_t[m][-1] - rec_t[m][0])
    tail_amp = np.ptp(tail_q)

    # wave direction: lag of distal vs proximal joint angle
    from numpy.fft import rfft, irfft
    a, b = rec_q[m, 4] - rec_q[m, 4].mean(), rec_q[m, 20] - rec_q[m, 20].mean()
    corr = irfft(rfft(a) * np.conj(rfft(b)), n=len(a))
    lag = int(np.argmax(corr))                              # samples a leads b
    if lag > len(a) / 2:
        lag -= len(a)
    lag_s = lag * model.opt.timestep

    # beat frequency via the published method (Wandler et al. 2025
    # analysis_functions.get_freq: autocorrelation peak); it hardcodes its
    # own sr=10000 (0.1 ms grid), so rescale for our 1 kHz kinematics sampling
    sys.path.insert(0, str(_ROOT / "_work" / "ZFswim"))
    from analysis_functions import get_freq as _pub_get_freq
    try:
        beat_freq_pub = _pub_get_freq(tail_q) \
            / (1.0 / (10000 * model.opt.timestep))
    except ValueError:
        beat_freq_pub = float("nan")

    print(f"\n===== NeuroMechFish results ({rec_t[-1]:.1f} s sim) =====")
    print(f"displacement      : {displacement*1e3:.2f} mm  ({displacement/bl:.2f} BL)")
    print(f"mean speed        : {mean_speed*1e3:.2f} mm/s ({mean_speed/bl:.2f} BL/s)")
    print(f"tail-beat freq    : {beat_freq:.2f} Hz (crossings) / "
          f"{beat_freq_pub:.2f} Hz (published autocorr)")
    print(f"tail joint amp    : {np.degrees(tail_amp):.1f} deg")
    print(f"wave lag j5->j21  : {lag_s*1e3:+.1f} ms ({'head->tail OK' if lag_s > 0 else 'REVERSED'})")
    # acceptance thresholds, at documented operating point: the published
    # circuit drives the larval gait at ~40 Hz (Wandler et al. 2025 run;
    # larval zebrafish beats 20-80 Hz, Budick & O'Malley 2000) and must
    # translate the body (glide between beats is the published phenotype of
    # beat-and-glide, Roussel et al. 2021)
    ok = (mean_speed > 0.05 * bl) and (15.0 < beat_freq_pub < 85.0)
    print("VALIDATION        :", "PASS" if ok else "FAIL (drive/sensory gain?)")

    # ------------------------------------------------------------------ plots
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    im0 = axs[0].pcolormesh(rec_t, np.arange(1, 26), rec_q.T, cmap="RdBu_r",
                            vmin=-0.35, vmax=0.35, shading="auto")
    axs[0].set_ylabel("segment"); axs[0].set_title("joint angles [rad]")
    plt.colorbar(im0, ax=axs[0])
    im1 = axs[1].pcolormesh(rec_t, np.arange(1, 26), rec_act[:, 0, :].T,
                            cmap="magma", vmin=0, vmax=1, shading="auto")
    axs[1].set_ylabel("segment"); axs[1].set_title("left muscle activation")
    plt.colorbar(im1, ax=axs[1])
    axs[2].plot(rec_t, rec_speed * 1e3, label="speed")
    axs[2].set_ylabel("speed [mm/s]"); axs[2].set_xlabel("t [s]")
    axs[2].legend(); axs[2].set_title("swim speed")
    fig.tight_layout()
    fig.savefig(OUT / "sim_kinematics.png", dpi=140)

    fig2, axs2 = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axs2[0].plot(rec_t[m], np.degrees(rec_q[m, -1]), label="seg 25")
    axs2[0].plot(rec_t[m], np.degrees(rec_q[m, -8]), label="seg 18")
    axs2[0].set_ylabel("angle [deg]"); axs2[0].legend()
    axs2[0].set_title("tail angle")
    axs2[1].plot(rec_t, rec_speed * 1e3)
    axs2[1].set_ylabel("speed [mm/s]"); axs2[1].set_xlabel("t [s]")
    fig2.tight_layout(); fig2.savefig(OUT / "sim_tail.png", dpi=140)
    print(f"wrote {OUT}/sim_kinematics.png, {OUT}/sim_tail.png")

    if viewer is not None:
        if os.environ.get("NMF_NOWAIT") == "1":      # for scripted runs/tests
            viewer.close()
            return
        print("sim done — close the viewer window to exit.", flush=True)
        hold = float(os.environ.get("NMF_AUTOCLOSE", "1e9"))
        t_end = time.time() + hold
        while viewer.is_running() and time.time() < t_end:
            viewer.sync()
            time.sleep(0.02)
        viewer.close()


if __name__ == "__main__":
    main()
