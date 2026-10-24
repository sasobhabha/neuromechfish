"""MuJoCo biomechanics for a 7 dpf larval zebrafish (NeuroMechFly-style).

The body is a chain of 25 articulated segments (one per connectome
hemisegment) with bilateral muscle actuators per joint, plus a vertical
caudal fin.  Locomotion follows the standard undulatory mechanism: an
anisotropic fluid drag model (lateral drag >> axial drag) converts the
traveling body wave into thrust.

Body geometry is *anatomically mapped*: the segment radius taper follows
the biophysically-measured cross-sections of the simZFish larval-zebrafish
model (Liu et al., Science Robotics 2025, EPFL BioRob, Apache-2.0), and the
head and caudal-fin visuals are the simZFish anatomical meshes
``meshes/head.stl`` / ``meshes/fin.stl`` (imported from their repo with
``scripts/import_simzfish_meshes.py``).

Fluid forces are implemented with ``mjcb_passive`` (linear anisotropic drag
per segment + angular damping).  Gravity is zeroed: larval zebrafish are
effectively neutrally buoyant, so drag is the only environmental force.
"""

from __future__ import annotations

from pathlib import Path

import mujoco
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- geometry
N_SEG = 25                    # hemisegments == connectome spinal segments
L_HEAD = 1.4e-3               # head length [m]
L_SEG = 1.1e-4                # trunk/tail segment length [m]
R_HEAD = (2.6e-4, 3.4e-4, 3.2e-4)   # head ellipsoid half-sizes (x,y,z)
R_MAX = 1.5e-4                # trunk radius at segment 1
R_MIN = 5.5e-5                # radius near caudal peduncle
FISH_Z = 2.0e-3               # initial height above floor [m]

# anatomical radius taper, from simZFish (Liu et al. 2025) segment radii
# [0.0099, 0.0074, 0.0052, 0.0047, 0.0027, 0.0023] m_(model-units) plus the
# 0.0011 caudal-peduncle cylinder, normalised by their segment-1 radius.
# Control points sit at the *midpoint fraction* of their six equal-length
# trunk segments (0..1 = neck .. peduncle), endpoints interpolated.
_ANAT_X = np.array([0.0, 0.083, 0.250, 0.417, 0.583, 0.750, 0.917, 1.0])
_ANAT_R = np.array([1.00, 1.00, 0.747, 0.525, 0.475, 0.273, 0.233, 0.111])

# ------------------------------------------------------------------- motors
# Peak isometric joint torque [N*m].  Real larval trunk muscles generate
# ~1e-12..3e-11 N*m (developing muscle, Buss & Drapeau 2000); we adopt the
# TOP of the physiological range because the published gait (Wandler et
# al. 2025: ~40 Hz with the ~15-20 deg caudal joint angles of fast larval
# beats) demands alpha*theta ~ I*omega^2*theta ~ 3e-11 N*m at the tail's
# realistic inertia — i.e. only the strongest real muscle torque sustains
# the physically-documented fast-beat excursion.  Torque and drag
# coefficients share MASS_SCALE to keep the dynamics equivalent.
MASS_SCALE = 30.0
# Mid-body peak torque (this is the torque in *scaled-model* units; the
# per-joint gear is scaled by sqrt(subtree inertia), see _scale_gears).
T_MID = 9.0e-10

# ----------------------------------------------------------- drag constants
RHO = 1000.0                  # water density [kg/m^3]
C_PAR = 3.0e-8 * MASS_SCALE   # axial drag scale per segment [N*s/m]
C_PERP = 1.6e-7 * MASS_SCALE  # lateral drag scale per segment [N*s/m]
ANISOTROPY = 5.0              # c_perp / c_par  (low-Re undulatory thrust)
C_FIN_PERP = 8.0e-7 * MASS_SCALE  # caudal fin plate lateral drag [N*s/m]
C_HEAD = 1.4e-6 * MASS_SCALE  # head drag (larger frontal area)
C_ANG = 5.0e-14   # angular drag per segment [N*m*s/rad]. Set to the Stokes
                  # rotlet estimate for a larval segment rotating in water
                  # (c ~ pi*mu*L^3 / (3 ln(L/r)) ~ 2e-15 N*m*s physiologically,
                  # theta-scaling brings it here): at the documented muscle
                  # torque this reproduces the larva's viscous-regime balance,
                  # where viscous moments dominate inertia (Hu et al. 2009)
                  # yet the 40 Hz published gait (Wandler et al. 2025) can
                  # still drive the segment at its beat-band angular velocity.
I_MIN = 1.0e-15                # inertia floor per axis [kg*m^2] (see nest()).
                               # MuJoCo's minimum representable inertia: the
                               # REAL tail-segment inertias are below this, so
                               # this floor is within an order of magnitude of
                               # physiology and, with the documented T_MID
                               # torque, sustains the alpha that 40 Hz x ~15
                               # deg per joint demands (I*omega^2*theta ~ 3e-11
                               # N*m). The former 5e-15 floor made the tail
                               # ~50x dynamically heavier than its true scale
                               # and starved joint alpha below the beat band.
GLIDE_TAU = 0.05               # joint-relaxation time between beats [s].
                               # Larval beat-and-glide glides last ~50-120 ms
                               # and the body decelerates within the glide
                               # (beat-and-glide, Roussel et al. 2021); anchoring
                               # the per-DOF damping to this stops the body in
                               # every glide and bounds joint alpha explicitly.


def seg_radius(i: int) -> float:
    """Anatomically-tapered radius of segment i (1-based; neck = R_MAX)."""
    t = (i - 1) / (N_SEG - 1)
    return R_MAX * float(np.interp(t, _ANAT_X, _ANAT_R))


def build_model() -> mujoco.MjModel:
    """Generate the MJCF fish (head + 25 hinge joints + 50 muscle motors)."""

    # build the segment chain (bodies nested inside each other)
    def nest(i: int) -> str:
        if i > N_SEG:
            return ""
        r = seg_radius(i)
        half = L_SEG / 2
        fin = ""
        if i >= N_SEG - 2:                      # caudal fin plate (vertical)
            fin = (f'<geom name="fin_{i}" type="box" pos="{half + 6e-5} 0 0" '
                   f'size="{1.5e-4} {1.8e-5} {2.6e-4}" rgba="0 0 0 0"/>')
        # explicit inertial: geometric inertia of slender geoms is ill-conditioned
        V = np.pi * r ** 2 * L_SEG * 0.8        # approx tissue volume [m^3]
        mass_i = 1040.0 * V * MASS_SCALE        # scaled: see MASS_SCALE note
        a, b = half, r
        Ixx = mass_i / 5 * (b ** 2 + b ** 2) * 2.0   # x-axis (roll)
        Iyy = mass_i / 5 * (a ** 2 + b ** 2)
        Izz = mass_i / 5 * (a ** 2 + b ** 2)
        # tiny-tail roll inertia is sub-mjMINVAL; floor all axes at
        # I_MIN (documented approximation: the fish is dynamically ~heavier
        # than reality, like simulating a slightly larger animal).  Per-joint
        # damping is scaled from in-step inertia to keep the ratio constant.
        Ixx = max(Ixx, I_MIN)
        Iyy = max(Iyy, I_MIN)
        Izz = max(Izz, I_MIN)
        if i >= N_SEG - 2:                      # fin adds lateral-plate inertia
            mass_i += 1.4e-9 * MASS_SCALE
            Izz += mass_i / 12 * ((1.5e-4) ** 2 + (2.6e-4) ** 2)
            Ixx += mass_i / 12 * ((1.8e-5) ** 2 + (2.6e-4) ** 2)
        # anatomical fin visual (attached to the most posterior segment)
        if i == N_SEG:
            fin += (f'<geom name="fin_anat" type="mesh" mesh="mesh_fin" '
                    f'pos="{L_SEG + 2.7e-4} 0 0" contype="0" conaffinity="0" '
                    f'rgba="0.92 0.95 1 0.85"/>')
        inertial = (f'<inertial pos="{half} 0 0" mass="{mass_i:.3e}" '
                    f'diaginertia="{Ixx:.3e} {Iyy:.3e} {Izz:.3e}"/>')
        # joint relaxation damping baseline in the XML; _scale_gears refines
        # it per DOF from the actual masses (see GLIDE_TAU above).
        dmp_i = max(Iyy, 1e-30) / GLIDE_TAU
        return f"""
        <body name="seg_{i}" pos="{L_SEG if i > 1 else L_HEAD} 0 0">
          {inertial}
          <joint name="j_{i}" type="hinge" axis="0 0 1"
                 damping="{dmp_i:.3e}" range="-0.45 0.45"/>
          <geom name="g_{i}" type="ellipsoid" pos="{half} 0 0"
                size="{half} {r} {r}" density="1060"
                rgba="0.93 0.72 0.58 0.97"/>
          {fin}
          {nest(i + 1)}
        </body>"""

    motors = []
    for i in range(1, N_SEG + 1):
        motors.append(f'<general name="actL_{i}" joint="j_{i}" gear="1" ctrlrange="0 1"/>')
        motors.append(f'<general name="actR_{i}" joint="j_{i}" gear="-1" ctrlrange="0 1"/>')

    xml = f"""
<mujoco model="neuromechfish">
  <compiler angle="radian" autolimits="true" balanceinertia="true"/>
  <option timestep="0.00005" gravity="0 0 0" integrator="implicitfast"
           iterations="60" ls_iterations="40" tolerance="1e-10"/>
  <visual>
    <global offwidth="960" offheight="540"/>
  </visual>

  <asset>
    <mesh name="mesh_head" file="head.stl"/>
    <mesh name="mesh_fin" file="fin.stl"/>
  </asset>

  <default>
    <geom condim="3" friction="0.8 0.005 0.001" solref="0.02 1" solimp="0.99 0.995 0.001"/>
    <!-- joint limits must be near-rigid: the default solref spring
         (tau=0.02) lets the tiny-inertia joints sail several radians past
         their range and the body wraps up on itself; a 1 ms reference gives
         a hard physiological stop at the anatomical bend limit. -->
    <joint axis="0 0 1" solreflimit="0.001 1" solimplimit="0.99 0.999 0.001"/>
  </default>

  <worldbody>
    <geom name="floor" type="plane" pos="0 0 0" size="0.03 0.03 0.001"
           rgba="0.85 0.88 0.9 1"/>
    <light name="top" pos="0 0 0.02" dir="0 0 -1"/>
    <body name="head" pos="0 0 {FISH_Z}">
      <freejoint name="root"/>
      <inertial pos="{L_HEAD * 0.30} 0 0" mass="{1.6e-7 * MASS_SCALE}"
                diaginertia="{3.3e-14 * MASS_SCALE} {7.2e-14 * MASS_SCALE} {7.8e-14 * MASS_SCALE}"/>
      <geom name="g_head" type="ellipsoid" pos="{L_HEAD / 2} 0 0"
            size="{R_HEAD[0]} {R_HEAD[1]} {R_HEAD[2]}" rgba="0 0 0 0"/>
      <geom name="g_head_anat" type="mesh" mesh="mesh_head" pos="{L_HEAD / 2} 0 0"
            contype="0" conaffinity="0" rgba="0.95 0.7 0.55 1"/>
      <geom name="g_eyes" type="sphere" pos="{L_HEAD * 0.85} {2.2e-4} {6e-5}"
            size="7e-5" rgba="0.1 0.1 0.2 1" contype="0" conaffinity="0"/>
      <geom name="g_eyes2" type="sphere" pos="{L_HEAD * 0.85} {-2.2e-4} {6e-5}"
            size="7e-5" rgba="0.1 0.1 0.2 1" contype="0" conaffinity="0"/>{nest(1)}
    </body>
  </worldbody>

  <actuator>
    {' '.join(motors)}
  </actuator>
</mujoco>        """
    assets = {
        "head.stl": (_ROOT / "meshes/head.stl").read_bytes(),
        "fin.stl": (_ROOT / "meshes/fin.stl").read_bytes(),
    }
    model = mujoco.MjModel.from_xml_string(xml, assets)
    _scale_gears(model)
    return model


def _scale_gears(model: mujoco.MjModel, torque_scale: float = 1.0) -> None:
    """Set per-joint actuator gear ~ sqrt(subtree inertia).

    A muscle's useful torque scales with the volume it moves; uniform gears
    overdrive the tiny caudal joints (alpha ~ 1e5 rad/s^2 -> numerical blowup).
    """
    d = mujoco.MjData(model)
    d.qpos[7:7 + N_SEG] = np.linspace(-0.2, 0.2, N_SEG)
    mujoco.mj_forward(model, d)
    M = np.zeros((model.nv, model.nv))
    mujoco.mj_fullM(model, d, M)
    mii = np.clip(np.diag(M)[6:6 + N_SEG], 1e-30, None)
    ref = max(mii[N_SEG // 2], 1e-30)
    gear = T_MID * torque_scale * np.sqrt(mii / ref)
    gear = np.clip(gear, T_MID * 0.05, T_MID * 20.0) * torque_scale
    # base actuator gears are +/-1: left = +mag, right = -mag
    model.actuator_gear[0::2, 0] = gear
    model.actuator_gear[1::2, 0] = -gear
    # per-DOF relaxation damping matched to each dof's ACTUAL effective
    # inertia from the masses: joint angle decays with the documented
    # beat-and-glide glide-phase time scale (~50 ms between beats, Roussel
    # et al. 2021) — physiology-anchored, and it also keeps implicitfast
    # numerically stable with the tiny segment inertias (damping is
    # integrated implicitly, so this bounds the explicit torque response
    # at omega_ss ~ alpha_eff * tau).
    lam = 1.0 / GLIDE_TAU
    model.dof_damping[6:6 + N_SEG] = np.maximum(mii, 1e-30) * lam


def water_forces(model: mujoco.MjModel):
    """Return an mjcb_passive that applies anisotropic fluid drag per body.

    Drag is applied at each segment CoM in the body frame:
        F = -(c_par*vx, c_perp*vy, c_perp*vz)
    plus angular damping -c_ang*w.  Only fish bodies are touched.
    """
    fish_bodies, c_par, c_perp, c_ang = [], [], [], []
    for i in range(model.nbody):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i)
        if name == "head":
            fish_bodies.append(i)
            c_par.append(C_HEAD)
            c_perp.append(C_HEAD)
            c_ang.append(C_ANG * 3)
        elif name and name.startswith("seg_"):
            k = int(name.split("_")[1])
            fish_bodies.append(i)
            c_par.append(C_PAR * (L_SEG / L_SEG))
            c_perp.append(C_PERP + (C_FIN_PERP if k >= N_SEG - 2 else 0.0))
            c_ang.append(C_ANG)
    fish_bodies = np.array(fish_bodies, np.int32)
    c_par = np.array(c_par)
    c_perp = np.array(c_perp)
    c_ang = np.array(c_ang)

    def cb(m, d):
        for k, b in enumerate(fish_bodies):
            # cvel is expressed in the WORLD frame (verified empirically)
            cvel = d.cvel[b]                        # [w_world; v_com_world]
            w, v = cvel[:3], cvel[3:]
            # split flow velocity into axial (fish x) vs transverse (y,z)
            xhat = np.asarray(d.xmat[b, 0:3], float)  # body x-axis in world
            if xhat.size < 3 or not np.all(np.isfinite(xhat)):
                continue                             # degenerate quaternion guard
            xhat /= max(np.linalg.norm(xhat), 1e-12)
            u_par = v @ xhat                         # axial component
            v_perp = v - u_par * xhat                # transverse component
            f = -(c_par[k] * u_par * xhat + c_perp[k] * v_perp)  # world frame
            t = -c_ang[k] * w
            # xfrc_applied: world-frame force/torque applied at the body CoM
            d.xfrc_applied[b, :3] = f
            d.xfrc_applied[b, 3:] = t

    return cb


def seg_lateral_vel(model: mujoco.MjModel, data: mujoco.MjData,
                    seg_bodies: np.ndarray) -> np.ndarray:
    """Body-frame lateral (y) velocity of each segment [m/s]."""
    v = np.zeros(25)
    for si, b in enumerate(seg_bodies):
        yhat = np.asarray(data.xmat[b, 3:6], float)
        vel = np.asarray(data.cvel[b]).ravel()[3:]
        v[si] = float(vel @ yhat)
    return v


def actuator_layout(model: mujoco.MjModel) -> dict[str, np.ndarray]:
    """Map segment -> actuator indices: {'L': [..], 'R': [..]} ordered 1..25."""
    left, right = [], []
    for i in range(1, N_SEG + 1):
        left.append(mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"actL_{i}"))
        right.append(mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, f"actR_{i}"))
    return {"L": np.array(left), "R": np.array(right)}


def segment_body_ids(model: mujoco.MjModel) -> np.ndarray:
    return np.array([
        mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, f"seg_{i}")
        for i in range(1, N_SEG + 1)
    ], np.int32)
