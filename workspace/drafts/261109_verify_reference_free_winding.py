# 261109_verify_reference_free_winding.py
"""
Numerical verification of the reference-free winding-number statements
derived in 261109_winding_angle_unified.md, starting from the two-band
definitions of notes/Note.tex (Sec. "Winding number").

Models are chiral two-band systems d(k) = r(k) (cos phi_k e1 + sin phi_k e2)
in the plane perpendicular to a RANDOM chiral axis n_Gamma, with a random
positive radius r(k) (the radius must drop out of every statement).

Checks (all must pass):
  c1  chiral-axis identity:  n_Gamma . (dhat x d/dk dhat) = d(phi)/dk
      (random axis, random radius -> identity is radius-independent)
  c2  W from the Note.tex integral = lift difference = degree of the
      relative-angle map k -> phi_k - phi_{k0} for several k0
      (reference independence)
  c3  W != 0  =>  for every k0 the relative angle covers [0, 2pi):
      max circular gap of {Delta phi mod 2pi} ~ 0, and the pi point
      (Delta phi = pi, i.e. dhat(k) = -dhat(k0)) is hit
  c4  retrace model: W = 0 while the image is still the full circle
      (coverage and pi point hold for every k0) -> coverage is necessary,
      not sufficient
  c5  ARPES shadow: for a random external axis z_hat with in-plane
      projection rho, dhat(k).z_hat = rho cos(phi_k - phi_z); the pi point
      dhat = -z_hat is unreachable unless rho = 1 (in-plane axis)
"""

import numpy as np

rng = np.random.default_rng(20261109)

N = 400001
ks = np.linspace(-np.pi, np.pi, N)


def trapz(y, x):
    return np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)


def random_unit(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def frame(n):
    """Right-handed orthonormal frame (e1, e2, n): e1 x e2 = n."""
    a = np.array([1.0, 0.0, 0.0])
    if abs(n @ a) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    e1 = a - (a @ n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return e1, e2


# ---------------- lifts phi(k) for the three test trajectories ----------------
def phi_T(k):  # monotone, W = 1
    return k + 0.5 * np.sin(k)


def phi_A(k):  # monotone, W = 2
    return 2.0 * k


def phi_R(k):  # retrace (go around once, come back), W = 0, image = full circle
    return np.where(k <= 0.0, 2.0 * k + np.pi, np.pi - 2.0 * k)


MODELS = [
    ("T  (monotone, W=1)", phi_T, 1.0),
    ("A  (monotone, W=2)", phi_A, 2.0),
    ("R  (retrace,   W=0)", phi_R, 0.0),
]

K0_LIST = [-2.7, -1.3, -0.4, 0.9, 2.2]
TOL = 1e-6

results = {}
for name, phi_fn, W_true in MODELS:
    n = random_unit(rng)
    e1, e2 = frame(n)
    # random positive radius: must not affect any statement below
    r = 1.0 + 0.3 * np.cos(2.0 * ks)

    ph = phi_fn(ks)
    d = r[:, None] * (np.cos(ph)[:, None] * e1[None, :]
                      + np.sin(ph)[:, None] * e2[None, :])
    dh = d / np.linalg.norm(d, axis=1, keepdims=True)

    # --- c1: n . (dhat x d/dk dhat) = d(phi)/dk (radius-independent) ---
    ddh = np.gradient(dh, ks, axis=0)
    integrand = np.cross(dh, ddh) @ n
    dph_num = np.gradient(np.unwrap(np.arctan2(dh @ e2, dh @ e1)), ks)
    c1 = np.max(np.abs(integrand - dph_num))

    # --- W from the Note.tex integral, and from the lift ---
    W_int = trapz(integrand, ks) / (2 * np.pi)
    ph_lift = np.unwrap(np.arctan2(dh @ e2, dh @ e1))
    W_lift = (ph_lift[-1] - ph_lift[0]) / (2 * np.pi)

    # --- c2/c3/c4: relative angle w.r.t. several trajectory references ---
    deg_dev = 0.0
    max_gap = 0.0
    pi_miss = 0.0
    for k0 in K0_LIST:
        i0 = np.argmin(np.abs(ks - k0))
        dph = ph_lift - ph_lift[i0]
        deg_dev = max(deg_dev, abs((dph[-1] - dph[0]) / (2 * np.pi) - W_int))
        v = np.sort(np.mod(dph, 2 * np.pi))
        gaps = np.diff(np.concatenate([v, [v[0] + 2 * np.pi]]))
        max_gap = max(max_gap, gaps.max())
        # pi point: cos(Delta phi) = -1 somewhere  <=>  dhat(k) = -dhat(k0)
        pi_miss = max(pi_miss, np.min(np.abs(np.cos(dph) + 1.0)))

    # --- c6: arbitrary FIXED in-plane reference (not on the trajectory) ---
    # any constant in-plane angle phi_a must give the same degree, the same
    # coverage, and the same pi point as a trajectory reference
    deg_dev_a = 0.0
    max_gap_a = 0.0
    pi_miss_a = 0.0
    for _ in range(5):
        phi_a = rng.uniform(-np.pi, np.pi)
        dph = ph_lift - phi_a
        deg_dev_a = max(deg_dev_a,
                        abs((dph[-1] - dph[0]) / (2 * np.pi) - W_int))
        v = np.sort(np.mod(dph, 2 * np.pi))
        gaps = np.diff(np.concatenate([v, [v[0] + 2 * np.pi]]))
        max_gap_a = max(max_gap_a, gaps.max())
        pi_miss_a = max(pi_miss_a, np.min(np.abs(np.cos(dph) + 1.0)))

    results[name] = dict(c1=c1, W_int=W_int, W_lift=W_lift, W_true=W_true,
                         deg_dev=deg_dev, max_gap=max_gap, pi_miss=pi_miss,
                         deg_dev_a=deg_dev_a, max_gap_a=max_gap_a,
                         pi_miss_a=pi_miss_a)

# ---------------- c5: ARPES shadow for an external fixed axis ----------------
name, phi_fn, W_true = MODELS[0]  # use T (W = 1)
n = random_unit(rng)
e1, e2 = frame(n)
ph = phi_fn(ks)
r = 1.0 + 0.3 * np.cos(2.0 * ks)
d = r[:, None] * (np.cos(ph)[:, None] * e1[None, :]
                  + np.sin(ph)[:, None] * e2[None, :])
dh = d / np.linalg.norm(d, axis=1, keepdims=True)

arp = {}
for label, z in [("random z (rho<1)", random_unit(rng)),
                 ("in-plane z  (rho=1)", np.cos(1.1) * e1 + np.sin(1.1) * e2)]:
    rho = np.linalg.norm(z - (z @ n) * n)
    zp = (z - (z @ n) * n)
    if rho > 1e-12:
        phi_z = np.arctan2(zp @ e2, zp @ e1)
        shadow_dev = np.max(np.abs(dh @ z - rho * np.cos(ph - phi_z)))
    else:
        shadow_dev = np.max(np.abs(dh @ z))
    # pi point for ARPES would need dhat(k) = -z, i.e. min_k dhat.z = -1
    min_dot = np.min(dh @ z)
    arp[label] = dict(rho=rho, shadow_dev=shadow_dev, min_dot=min_dot)

# ---------------- report ----------------
print("model                 |  c1 integrand  |  W(int)  W(lift)  W(true)"
      "  |  c2 deg dev  |  c3/c4 max gap  |  pi miss"
      "  |  c6 deg dev  |  c6 max gap  |  c6 pi miss")
ok = True
for name, res in results.items():
    line_ok = (res["c1"] < TOL
               and abs(res["W_int"] - res["W_true"]) < TOL
               and abs(res["W_lift"] - res["W_true"]) < TOL
               and res["deg_dev"] < TOL
               and res["deg_dev_a"] < TOL)
    print(f"{name:22s}|  {res['c1']:.2e}   |  {res['W_int']:+.6f}"
          f"  {res['W_lift']:+.6f}  {res['W_true']:+.1f}"
          f"  |  {res['deg_dev']:.2e}   |  {res['max_gap']:.2e}"
          f"     |  {res['pi_miss']:.2e}"
          f"  |  {res['deg_dev_a']:.2e}   |  {res['max_gap_a']:.2e}"
          f"     |  {res['pi_miss_a']:.2e}")
    ok &= line_ok

print("\n[c5] ARPES shadow, model T (W=1):")
for label, res in arp.items():
    print(f"  {label:20s}  rho = {res['rho']:.4f}   "
          f"|d.z - rho cos(phi-phi_z)| = {res['shadow_dev']:.2e}   "
          f"min_k d.z = {res['min_dot']:.4f}   (pi point needs -1)")
    ok &= res["shadow_dev"] < TOL
ok &= abs(arp["random z (rho<1)"]["min_dot"]
          + arp["random z (rho<1)"]["rho"]) < 1e-3
ok &= abs(arp["in-plane z  (rho=1)"]["min_dot"] + 1.0) < 1e-3

# c3/c4 semantic checks
ok &= results["T  (monotone, W=1)"]["max_gap"] < 1e-3
ok &= results["A  (monotone, W=2)"]["max_gap"] < 1e-3
ok &= results["R  (retrace,   W=0)"]["max_gap"] < 1e-3   # full circle, W=0!
ok &= results["R  (retrace,   W=0)"]["pi_miss"] < 1e-3   # pi point hit, W=0!

# c6 semantic checks: fixed in-plane (non-trajectory) references behave
# exactly like trajectory references
ok &= results["T  (monotone, W=1)"]["max_gap_a"] < 1e-3
ok &= results["A  (monotone, W=2)"]["max_gap_a"] < 1e-3
ok &= results["R  (retrace,   W=0)"]["max_gap_a"] < 1e-3
ok &= results["T  (monotone, W=1)"]["pi_miss_a"] < 1e-3
ok &= results["R  (retrace,   W=0)"]["pi_miss_a"] < 1e-3

print(f"\nAll checks passed: {ok}")
