"""Numerical verification of the derivation chain: intensity cosine forms -> winding number.

Checks (all to machine precision unless noted):
 c1  SPSF/Bragg cosine formulas, incl. delta-string general formula
 c2  mirror invariance (d_y -> -d_y) of the full delta=0 data set
 c3  winding W flips sign under the mirror
 c4  exact identity  data(mirror S; delta) = data(S; -delta)
 c5  chiral full-transfer theorem: W != 0 -> W_max(k0) = 1 for all k0; antipode-hit count = |W|
 c6  DC double-cover: W = 0 but W_max = 1 for all k0 (false positive), hit count = 2
 c7  pole reference (Kitaev k0 = pi): Bragg weights pointwise equal SPSF weights
 c8  delta = 1/2 erasure (uniform winding pair W=1 vs W=0; real SSH pair max w both 1/2)
 c9  SPSF pi-crossing count = 2|W| (monotone), = 2 for the folded (W=0) false positive
 c10 folded-Kitaev false positive: W = 0, W_max < 1 strictly for all k0
 c11 minimal false-positive pair T (W=1) vs R (retrace, W=0): SPSF identical,
      pole-reference Bragg identical, generic-reference Bragg separates (1 vs 0.5)
 c13 delta-string rotation form w = (1 -/+ dhat(k).R_z(q*delta).dhat(k0))/2
      (incl. non-chiral); pole reference (k0=pi): Bragg(delta) = SPSF for all delta
"""
import numpy as np

# ---------------- core routines ----------------
def H_of(d):
    dx, dy, dz = d
    return np.array([[dz, dx - 1j*dy], [dx + 1j*dy, -dz]], dtype=complex)

def eig_states(d):
    E, V = np.linalg.eigh(H_of(d))          # E ascending: E[0]=-|d| (lower), E[1]=+|d|
    return V[:, 0], V[:, 1]                  # u_-, u_+

def spsf(d):
    u_m, u_p = eig_states(d)
    A = np.array([1.0, 0.0 + 0j])
    return abs(np.vdot(A, u_p))**2, abs(np.vdot(A, u_m))**2

def bragg(d0, dk):
    u_m0, _ = eig_states(d0)
    _, u_pk = eig_states(dk)
    u_mk, _ = eig_states(dk)
    return abs(np.vdot(u_pk, u_m0))**2, abs(np.vdot(u_mk, u_m0))**2

def bragg_delta(d0, dk, q, delta):
    u_m0, _ = eig_states(d0)
    _, u_pk = eig_states(dk)
    F = np.diag([1.0, np.exp(1j*q*delta)])
    return abs(np.vdot(u_pk, F @ u_m0))**2

def winding(z_loop):
    z_ext = np.concatenate([z_loop, z_loop[:1]])
    ang = np.unwrap(np.angle(z_ext))
    return int(np.round((ang[-1] - ang[0]) / (2*np.pi)))

def mirror(d_func):
    def mf(k):
        d = d_func(k)
        return np.array([d[0], -d[1], d[2]])
    return mf

# ---------------- models (all: meridian models wind about x via d_z + i d_y) ----------------
def m_random(k):  r = 1.2 + 0.5*np.cos(k);  return np.array([r*np.cos(k), r*np.sin(k), 0.7*np.cos(2*k)])  # non-chiral, W_z = 1
def m_A1(k):      return np.array([np.cos(k), np.sin(k), 0.0])                    # equator, W = 1
def m_A2(k):      return np.array([np.cos(2*k), np.sin(2*k), 0.0])                # equator, W = 2
def m_A0(k):      return np.array([1.0, 0.0, 0.0])                                # constant, W = 0
def m_B1(k):      return np.array([0.0, np.sin(k), np.cos(k)])                    # meridian T, W = 1
def m_R(k):       psi = k if k <= 0 else -k;  return np.array([0.0, np.sin(psi), np.cos(psi)])   # retrace, W = 0
def m_DC(k):      psi = 2*k + np.pi if k < 0 else np.pi - 2*k;  return np.array([0.0, np.sin(psi), np.cos(psi)])  # W = 0, image = S^1
def m_fold(k):    return np.array([0.0, np.cos(k), -0.5 - np.cos(k)])             # folded Kitaev, W = 0, SPSF false pos.
def m_sshT(k):    return np.array([1.0 + 2.0*np.cos(k), 2.0*np.sin(k), 0.0])       # SSH topological, W = 1
def m_ssh0(k):    return np.array([2.0 + np.cos(k), np.sin(k), 0.0])               # SSH trivial, W = 0

K = np.linspace(-np.pi, np.pi, 2001)[:-1] + np.pi/1000          # midpoints (avoid poles/endpoints)
K0S = np.linspace(-np.pi*0.9, np.pi*0.9, 37)
KEDGE = np.linspace(-np.pi, np.pi, 4001)                        # includes window edges (for c8)

def full_data(d_func, delta=0.0):
    sp = np.array([spsf(d_func(k)) for k in K])
    bg = np.empty((len(K0S), len(K)))
    for i, k0 in enumerate(K0S):
        d0 = d_func(k0)
        for j, k in enumerate(K):
            bg[i, j] = bragg_delta(d0, d_func(k), k - k0, delta)
    return sp, bg

def sph(d):
    dh = d/np.linalg.norm(d)
    return np.arccos(np.clip(dh[2], -1, 1)), np.angle(dh[0] + 1j*dh[1])

def wmax_hits(f, k0_list):
    out = []
    for k0 in k0_list:
        d0 = f(k0); dh0 = d0/np.linalg.norm(d0)
        w = np.array([bragg(d0, f(k))[0] for k in K])
        dh = np.array([f(k)/np.linalg.norm(f(k)) for k in K])
        near = (dh @ dh0) < -1 + 1e-4          # near-antipode
        clusters = int(near[0]) + int(np.sum(near[1:] & ~near[:-1]))
        if near[0] and near[-1]:                # merge wrap-around seam
            clusters -= 1
        out.append((w.max(), clusters))
    return out

def wind_of(f, about_x=False):
    zs = np.concatenate([K, K[:1]])
    if about_x:
        z = np.array([f(k)[2] + 1j*f(k)[1] for k in zs])
    else:
        z = np.array([f(k)[0] + 1j*f(k)[1] for k in zs])
    return winding(z)

R = []
# c1
for name, f in [("random", m_random), ("A1", m_A1), ("B1", m_B1), ("R", m_R), ("fold", m_fold), ("DC", m_DC)]:
    e = 0.0
    for k in K[::7]:
        d = f(k); dh = d/np.linalg.norm(d)
        w_p, w_m = spsf(d)
        e = max(e, abs(w_p - 0.5*(1+dh[2])), abs(w_m - 0.5*(1-dh[2])))
    for k0 in K0S[::4]:
        d0 = f(k0)
        for k in K[::13]:
            d = f(k); dh0 = d0/np.linalg.norm(d0); dh = d/np.linalg.norm(d)
            w_i, w_t = bragg(d0, d)
            e = max(e, abs(w_i - 0.5*(1 - dh@dh0)), abs(w_t - 0.5*(1 + dh@dh0)))
            q = k - k0; beta = q*0.37
            th0, ph0 = sph(d0); th, ph = sph(d)
            wds = bragg_delta(d0, d, q, 0.37)
            e = max(e, abs(wds - 0.5*(1 - np.cos(th0)*np.cos(th) - np.sin(th0)*np.sin(th)*np.cos((ph-ph0) - beta))))
    R.append((f"c1 formulas {name}", f"{e:.2e}"))

# c2 + c3
for name, f, ax in [("random W_z", m_random, False), ("A1", m_A1, False), ("A2", m_A2, False),
                    ("B1", m_B1, True), ("R", m_R, True), ("fold", m_fold, True), ("DC", m_DC, True)]:
    sp, bg = full_data(f, 0.0)
    spm, bgm = full_data(mirror(f), 0.0)
    dmax = max(np.abs(sp - spm).max(), np.abs(bg - bgm).max())
    R.append((f"c2/c3 mirror {name}: data diff, W(S), W(MS)",
              f"{dmax:.2e}, {wind_of(f, ax)}, {wind_of(mirror(f), ax)}"))

# c4
for name, f in [("A1", m_A1), ("B1", m_B1), ("random", m_random)]:
    dmax = 0.0
    for delta in (0.3, 0.7):
        _, bg = full_data(f, delta)
        _, bgm = full_data(mirror(f), delta)
        _, bgn = full_data(f, -delta)
        dmax = max(dmax, np.abs(bgm - bgn).max())
    R.append((f"c4 delta-mirror identity {name}", f"{dmax:.2e}"))

# c5 / c6 / c10
for name, f in [("A1 (W=1)", m_A1), ("A2 (W=2)", m_A2), ("B1 (W=1)", m_B1),
                ("DC (W=0)", m_DC), ("fold (W=0)", m_fold)]:
    vals = wmax_hits(f, K0S[::3])
    wmin, wmax = min(v[0] for v in vals), max(v[0] for v in vals)
    hits = sorted(set(v[1] for v in vals))
    R.append((f"c5/c6/c10 {name}: W_max(k0) in [{wmin:.6f},{wmax:.6f}], antipode clusters {hits}", ""))

# c7: pole reference k0 = pi on B1
d0 = m_B1(np.pi)
sp_arr = np.array([spsf(m_B1(k)) for k in K])
bg_arr = np.array([bragg(d0, m_B1(k)) for k in K])
R.append(("c7 pole k0=pi: Bragg vs SPSF pointwise (inter vs +, intra vs -)",
          f"{np.abs(bg_arr[:,0]-sp_arr[:,0]).max():.2e}, {np.abs(bg_arr[:,1]-sp_arr[:,1]).max():.2e}"))

# c8
wA1 = np.array([bragg_delta(m_A1(0.0), m_A1(k), k, 0.5) for k in K])
wA0 = np.array([bragg_delta(m_A0(0.0), m_A0(k), k, 0.5) for k in K])
R.append(("c8 uniform pair W=1 vs W=0 at delta=1/2: pattern diff", f"{np.abs(wA1-wA0).max():.2e}"))
for name, f in [("sshT(W=1)", m_sshT), ("ssh0(W=0)", m_ssh0)]:
    wmax = max(bragg_delta(f(0.0), f(k), k, 0.5) for k in KEDGE)
    R.append((f"c8 real {name} at delta=1/2 (edge-inclusive grid): max w", f"{wmax:.6f}"))

# c9
for name, f, exp in [("B1 (W=1)", m_B1, 2), ("R (W=0)", m_R, 2), ("fold (W=0)", m_fold, 2)]:
    dz = np.array([f(k)[2] for k in K])
    cross = int(np.sum(np.diff(np.sign(dz)) != 0))
    R.append((f"c9 SPSF pi-crossings {name}", f"{cross} (criterion reads nu={cross}/2; expect {exp})"))

# c11: T vs R minimal false-positive pair
spT = np.array([spsf(m_B1(k)) for k in K]); spR = np.array([spsf(m_R(k)) for k in K])
bgT_pole = np.array([bragg(m_B1(np.pi), m_B1(k))[0] for k in K])
bgR_pole = np.array([bragg(m_R(np.pi), m_R(k))[0] for k in K])
T_gen = max(bragg(m_B1(-np.pi/2), m_B1(k))[0] for k in K)
R_gen = max(bragg(m_R(-np.pi/2), m_R(k))[0] for k in K)
R.append(("c11 T vs R: SPSF diff, pole-Bragg diff, generic-k0 W_max (T vs R)",
          f"{np.abs(spT-spR).max():.2e}, {np.abs(bgT_pole-bgR_pole).max():.2e}, {T_gen:.6f} vs {R_gen:.6f}"))

# c12: rotation about z is an exact data symmetry; W is rotation-invariant
def rot_z(f, alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    def rf(k):
        d = f(k)
        return np.array([c*d[0]-s*d[1], s*d[0]+c*d[1], d[2]])
    return rf

for name, f in [("random", m_random), ("B1", m_B1)]:
    sp, bg = full_data(f, 0.0)
    spr, bgr = full_data(rot_z(f, 0.7), 0.0)
    dmax = max(np.abs(sp-spr).max(), np.abs(bg-bgr).max())
    extra = ""
    if name == "random":
        wS = wind_of(f, False); wR = wind_of(rot_z(f, 0.7), False)
        extra = f", W {wS} -> {wR}"
    R.append((f"c12 rotation R_z(0.7) {name}: data diff{extra}", f"{dmax:.2e}"))

# c13: delta-string rotation form + pole inertia for all delta
def rotz(v, beta):
    c, s = np.cos(beta), np.sin(beta)
    return np.array([c*v[0]-s*v[1], s*v[0]+c*v[1], v[2]])

e = 0.0
for name, f in [("random", m_random), ("A1", m_A1), ("B1", m_B1), ("fold", m_fold)]:
    for k0 in K0S[::4]:
        dh0 = f(k0)/np.linalg.norm(f(k0))
        for k in K[::13]:
            dh = f(k)/np.linalg.norm(f(k))
            w1 = bragg_delta(f(k0), f(k), k - k0, 0.37)
            w2 = 0.5*(1 - dh @ rotz(dh0, (k - k0)*0.37))
            e = max(e, abs(w1 - w2))
R.append(("c13 delta-string Rz form (4 models, delta=0.37): w_+ vs (1-dh.Rz(qd)dh0)/2", f"{e:.2e}"))

d0 = m_B1(np.pi)
e = 0.0
for delta in (0.11, 0.37, 0.73):
    for k in K[::13]:
        w = bragg_delta(d0, m_B1(k), k - np.pi, delta)
        wp, _ = spsf(m_B1(k))
        e = max(e, abs(w - wp))
R.append(("c13 pole k0=pi: Bragg(delta=0.11/0.37/0.73) = SPSF upper weight", f"{e:.2e}"))

for tag, val in R:
    print(f"{tag:76s} {val}")
