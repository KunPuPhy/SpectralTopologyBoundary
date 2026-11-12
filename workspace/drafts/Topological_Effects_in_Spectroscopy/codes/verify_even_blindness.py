"""verify_even_blindness.py — numerical confirmation of the cos-evenness
blindness of intensity spectroscopy (SPSF and Bragg) in 1D chiral two-level
systems.  Companion note: notes/cos_even_blindness.tex.

Both intensity signals are cos-forms of SIGNED great-circle angles,
  SPSF : I_pm(k)     = pi (1 +- cos Phi(k))               [fixed z axis]
  Bragg: I_pm(k, k0) = (N/2)(1 +- cos(Phi(k) - Phi(k0)))  [ground-state axis]
and cos is even: the sign/orientation of the traversal never enters.

  (T1) z-fold blindness : a winding Phi(k) and its retrace |Phi(k)| (same
       d_z and r) have IDENTICAL diagonal SPSF — weights and pi-crossing
       readout (false positive: reads nu=1 at true nu=0).
  (T2) pole degeneracy  : from a pole reference (d_y(k0)=0, d_hat(k0)=+-z),
       cos Theta = +- d_hat_z(k): Bragg degenerates to the SPSF pattern and
       inherits its blindness (the natural BEC point k0=0 of the Kitaev
       family IS a pole).  Off the pole the z-fold pair IS resolved.
  (T3) global mirror    : d_y -> -d_y (nu -> -nu) leaves every intensity
       invariant: sign(nu) invisible to ALL intensity data, any reference.
  (T4) image != winding : the antipode test reads the image of d_hat, not
       the signed sweep: a full-circle forth-and-back path (nu=0, image=S^1)
       has W_max(k0) = 1 for every k0 (and 4 pi-crossings -> reads nu=2).

Unit-radius models on the y-z great circle (d_x = 0, r = 1, E_pm = +-1):
  T   Phi=-k               winding nu=-1    [pi -> 0 -> -pi]
  R   Phi=|k|              retrace nu=0     [pi -> 0 -> +pi]  (z-fold of T)
  T*  Phi=+k               mirror  nu=+1
  W2  Phi=2k+pi            genuine winding nu=+2
  DC  Phi=2k+pi | pi-2k    full circle twice: nu=0, image=S^1 (z-fold of W2)
  R2  d=(0, sin^2 k, cos k)  smooth fold of T: same d_z, r differs

Outputs: console verdicts + fig/even_blindness.pdf (six panels).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PI = np.pi
TWOPI = 2.0 * np.pi
K = np.linspace(-PI, PI, 1601)          # nodes at 0, +-pi/2, +-pi (1-ulp exact)
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fig")
TOL = 1.0e-12                            # float-rounding threshold for "identical"


# ----------------------------------------------------------------- models

def d_T(k):
    return np.stack([np.zeros_like(k), -np.sin(k), np.cos(k)], axis=-1)


def d_R(k):
    return np.stack([np.zeros_like(k), np.abs(np.sin(k)), np.cos(k)], axis=-1)


def d_Ts(k):
    return np.stack([np.zeros_like(k), np.sin(k), np.cos(k)], axis=-1)


def d_W2(k):
    return np.stack([np.zeros_like(k), np.sin(2 * k + PI), np.cos(2 * k + PI)], axis=-1)


def d_DC(k):
    phi = np.where(k < 0, 2 * k + PI, PI - 2 * k)
    return np.stack([np.zeros_like(k), np.sin(phi), np.cos(phi)], axis=-1)


def d_R2(k):
    return np.stack([np.zeros_like(k), np.sin(k) ** 2, np.cos(k)], axis=-1)


# ------------------------------------------------------------- utilities

def winding(d):
    """Signed sweep of Phi = atan2(d_y, d_z) along the great circle / 2pi."""
    th = np.unwrap(np.arctan2(d[:, 1], d[:, 2]))
    return (th[-1] - th[0]) / TWOPI


def dhat(d):
    return d / np.linalg.norm(d, axis=-1, keepdims=True)


def spsf_w(dfun):
    """(w_up, w_lo) = (1 +- d_hat_z)/2, resolved axis = z (Nambu particle)."""
    dh = dhat(dfun(K))
    return 0.5 * (1 + dh[:, 2]), 0.5 * (1 - dh[:, 2])


def bragg_w(dfun, k0):
    """(w_intra, w_inter) from the reference d_hat(k0): (1 +- cos Theta)/2.

    The reference is evaluated at the exact k0 (not at the nearest grid
    point), so pole references k0 = 0 are exact to the last bit."""
    dh = dhat(dfun(K))
    dh0 = dhat(dfun(np.atleast_1d(float(k0))))[0]
    c = dh @ dh0
    return 0.5 * (1 + c), 0.5 * (1 - c)


def wmax_scan(dfun):
    """W_max(k0) = max_k (1 - d_hat(k).d_hat(k0))/2  (antipode test)."""
    dh = dhat(dfun(K))
    return 0.5 * (1 - (dh @ dh.T).min(axis=1))


def pi_crossings(dfun):
    """Sign changes of d_z around the closed BZ loop (SPSF criterion)."""
    d = dfun(K)
    s = np.sign(d[:, 2])
    return int(np.count_nonzero(s[:-1] * s[1:] < 0) + (s[0] * s[-1] < 0))


def maxdev(a, b):
    return float(max(np.abs(a[0] - b[0]).max(), np.abs(a[1] - b[1]).max()))


def verdict(dev):
    return "IDENTICAL" if dev < TOL else "differ"


def hdr(msg):
    print("\n" + "=" * 74 + "\n" + msg + "\n" + "=" * 74)


# ------------------------------------------------------------------ main

def run_checks():
    hdr("ground truth: windings (signed sweep of Phi / 2pi)")
    for name, df in [("T  (Phi=-k,  pi -> 0 -> -pi)", d_T),
                     ("R  (Phi=|k|, pi -> 0 -> +pi)", d_R),
                     ("T* (Phi=+k,  mirror of T)", d_Ts),
                     ("W2 (Phi=2k+pi)", d_W2),
                     ("DC (full circle twice)", d_DC),
                     ("R2 (0, sin^2 k, cos k)", d_R2)]:
        print(f"  {name:30s}: nu = {winding(df(K)):+.1f}")

    hdr("T1  SPSF: winding T vs retrace R (same d_z, r=1) — blind")
    dev = maxdev(spsf_w(d_T), spsf_w(d_R))
    print(f"  max|Delta weights| = {dev:.1e}  -> {verdict(dev)}")
    print(f"  pi-crossings: T = {pi_crossings(d_T)}, R = {pi_crossings(d_R)}"
          f"  -> readout nu = 1 for BOTH (true: -1 vs 0): FALSE POSITIVE")

    hdr("T2  Bragg from pole references (d_y(k0)=0): k0 = 0, +-pi — blind")
    for k0 in (0.0, PI, -PI):
        dev = maxdev(bragg_w(d_T, k0), bragg_w(d_R, k0))
        extra = " (exact)" if dev == 0.0 else " (0 up to float rounding)"
        print(f"  k0 = {k0:+.2f}: max|Delta| = {dev:.1e}  -> {verdict(dev)}"
              f"{extra}: Bragg == SPSF pattern")

    hdr("T2' Bragg from generic reference k0 = -pi/2 — the pair IS resolved")
    bT, bR = bragg_w(d_T, -PI / 2), bragg_w(d_R, -PI / 2)
    i = int(np.argmax(np.abs(bT[1] - bR[1])))
    print(f"  max|Delta w_inter| = {np.abs(bT[1] - bR[1]).max():.3f} at k = {K[i]:+.3f}:"
          f"  T = {bT[1][i]:.2f} (Theta = pi, full transfer), R = {bR[1][i]:.2f} (Theta = 0)")
    print("  (weight difference ~ sin Phi(k) sin Phi(k0): vanishes iff k0 is a pole)")

    hdr("T3  global mirror T (nu=-1) vs T* (nu=+1): every intensity identical")
    dev_s = maxdev(spsf_w(d_T), spsf_w(d_Ts))
    dev_b = 0.0
    for k0 in K[::40]:
        dev_b = max(dev_b, maxdev(bragg_w(d_T, k0), bragg_w(d_Ts, k0)))
    print(f"  SPSF max|Delta| = {dev_s:.1e};  Bragg max|Delta| over k0-scan = {dev_b:.1e}")
    print("  -> sign(nu) invisible to all intensity data (mirror d_y -> -d_y)")

    hdr("T4  image != winding: DC (nu=0, image = full circle)")
    Wm = wmax_scan(d_DC)
    print(f"  DC: W_max(k0) in [{Wm.min():.8f}, {Wm.max():.8f}] for all k0"
          f"  -> = 1 up to grid resolution: antipode test fires at nu = 0")
    print(f"  DC: pi-crossings = {pi_crossings(d_DC)}  -> readout 'nu = 2'; true nu = 0")
    dev = maxdev(spsf_w(d_DC), spsf_w(d_W2))
    print(f"  DC vs W2 (genuine nu=2; DC is its z-fold): SPSF max|Delta| = {dev:.1e}"
          f"  -> {verdict(dev)}")
    for k0 in (0.0, PI / 2, -PI / 2, PI, -PI):
        dev = maxdev(bragg_w(d_DC, k0), bragg_w(d_W2, k0))
        print(f"    Bragg from pole k0 = {k0:+.2f}: max|Delta| = {dev:.1e}  -> {verdict(dev)}")
    dev = maxdev(bragg_w(d_DC, -PI / 4), bragg_w(d_W2, -PI / 4))
    print(f"    Bragg from generic k0 = -pi/4: max|Delta| = {dev:.3f}  (resolved)")

    hdr("T2'' retrace R: W_max(k0) = 1 exactly at the poles (blind == fooled)")
    WmR = wmax_scan(d_R)
    i0, ip, im, ig = (int(np.argmin(np.abs(K - x))) for x in (0.0, PI, -PI, -PI / 2))
    print(f"  R: W_max = {WmR[i0]:.6f} (k0=0), {WmR[ip]:.6f} (k0=+pi), {WmR[im]:.6f} (k0=-pi)"
          f"  |  {WmR[ig]:.6f} (k0=-pi/2, generic)")
    print("  -> R fires (=1) ONLY where it is blind (the poles): the false-fire")
    print("     points of the antipode test coincide with the blind points of T2")

    hdr("combined probes: SPSF + Bragg@generic k0 reconstruct sin Phi up to ONE global sign")
    k0 = -PI / 3
    dhT, dhTs = dhat(d_T(K)), dhat(d_Ts(K))
    i0 = int(np.argmin(np.abs(K - k0)))
    c1 = dhT[:, 2]                        # cos Phi(k)          (SPSF data)
    c2 = dhT @ dhT[i0]                    # cos(Phi - Phi(k0))  (Bragg data)
    sPhi0 = np.sqrt(1.0 - dhT[i0, 2] ** 2)
    rec = (c2 - c1 * dhT[i0, 2]) / sPhi0  # = sin Phi(k), branch sign fixed by hand
    err_T = np.abs(rec - dhT[:, 1]).max()
    c1s = dhTs[:, 2]
    c2s = dhTs @ dhTs[i0]
    same_data = max(np.abs(c1 - c1s).max(), np.abs(c2 - c2s).max())
    err_Ts = np.abs(rec - dhTs[:, 1]).max()
    print(f"  reconstruction error for T          : {err_T:.1e}  (exact)")
    print(f"  mirrored system T* gives same data  : max|Delta(c1,c2)| = {same_data:.1e},"
          f"  but |rec - sin Phi*| = {err_Ts:.3f}")
    print("  -> |nu| determinable in principle, sign(nu) NOT (the irreducible blind spot)")

    hdr("smooth fold R2 = (0, sin^2 k, cos k): criteria co-fire, weights differ")
    wS, wT_ = spsf_w(d_R2), spsf_w(d_T)
    dev = max(np.abs(wS[0] - wT_[0]).max(), np.abs(wS[1] - wT_[1]).max())
    r2 = np.linalg.norm(d_R2(K), axis=-1)
    iz = int(np.argmin(np.abs(K)))
    print(f"  SPSF weights vs T: max|Delta| = {dev:.3f}"
          f"  (r differs: max|1-r| = {np.abs(1 - r2).max():.3f} -> line positions differ)")
    print(f"  but pi-crossings: R2 = {pi_crossings(d_R2)} = T's {pi_crossings(d_T)};"
          f"  antipode test from k0=0: W_max = {wmax_scan(d_R2)[iz]:.6f}"
          f" = T's {wmax_scan(d_T)[iz]:.6f} — both fire")
    print("  -> smooth folds are not literally iso-spectral, but every intensity")
    print("     CRITERION (crossings, antipode test) co-fires: false positive persists")


# ---------------------------------------------------------------- figure

def style_k(ax, ylab=None, xlab=r"$k$"):
    ax.set_xlim(-PI, PI)
    ax.set_xticks([-PI, -PI / 2, 0, PI / 2, PI])
    ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(xlab, fontsize=11)
    if ylab:
        ax.set_ylabel(ylab, fontsize=11)
    ax.grid(alpha=0.3)


def path_arrow(ax, xy0, xy1, color, ls="-"):
    ax.annotate("", xy=xy1, xytext=xy0,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.9, ls=ls))


def panel_a(ax):
    th = np.linspace(0, TWOPI, 400)
    ax.plot(np.cos(th), np.sin(th), color="lightgray", lw=1.0, zorder=1)
    ax.axvline(0.0, color="gray", lw=0.9, ls="--", alpha=0.6, zorder=2)
    kk = np.linspace(-PI, PI, 400)
    ax.plot(-np.sin(kk), np.cos(kk), color="tab:blue", lw=2.8, zorder=5)
    for kc in (-0.85 * PI, -0.3, 0.55 * PI):
        path_arrow(ax, (-np.sin(kc), np.cos(kc)), (-np.sin(kc + 0.3), np.cos(kc + 0.3)), "tab:blue")
    kf = np.linspace(-PI, 0, 200)
    kb = np.linspace(0, PI, 200)
    ax.plot(np.abs(np.sin(kf)), np.cos(kf), color="tab:red", lw=2.8, zorder=6)
    ax.plot(np.abs(np.sin(kb)) + 0.06, np.cos(kb), color="tab:red", lw=1.7, ls="--", zorder=6)
    path_arrow(ax, (np.abs(np.sin(-PI / 2)), np.cos(-PI / 2)),
               (np.abs(np.sin(-PI / 2 + 0.3)), np.cos(-PI / 2 + 0.3)), "tab:red")
    path_arrow(ax, (np.abs(np.sin(PI / 2)) + 0.06, np.cos(PI / 2)),
               (np.abs(np.sin(PI / 2 - 0.3)) + 0.06, np.cos(PI / 2 - 0.3)), "tab:red", ls="--")
    ax.plot([0], [0], "+", color="k", ms=12, mew=1.6, zorder=8)
    ax.plot([0], [1], "o", color="k", ms=5, zorder=8)
    ax.plot([0], [-1], "o", color="k", ms=5, zorder=8)
    ax.plot([1], [0], "*", color="k", ms=13, zorder=8)
    ax.annotate(r"N: $\hat{\mathbf{d}}(0)$ (pole)", xy=(0.04, 1.02), fontsize=9)
    ax.annotate("S (pole)", xy=(0.05, -1.13), fontsize=9)
    ax.annotate(r"generic $k_0=-\pi/2$", xy=(0.62, 0.16), fontsize=9)
    ax.annotate("z-diameter\n(fold axis)", xy=(0.06, -0.66), fontsize=8.5, color="gray")
    ax.annotate(r"T: $\Phi=-k$, $\pi\to 0\to -\pi$, $\nu=-1$", xy=(-1.38, -0.30), fontsize=9, color="tab:blue")
    ax.annotate(r"R: $\Phi=|k|$, $\pi\to 0\to \pi$, $\nu=0$" + "\n(return offset for visibility)",
                xy=(0.30, -1.28), fontsize=9, color="tab:red")
    ax.set_xlabel(r"$d_y$", fontsize=11)
    ax.set_ylabel(r"$d_z$", fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-1.45, 1.62)
    ax.set_ylim(-1.42, 1.30)
    ax.set_title(r"(a) Bloch trajectories: T winds, R retraces (z-fold)", fontsize=10.5)


def panel_b(ax):
    ax.plot(K, PI * (1 + np.cos(K)), color="tab:red", lw=2.4, label=r"$I_+$ (R)")
    ax.plot(K, PI * (1 - np.cos(K)), color="tab:red", lw=2.4, ls="--", label=r"$I_-$ (R)")
    ax.plot(K, PI * (1 + np.cos(K)), color="tab:blue", lw=1.1, ls=(0, (4, 2)), label="T (on top)")
    ax.plot(K, PI * (1 - np.cos(K)), color="tab:blue", lw=1.1, ls=(0, (4, 2)))
    ax.axhline(PI, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.set_ylim(-0.3, TWOPI + 0.4)
    ax.set_yticks([0, PI, TWOPI])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    ax.text(0.0, 0.62, r"T $\equiv$ R: identical (max $|\Delta| = 0$)" + "\n"
            + "2 $\\pi$-crossings per band for both" + "\n"
            + r"$\Rightarrow$ SPSF reads $\nu=1$ (true: $-1$ vs $0$)",
            ha="center", fontsize=9)
    ax.legend(fontsize=8.5, loc="upper center", ncol=3)
    style_k(ax, ylab=r"$I_\pm(k)$")
    ax.set_title(r"(b) SPSF $I_\pm=\pi(1\pm\cos\Phi(k))$: even in $\Phi$", fontsize=10.5)


def panel_c(ax):
    ax.plot(K, 0.5 * (1 + np.cos(K)), color="tab:red", lw=2.4, label=r"$I_{\rm intra}$ (R)")
    ax.plot(K, 0.5 * (1 - np.cos(K)), color="tab:red", lw=2.4, ls="--", label=r"$I_{\rm inter}$ (R)")
    ax.plot(K, 0.5 * (1 + np.cos(K)), color="tab:blue", lw=1.1, ls=(0, (4, 2)), label="T (on top)")
    ax.plot(K, 0.5 * (1 - np.cos(K)), color="tab:blue", lw=1.1, ls=(0, (4, 2)))
    ax.set_ylim(-0.06, 1.12)
    ax.set_yticks([0, 0.5, 1.0])
    ax.text(0.0, 0.80, r"$\cos\Theta(k,0)=\hat d_z(k)$:" + "\n"
            + "Bragg @ pole $\\equiv$ SPSF pattern" + "\n"
            + r"$\Rightarrow$ same blindness as (b)",
            ha="center", fontsize=9)
    ax.legend(fontsize=8.5, loc="center left")
    style_k(ax, ylab=r"$I(k)/N$")
    ax.set_title(r"(c) Bragg from pole $k_0=0$: T $\equiv$ R (blind)", fontsize=10.5)


def panel_d(ax):
    bT, bR = bragg_w(d_T, -PI / 2), bragg_w(d_R, -PI / 2)
    ax.plot(K, bT[1], color="tab:blue", lw=2.4, label=r"$I_{\rm inter}$, T ($\nu=-1$)")
    ax.plot(K, bR[1], color="tab:red", lw=2.4, label=r"$I_{\rm inter}$, R ($\nu=0$)")
    ax.plot([PI / 2], [1.0], "o", color="tab:blue", ms=7)
    ax.plot([PI / 2], [0.0], "o", color="tab:red", ms=7)
    ax.annotate("", xy=(PI / 2 - 0.25, 0.86), xytext=(PI / 2 - 0.25, 0.18),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))
    ax.text(PI / 2 - 0.38, 0.50, r"$k=\pi/2$: T $\Theta=\pi$ (full)" + "\n" + r"R $\Theta=0$ (none)",
            ha="right", fontsize=9)
    ax.set_ylim(-0.06, 1.12)
    ax.set_yticks([0, 0.5, 1.0])
    ax.legend(fontsize=8.5, loc="center left")
    style_k(ax, ylab=r"$I_{\rm inter}(k)/N$")
    ax.set_title(r"(d) Bragg from generic $k_0=-\pi/2$: resolved", fontsize=10.5)


def panel_e(ax):
    WmT, WmR, WmDC = wmax_scan(d_T), wmax_scan(d_R), wmax_scan(d_DC)
    ax.plot(K, WmT, color="tab:blue", lw=2.6, label=r"T ($|\nu|=1$)")
    ax.plot(K, WmDC, color="dimgray", lw=1.5, ls="--", label=r"DC ($\nu=0$, image $=S^1$)")
    ax.plot(K, WmR, color="tab:red", lw=2.4, label=r"R ($\nu=0$, retrace)")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.7)
    for x in (0.0, PI, -PI):
        ax.axvline(x, color="gray", lw=0.8, ls=":", alpha=0.5)
    ax.text(0.0, 0.42, "R fires ($=1$) only at the poles" + "\n"
            + r"$k_0=0,\pm\pi$ — its blind points" + "\n"
            + "DC ($\\nu=0$) fires everywhere:" + "\n"
            + r"image $\neq$ winding",
            ha="center", fontsize=9)
    ax.set_ylim(0.0, 1.15)
    ax.set_yticks([0, 0.5, 1.0])
    ax.legend(fontsize=8.5, loc="lower center")
    style_k(ax, ylab=r"$W_{\max}(k_0)$", xlab=r"$k_0$")
    ax.set_title(r"(e) antipode test $W_{\max}(k_0)$: image $\neq$ winding", fontsize=10.5)


def panel_f(ax):
    bT, bTs = bragg_w(d_T, -PI / 2), bragg_w(d_Ts, -PI / 2)
    ax.plot(K, bT[1], color="tab:blue", lw=2.4, label=r"$I_{\rm inter}$, T ($\nu=-1$)")
    ax.plot(K, bTs[1], color="tab:green", lw=1.2, ls=(0, (4, 2)), label=r"T$^{*}$ ($\nu=+1$, on top)")
    ax.set_ylim(-0.06, 1.12)
    ax.set_yticks([0, 0.5, 1.0])
    ax.text(0.0, 0.78, r"T $\equiv$ T$^{*}$: identical from every $k_0$" + "\n"
            + r"$\Rightarrow$ sign$(\nu)$ invisible (global mirror)",
            ha="center", fontsize=9)
    ax.legend(fontsize=8.5, loc="center left")
    style_k(ax, ylab=r"$I_{\rm inter}(k)/N$")
    ax.set_title(r"(f) mirror pair $\nu=\pm1$: indistinguishable", fontsize=10.5)


def make_figure():
    fig = plt.figure(figsize=(13.2, 8.6))
    gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.32)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[0, 2]))
    panel_d(fig.add_subplot(gs[1, 0]))
    panel_e(fig.add_subplot(gs[1, 1]))
    panel_f(fig.add_subplot(gs[1, 2]))
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, "even_blindness.pdf")
    fig.savefig(path, dpi=300)
    print("\nsaved: " + os.path.normpath(path))
    plt.close(fig)


if __name__ == "__main__":
    run_checks()
    make_figure()
