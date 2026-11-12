"""Six-panel Bloch-sphere figure set for the 1D two-band models in the notes.

Each model produces one PDF (six stitched panels, 2 per row x 3 rows, PDF only):
  ssh_six_panels.pdf            chiral SSH chain, d = (t1 + t2 cos k, t2 sin k, 0)
  kitaev_six_panels.pdf         Kitaev chain BdG (Estake-Roy Eq. 9, J = Delta = 1),
                                 d = (0, Delta sin k, -mu/2 - J cos k)
  falsepositive_six_panels.pdf  Kitaev chain with d_y: sin k -> cos k (1D analog of
                                  Estake-Roy's construction, cf. Kitaev_dy2cosk_FalsePositive):
                                  winding destroyed (nu = 0 for all mu, loop = line segment,
                                  gapless only at mu = 0); d_z untouched -> SPSF pi-crossings
                                  identical to nu = 1 (false positive); d_hat image = arc,
                                  antipode NOT reached (Theta_max < pi) -> Bragg resolves it

Panels: (a) Bloch sphere, phase 1: raw trajectory d(k) + normalized d_hat(k)
        (b) Bloch sphere, phase 2
        (c) theta_k = angle between d_hat(k) and the d_z axis, vs k
        (d) Theta(k, k0) = angle between d_hat(k) and d_hat(k0), vs k   (k0 = 0)
        (e) SPSF line intensities   I_pm(k) = pi (1 +- cos theta_k)
        (f) Bragg line intensities  I(k) = N (1 +- cos Theta(k, k0)) / 2   (N = 1)
"""

import os

import matplotlib.pyplot as plt
import numpy as np

K0 = 0.0
NK = 801
K = np.linspace(-np.pi, np.pi, NK)

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fig")


def ssh_d(k, t1, t2):
    k = np.asarray(k, dtype=float)
    return np.column_stack([t1 + t2 * np.cos(k), t2 * np.sin(k), np.zeros_like(k)])


def kitaev_d(k, mu, J=1.0, Delta=1.0):
    k = np.asarray(k, dtype=float)
    return np.column_stack([np.zeros_like(k), Delta * np.sin(k), -mu / 2.0 - J * np.cos(k)])


def falsepositive_d(k, mu, J=1.0, Delta=1.0):
    k = np.asarray(k, dtype=float)
    return np.column_stack(
        [np.zeros_like(k), Delta * np.cos(k), -mu / 2.0 - J * np.cos(k)]
    )


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def angle_with_z(d):
    return np.arccos(np.clip(unit(d)[:, 2], -1.0, 1.0))


def angle_between(d, d_ref):
    return np.arccos(np.clip(unit(d) @ unit(d_ref), -1.0, 1.0))


def save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name + ".pdf")
    fig.savefig(path, dpi=300)
    print("saved: " + os.path.normpath(path))
    plt.close(fig)


def style_k_axis(ax, fs=14.0, xlabel=True):
    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks([-np.pi, -np.pi / 2, 0.0, np.pi / 2, np.pi])
    ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    if xlabel:
        ax.set_xlabel(r"$k$", fontsize=fs)
    ax.grid(alpha=0.3)


def draw_bloch_reference(ax):
    u = np.linspace(0.0, 2.0 * np.pi, 48)
    v = np.linspace(0.0, np.pi, 24)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color="lightsteelblue", alpha=0.10, linewidth=0, antialiased=False)
    th = np.linspace(0.0, 2.0 * np.pi, 240)
    ax.plot(np.cos(th), np.sin(th), 0.0 * th, color="gray", lw=0.8, ls="--", alpha=0.45)
    ax.plot(np.cos(th), 0.0 * th, np.sin(th), color="gray", lw=0.6, ls=":", alpha=0.30)
    ax.plot(0.0 * th, np.cos(th), np.sin(th), color="gray", lw=0.6, ls=":", alpha=0.30)


def draw_bloch_trajectory(ax, phase, lim, s, elev, azim):
    d = phase["d"]
    dhat = unit(d)
    draw_bloch_reference(ax)

    ax.plot(d[:, 0], d[:, 1], d[:, 2], color="tab:orange", lw=2.4, zorder=6, label=r"$\mathbf{d}(k)$")
    idx = np.linspace(0, NK - 1, 14).round().astype(int)
    ax.quiver(
        np.zeros(len(idx)), np.zeros(len(idx)), np.zeros(len(idx)),
        d[idx, 0], d[idx, 1], d[idx, 2],
        color="dimgray", alpha=0.5, lw=0.9, arrow_length_ratio=0.06, zorder=3,
    )
    ax.plot(
        dhat[:, 0], dhat[:, 1], dhat[:, 2],
        color="tab:blue", lw=3.2, zorder=10,
        label=r"$\hat{\mathbf{d}}(k)=\mathbf{d}(k)/|\mathbf{d}(k)|$",
    )
    tang = np.gradient(dhat, K, axis=0)
    for i in (NK // 8, 3 * NK // 8, 5 * NK // 8, 7 * NK // 8):
        t = 0.22 * tang[i] / np.linalg.norm(tang[i])
        ax.quiver(
            dhat[i, 0], dhat[i, 1], dhat[i, 2], t[0], t[1], t[2],
            color="tab:blue", lw=1.8, arrow_length_ratio=0.28, zorder=11,
        )

    ax.scatter([0.0], [0.0], [0.0], color="black", s=18, zorder=4)
    i0 = np.argmin(np.abs(K - K0))
    p0 = 1.18 * dhat[i0]
    ax.scatter([dhat[i0, 0]], [dhat[i0, 1]], [dhat[i0, 2]], color="black", marker="*", s=70, zorder=12)
    ax.text(p0[0] + 0.12, p0[1], p0[2], r"$\hat{\mathbf{d}}(k_0)$", fontsize=11 * s)

    ax.text2D(0.02, 0.02, phase["note3d"], transform=ax.transAxes, fontsize=11 * s)
    ax.set_title(phase["title"], fontsize=13 * s)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel(r"$d_x$")
    ax.set_ylabel(r"$d_y$")
    ax.set_zlabel(r"$d_z$")
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_zticks([-1, 0, 1])
    ax.tick_params(labelsize=9 * s)
    ax.view_init(elev=elev, azim=azim)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.0)
    ax.legend(loc="upper left", fontsize=11 * s)


def draw_theta_z(ax, data, s, xlabel, title, note=None, note_pos=(0.0, 0.0), note_box=False, legend_loc="lower center"):
    for ph in data:
        ax.plot(
            K, ph["theta_z"], color=ph["color"], lw=2.2, ls=ph["ls"],
            label=ph["legend"],
        )
    ax.set_ylim(0.0, np.pi)
    ax.set_yticks([0.0, np.pi / 2, np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])
    ax.set_ylabel(r"$\theta_k=\arccos\,\hat{d}_z(k)$", fontsize=13 * s)
    ax.set_title(title, fontsize=12.5 * s)
    if note is not None:
        kw = dict(ha="center", fontsize=11 * s)
        if note_box:
            kw["bbox"] = dict(boxstyle="round,pad=0.3", fc="ivory", ec="gray", alpha=0.9)
        ax.text(note_pos[0], note_pos[1], note, **kw)
    ax.legend(fontsize=11 * s, loc=legend_loc)
    style_k_axis(ax, fs=14 * s, xlabel=xlabel)


def draw_theta_k0(ax, data, s, xlabel, note_d=None):
    for ph in data:
        ax.plot(
            K, ph["Theta"], color=ph["color"], lw=2.2, ls=ph["ls"],
            label=ph["legend"],
        )
    ax.axhline(np.pi, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.text(-2.4, 2.76, r"$\Theta=\pi$ (antipode)", fontsize=10.5 * s, va="top", ha="left")
    ax.axvline(K0, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.text(K0 + 0.12, np.pi * 0.48, r"$k_0$", fontsize=11 * s)
    if note_d is not None:
        ax.text(note_d[1][0], note_d[1][1], note_d[0], fontsize=10.5 * s, ha="center")
    ax.set_ylim(0.0, np.pi * 1.06)
    ax.set_yticks([0.0, np.pi / 2, np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])
    ax.set_ylabel(r"$\Theta(k,k_0)=\arccos\,[\hat{\mathbf{d}}(k)\cdot\hat{\mathbf{d}}(k_0)]$", fontsize=12.5 * s)
    ax.set_title(r"Geodesic angle to the reference $\hat{\mathbf{d}}(k_0)$ ($k_0=0$)", fontsize=12.5 * s)
    ax.legend(fontsize=11 * s, loc="upper center")
    style_k_axis(ax, fs=14 * s, xlabel=xlabel)


def draw_I_spsf(ax, data, s, xlabel, title, note=None, note_pos=(0.0, 0.0), legend_loc="upper right"):
    for ph in data:
        cos_theta = np.cos(ph["theta_z"])
        ax.plot(K, np.pi * (1.0 + cos_theta), color=ph["color"], ls="-", lw=2.2,
                label="%s: $I_+$" % ph["short"])
        ax.plot(K, np.pi * (1.0 - cos_theta), color=ph["color"], ls="--", lw=2.2,
                label="%s: $I_-$" % ph["short"])
    ax.set_ylim(0.0, 2.0 * np.pi)
    ax.set_yticks([0.0, np.pi, 2.0 * np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    ax.set_ylabel(r"$I_\pm(k)$", fontsize=14 * s)
    ax.set_title(title, fontsize=12.5 * s)
    if note is not None:
        ax.text(
            note_pos[0], note_pos[1], note,
            ha="center", fontsize=11 * s,
            bbox=dict(boxstyle="round,pad=0.3", fc="ivory", ec="gray", alpha=0.9),
        )
    ax.legend(fontsize=11 * s, loc=legend_loc, ncol=2)
    style_k_axis(ax, fs=14 * s, xlabel=xlabel)


def draw_I_bragg(ax, data, s, xlabel, note_f=None, note_y=0.24):
    for ph in data:
        cos_theta = np.cos(ph["Theta"])
        ax.plot(K, 0.5 * (1.0 + cos_theta), color=ph["color"], ls="-", lw=2.2,
                label="%s: $I_{\\rm intra}$" % ph["short"])
        ax.plot(K, 0.5 * (1.0 - cos_theta), color=ph["color"], ls="--", lw=2.2,
                label="%s: $I_{\\rm inter}$" % ph["short"])
    ax.axvline(K0, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.text(K0 + 0.10, 0.90, r"$k_0$", fontsize=11 * s)
    if note_f is not None:
        ax.text(0.0, note_y, note_f, ha="center", fontsize=10 * s)
    ax.set_ylim(0.0, 1.05)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_ylabel(r"$I(k)/N$", fontsize=14 * s)
    ax.set_title(r"Bragg line intensities $I=\frac{N}{2}\,(1\pm\cos\Theta)$, $N=1$, $k_0=0$", fontsize=12.5 * s)
    ax.legend(fontsize=10 * s, loc="center", ncol=2)
    style_k_axis(ax, fs=14 * s, xlabel=xlabel)


MODELS = [
    dict(
        key="ssh",
        lim=1.62, elev=21, azim=-57,
        phases=[
            dict(dfunc=ssh_d, params=(0.5, 1.0), color="tab:blue", ls="-", nu=1, short="topo",
                 legend="topological ($\\nu=1$)",
                 title="SSH: topological phase ($t_1=0.5$, $t_2=1.0$, $\\nu=1$)",
                 note3d="$t_1=0.5,\\ t_2=1.0$\n$\\nu=1$: $\\hat{\\mathbf{d}}$ wraps the equator once"),
            dict(dfunc=ssh_d, params=(1.0, 0.5), color="tab:red", ls="--", nu=0, short="trivial",
                 legend="trivial ($\\nu=0$)",
                 title="SSH: trivial phase ($t_1=1.0$, $t_2=0.5$, $\\nu=0$)",
                 note3d="$t_1=1.0,\\ t_2=0.5$\n$\\nu=0$: $\\hat{\\mathbf{d}}$ oscillates back and forth (no wrapping)"),
        ],
        title_c=r"Angle with the $d_z$ axis: $\hat{\mathbf{d}}(k)$ pinned to the equator",
        note_c=(r"$d_z(k)\equiv 0\ \Rightarrow\ \theta_k\equiv\pi/2$ (both phases coincide)", (0.0, np.pi / 2 + 0.25), True),
        note_d=(r"trivial: $\Theta_{\max}=\arcsin(t_2/t_1)=\pi/6$", (0.0, 2.1)),
        title_e=r"SPSF line intensities $I_\pm=\pi\,(1\pm\cos\theta_k)$, sublattice $\alpha$ resolved",
        note_e=("\n".join([
            r"$\theta_k\equiv\pi/2\ \Rightarrow\ I_\pm\equiv\pi$:",
            r"the fixed-axis SPSF carries no winding information",
        ]), (0.0, 1.35)),
        legend_e="upper right",
        note_f=(r"topo: $I_{\rm inter}\to 1$ at $k=\pm\pi$", 0.24),
    ),
    dict(
        key="kitaev",
        lim=2.8, elev=18, azim=-30,
        phases=[
            dict(dfunc=kitaev_d, params=(1.0,), color="tab:blue", ls="-", nu=1, short="topo",
                 legend="topological ($\\nu=1$)",
                 title="Kitaev: topological ($\\mu=1<2J$, $J=\\Delta=1$, $\\nu=1$)",
                 note3d="$\\mu=1,\\ J=\\Delta=1$\n$\\nu=1$: $\\hat{\\mathbf{d}}$ wraps the $d_y$-$d_z$ great circle once"),
            dict(dfunc=kitaev_d, params=(3.0,), color="tab:red", ls="--", nu=0, short="trivial",
                 legend="trivial ($\\nu=0$)",
                 title="Kitaev: trivial ($\\mu=3>2J$, $J=\\Delta=1$, $\\nu=0$)",
                 note3d="$\\mu=3,\\ J=\\Delta=1$\n$\\nu=0$: $\\hat{\\mathbf{d}}$ stays in $z<0$ (no wrapping)"),
        ],
        title_c="\n".join([
            r"Angle with the $d_z$ axis ($\sin\tilde\theta_k=\hat d_z$, $\tilde\theta_k$ from the $d_y$ axis)",
            r"topo: $\theta_k\in[0,\pi]$ (both $d_z$ poles); trivial: $\theta_k>\pi/2$",
        ]),
        note_c=None,
        note_d=(r"trivial: $\Theta=\pi-\theta_k<\pi/2$ (no antipode)", (0.0, 2.05)),
        title_e="\n".join([
            r"SPSF line intensities $I_\pm=\pi(1\pm\sin\tilde\theta_k)$ (Nambu particle axis)",
            r"topo: 2 $\pi$-crossings per band ($=2\nu$); trivial: $I_+<\pi<I_-$",
        ]),
        note_e=None,
        legend_e="center",
        note_f=(r"topo: $I_{\rm inter}\to 1$ at $k=\pm\pi$", 0.16),
    ),
    dict(
        key="falsepositive",
        lim=2.8, elev=18, azim=-30,
        phases=[
            dict(dfunc=falsepositive_d, params=(1.0,), color="tab:blue", ls="-", nu=0, short="$\\mu<2J$",
                 legend="$\\mu<2J$ (false positive, $\\nu=0$)",
                 title="False positive: $\\mu=1<2J$ -- $\\nu=0$, SPSF mimics $\\nu=1$",
                 note3d="$d_y=\\cos k,\\ \\mu=1$\n$\\nu=0$: line-segment loop (retraced); south pole only, $\\Theta_{\\max}<\\pi$"),
            dict(dfunc=falsepositive_d, params=(3.0,), color="tab:red", ls="--", nu=0, short="$\\mu>2J$",
                 legend="$\\mu>2J$ (trivial, $\\nu=0$)",
                 title="False positive: $\\mu=3>2J$ -- trivial, $\\nu=0$",
                 note3d="$d_y=\\cos k,\\ \\mu=3$\n$\\nu=0$: $\\hat{\\mathbf{d}}$ stays in $z<0$"),
        ],
        title_c="\n".join([
            r"Angle with the $d_z$ axis ($\sin\tilde\theta_k=\hat d_z$, $\tilde\theta_k$ from the $d_y$ axis)",
            r"$\mu<2J$: $\theta_k\in[63^\circ,\pi]$ -- south pole only, narrower than the $\nu=1$ chain",
        ]),
        note_c=None,
        note_d=(r"$\mu<2J$: $\Theta_{\max}\approx 150^\circ<\pi$ -- antipode NOT reached (image = arc)", (0.0, 2.05)),
        title_e="\n".join([
            r"SPSF line intensities $I_\pm=\pi(1\pm\sin\tilde\theta_k)$ (Nambu particle axis)",
            r"false positive: 2 $\pi$-crossings per band at the same $k=\pm 2\pi/3$ as $\nu=1$ ($d_z$ untouched)",
        ]),
        note_e=None,
        legend_e="center",
        note_f=(r"no full transfer: $I_{\rm inter,max}\approx 0.93<1$ (Bragg resolves the false positive)", 0.24),
    ),
]


def build_model_data(model):
    data = []
    for ph in model["phases"]:
        d = ph["dfunc"](K, *ph["params"])
        d0 = ph["dfunc"](np.array([K0]), *ph["params"])[0]
        entry = dict(ph)
        entry["d"] = d
        entry["theta_z"] = angle_with_z(d)
        entry["Theta"] = angle_between(d, d0)
        data.append(entry)
    return data


def six_panel_figure(model):
    data = build_model_data(model)
    fig = plt.figure(figsize=(12.6, 16.8))
    gs = fig.add_gridspec(3, 2, height_ratios=(1.3, 1.0, 1.0), hspace=0.26, wspace=0.24)
    ax1 = fig.add_subplot(gs[0, 0], projection="3d")
    draw_bloch_trajectory(ax1, data[0], model["lim"], 0.85, model["elev"], model["azim"])
    ax2 = fig.add_subplot(gs[0, 1], projection="3d")
    draw_bloch_trajectory(ax2, data[1], model["lim"], 0.85, model["elev"], model["azim"])
    ax3 = fig.add_subplot(gs[1, 0])
    draw_theta_z(
        ax3, data, 0.9, xlabel=False, title=model["title_c"],
        note=model["note_c"][0] if model["note_c"] else None,
        note_pos=model["note_c"][1] if model["note_c"] else (0.0, 0.0),
        note_box=model["note_c"][2] if model["note_c"] else False,
    )
    ax4 = fig.add_subplot(gs[1, 1])
    draw_theta_k0(ax4, data, 0.9, xlabel=False, note_d=model["note_d"])
    ax5 = fig.add_subplot(gs[2, 0])
    draw_I_spsf(
        ax5, data, 0.9, xlabel=True, title=model["title_e"],
        note=model["note_e"][0] if model["note_e"] else None,
        note_pos=model["note_e"][1] if model["note_e"] else (0.0, 0.0),
        legend_loc=model["legend_e"],
    )
    ax6 = fig.add_subplot(gs[2, 1])
    draw_I_bragg(ax6, data, 0.9, xlabel=True, note_f=model["note_f"][0], note_y=model["note_f"][1])
    for ax, tag in zip((ax1, ax2, ax3, ax4, ax5, ax6), "abcdef"):
        ax.title.set_text("(%s) %s" % (tag, ax.title.get_text()))
    save(fig, model["key"] + "_six_panels")


def main():
    for model in MODELS:
        six_panel_figure(model)


if __name__ == "__main__":
    main()
