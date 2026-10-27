"""Figures. Categorical colours follow a fixed slot order (validated palette);
outcome cells also carry a letter so identity never relies on colour alone."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.stats import spearmanr

from .model import OUTCOMES
from .parameters import ALLOYS, COATINGS, PERFUSION_LEVELS

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRIDC = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OUT_COLOR = dict(zip(OUTCOMES, SERIES))
OUT_LETTER = dict(zip(OUTCOMES, "ABCDEF"))

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRIDC, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
    "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
})


def _legend_outcomes(fig, **kw):
    handles = [Patch(color=OUT_COLOR[o], label=f"{OUT_LETTER[o]}  {o}") for o in OUTCOMES]
    fig.legend(handles=handles, ncol=3, loc="lower center", **kw)


def time_courses(traces, path):
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.6))
    panels = [("pH0", "Interface pH", False), ("mg0", "Interface Mg$^{2+}$ (mM)", False),
              ("h2", "H$_2$ evolution (mL cm$^{-2}$ day$^{-1}$)", False),
              ("tumour", "Residual tumour burden (x initial)", True),
              ("bone", "Bone formation (x Ti control)", False),
              ("remaining", "Implant remaining (fraction)", False)]
    for ax, (key, title, logy) in zip(axs.flat, panels):
        for i, (c, tr) in enumerate(traces.items()):
            y = np.maximum(tr[key], 1e-12) if logy else tr[key]
            ax.plot(tr["t"], y, color=SERIES[i], lw=2, label=c)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(1e-12, 1e4)
        if key == "pH0":
            ax.axhspan(8.0, 8.4, color=SERIES[2], alpha=0.10, lw=0)
            ax.text(5, 8.2, "osteoblast optimum (Galow 2017)", color=INK2, fontsize=7, va="center")
        if key == "bone":
            ax.axhline(0.8, color=INK2, lw=0.8, ls="--")
        ax.set_title(title, loc="left")
        ax.set_xlabel("days after implantation")
    axs[0, 0].legend(loc="upper right", fontsize=8)
    fig.suptitle("WE43 porous implant, low perfusion, confined margin, pH-sensitive residual tumour (R1)",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def grid_heatmap(rows, path):
    phenos = ["pH-sensitive", "pH-resistant", "Mg-stimulated"]
    cols = [(pf, cf, fm) for cf in ("open", "confined") for fm in ("dense", "porous")
            for pf in PERFUSION_LEVELS]
    rlabels = [(a, c) for a in ALLOYS for c in COATINGS]
    cmap = ListedColormap([OUT_COLOR[o] for o in OUTCOMES])
    fig, axs = plt.subplots(1, 3, figsize=(16, 9.5), sharey=True)
    look = {(r["alloy"], r["coating"], r["perfusion"], r["confinement"], r["form"], r["phenotype"]): r["outcome"]
            for r in rows if r["residual"] == "R1-marginal"}
    for ax, ph in zip(axs, phenos):
        M = np.zeros((len(rlabels), len(cols)), dtype=int)
        for i, (a, c) in enumerate(rlabels):
            for j, (pf, cf, fm) in enumerate(cols):
                M[i, j] = OUTCOMES.index(look[(a, c, pf, cf, fm, ph)])
        ax.imshow(M, cmap=cmap, vmin=-0.5, vmax=len(OUTCOMES) - 0.5, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, "ABCDEF"[M[i, j]], ha="center", va="center", fontsize=6.5,
                        color="white" if M[i, j] in (0, 5) else INK)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{pf}\n{cf[:4]}/{fm[:4]}" for pf, cf, fm in cols], rotation=90, fontsize=6.5)
        ax.set_yticks(range(len(rlabels)))
        ax.set_yticklabels([f"{a} | {c}" for a, c in rlabels], fontsize=7)
        ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(rlabels)), minor=True)
        ax.grid(which="minor", color=SURFACE, lw=1.5)
        ax.grid(which="major", visible=False)
        ax.tick_params(which="minor", length=0)
        ax.set_title(f"{ph} tumour", loc="left")
    _legend_outcomes(fig, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Deterministic scenario grid (R1 marginal residual, nominal parameters). "
                 "Columns: perfusion / confinement / implant form", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def mc_bars(rows, path, key="outcome"):
    labels = [(a, c) for a in ALLOYS for c in COATINGS]
    fig, ax = plt.subplots(figsize=(10, 8.5))
    for i, (a, c) in enumerate(labels):
        sub = [r[key] for r in rows if r["alloy"] == a and r["coating"] == c]
        left = 0.0
        for o in OUTCOMES:
            w = sub.count(o) / len(sub)
            if w > 0:
                ax.barh(i, w, left=left, color=OUT_COLOR[o], edgecolor=SURFACE, linewidth=2, height=0.8)
                if w >= 0.08:
                    ax.text(left + w / 2, i, f"{100 * w:.0f}", ha="center", va="center", fontsize=6.5,
                            color="white" if o in (OUTCOMES[0], OUTCOMES[5]) else INK)
            left += w
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([f"{a} | {c}" for a, c in labels], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("probability across parameter uncertainty (%)")
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.grid(axis="y", visible=False)
    ax.set_title("Monte Carlo outcome probabilities (Latin hypercube, perfusion & tumour biology sampled)"
                 + ("\ntumour/bone outcome only - H$_2$ and structural limits ignored" if key != "outcome" else ""),
                 loc="left")
    _legend_outcomes(fig, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def window_map(crs, ks, grids, ph, path):
    cmap = ListedColormap([OUT_COLOR[o] for o in OUTCOMES])
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
    for ax, conf in zip(axs, ("open", "confined")):
        M = np.vectorize(OUTCOMES.index)(grids[conf]).astype(float)
        X, Y = np.meshgrid(ks, crs)
        ax.pcolormesh(X, Y, M, cmap=cmap, vmin=-0.5, vmax=len(OUTCOMES) - 0.5, shading="nearest")
        cs = ax.contour(X, Y, ph[conf], levels=[7.8, 8.4, 9.0, 10.0], colors=INK, linewidths=0.8)
        ax.clabel(cs, fmt="pH %.1f", fontsize=7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("perfusion clearance k (1/s)")
        ax.set_title(f"{conf} margin", loc="left")
        ax.grid(False)
    axs[0].set_ylabel("effective corrosion rate (mm/y) = alloy rate x coating factor x area")
    _legend_outcomes(fig, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Therapeutic window: degradation rate vs local clearance (no depletion / gas limits; "
                 "contours = peak interface pH)", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def sensitivity(rows, path, top=12):
    keys = [k for k in rows[0] if k.startswith("in_")]
    tum = np.array([max(r["log10_tumour_ratio"], -12.0) for r in rows])
    bone = np.array([r["bone_score_180"] for r in rows])
    res = []
    for k in keys:
        v = np.array([float(r[k]) for r in rows])
        if np.ptp(v) == 0:
            continue
        res.append((k[3:], spearmanr(v, tum)[0], spearmanr(v, bone)[0]))
    res.sort(key=lambda t: -max(abs(t[1]), abs(t[2])))
    res = res[:top]
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    names = [r[0] for r in res][::-1]
    for ax, idx, title in ((axs[0], 1, "log10 residual tumour ratio (lower = better control)"),
                           (axs[1], 2, "bone formation score at day 180")):
        vals = [r[idx] for r in res][::-1]
        ax.barh(names, vals, color=[SERIES[0] if v >= 0 else SERIES[1] for v in vals], height=0.6)
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_title(title, loc="left")
        ax.set_xlabel("Spearman rho")
        ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return res
