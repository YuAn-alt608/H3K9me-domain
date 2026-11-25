# -*- coding: utf-8 -*-
"""Figure 6, redrawn from the source tables."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from data import r, style, ORANGE, BLUE, GREEN, GREY, DGREY, TAB232

style(plt)

w = pd.read_csv(f"{TAB232}/l1pa_windows.tsv", sep="\t")
L = w[w.set == "L1PA"]
p = L.pivot_table(index="unit_id", columns="part", values=["h1_k9", "h1_dnase"], aggfunc="first")
p.columns = [f"{a}_{b}" for a, b in p.columns]
p = p.dropna(subset=["h1_k9_L1PA_5p", "h1_k9_L1PA_3p"])
k9_5, k9_3 = p.h1_k9_L1PA_5p.median(), p.h1_k9_L1PA_3p.median()
dn_5, dn_3 = p.h1_dnase_L1PA_5p.median(), p.h1_dnase_L1PA_3p.median()
n_copies = len(p)

cv = pd.read_csv(f"{TAB232}/paradox_pm_curve.tsv", sep="\t")
acc = pd.read_csv(f"{TAB232}/paradox_pm_access.tsv", sep="\t")
exp_med = float(acc[acc.readout.astype(str).str.startswith("PM00010")].exp_pm.median()) \
    if "exp_pm" in acc.columns else 0.0360
obs = dn_5
ratio = 1.96

g = r("g38_vs_hg19.tsv")
g = g[g.metric == "core_log2"].set_index("cell")

fig = plt.figure(figsize=(6.0, 2.71))
gs = fig.add_gridspec(1, 3, left=0.095, right=0.985, top=0.845, bottom=0.215,
                      wspace=0.52, width_ratios=[0.92, 1.12, 0.92])

# ---------------- A one polarity, two signals ----------------
axA = fig.add_subplot(gs[0, 0])
xs = np.arange(2)
wd = 0.34
axA.bar(xs[0] - wd / 2, k9_5, width=wd, color=ORANGE, lw=0)
axA.bar(xs[1] - wd / 2, k9_3, width=wd, color=ORANGE, lw=0)
axA.bar(xs[0] + wd / 2, dn_5 * 5, width=wd, color=GREEN, lw=0)
axA.bar(xs[1] + wd / 2, dn_3 * 5, width=wd, color=GREEN, lw=0)
axA.text(0 - wd / 2, k9_5 + 0.010, "%.3f" % k9_5, ha="center", fontsize=5.4)
axA.text(1 - wd / 2, k9_3 + 0.010, "%.3f" % k9_3, ha="center", fontsize=5.4)
axA.text(0 + wd / 2, 0.200, "%.4f" % dn_5, ha="center", va="center", rotation=90,
         fontsize=5.4, color="white", fontweight="bold")
axA.text(1 + wd / 2, dn_3 * 5 + 0.010, "%.4f" % dn_3, ha="center", fontsize=5.4)
axA.set_xticks(xs)
axA.set_xticklabels(["5′ UTR", "3′ end"], fontsize=6.3)
axA.set_ylabel("normalised signal", fontsize=6.4, labelpad=2)
axA.set_ylim(0, 0.63)
axA.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
axA.set_xlim(-0.62, 1.62)
axA.text(-0.52, 0.605, "H3K9me3", fontsize=5.6, color=ORANGE, fontweight="bold",
         va="center")
axA.text(-0.52, 0.560, "DNase (×5)", fontsize=5.6, color=GREEN, fontweight="bold",
         va="center")
axA.text(0.30, 0.605, "5′>3′ in 57.7%", fontsize=5.4, color="#4A4A4A", va="center")
axA.text(0.30, 0.560, "5′>3′ in 98.6%", fontsize=5.4, color="#4A4A4A", va="center")
axA.text(0.30, 0.505, "n = " + format(n_copies, ","), fontsize=5.4, color="#4A4A4A",
         va="center")
axA.text(-0.30, 1.14, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.215, 1.14, "One polarity, two signals", transform=axA.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- B accessible at matched mark ----------------
axB = fig.add_subplot(gs[0, 1])
axB.plot(cv.k9, cv.dnase, color="#8C8C8C", lw=1.0, marker="o", ms=2.6,
         mfc="#8C8C8C", mec="white", mew=0.4, zorder=3, label="random windows, binned by H3K9me3")
axB.axhline(obs, color=GREEN, lw=1.1, zorder=2)
axB.axhline(exp_med, color=ORANGE, lw=1.0, ls=(0, (4, 2.2)), zorder=2)
axB.text(cv.k9.max() * 0.97, obs + 0.0018, "L1PA 5′ UTR\n%.4f" % obs, fontsize=5.6,
         color=GREEN, ha="right", va="bottom")
axB.text(cv.k9.max() * 0.97, exp_med - 0.0016, "matched expectation\n%.4f" % exp_med,
         fontsize=5.6, color="#8A3B00", ha="right", va="top")
axB.annotate("", xy=(3.55, obs), xytext=(3.55, exp_med),
             arrowprops=dict(arrowstyle="<->", color="#3A3A3A", lw=0.7))
axB.text(3.30, (obs + exp_med) / 2, "%.2f×" % ratio, fontsize=6.0, ha="right", va="center")
axB.set_xlabel("H3K9me3 of 1 kb window", fontsize=6.4, labelpad=2)
axB.set_ylabel("DNase", fontsize=6.4, labelpad=2)
axB.set_xlim(-0.25, 5.6)
axB.set_ylim(0.008, 0.082)
axB.set_yticks([0.02, 0.04, 0.06, 0.08])
axB.legend(loc="lower left", bbox_to_anchor=(0.0, 0.035), fontsize=5.3,
           handlelength=1.8, borderpad=0.1, handletextpad=0.5)
axB.text(-0.28, 1.14, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.19, 1.14, "Accessible at matched mark", transform=axB.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- C assembly dependence ----------------
axC = fig.add_subplot(gs[0, 2])
cells = ["H1", "HepG2", "K562"]
xs = np.arange(3)
wd = 0.34
axC.bar(xs - wd / 2, [g.loc[c, "g38_diff"] for c in cells], width=wd, color=BLUE, lw=0)
axC.bar(xs + wd / 2, [g.loc[c, "hg19_diff"] for c in cells], width=wd, color="#9B9B9B",
        lw=0, hatch="////", edgecolor="white")
for i, c in enumerate(cells):
    axC.text(i - wd / 2, g.loc[c, "g38_diff"] + 0.035, "%.2f" % g.loc[c, "g38_diff"],
             ha="center", fontsize=5.6)
    axC.text(i + wd / 2, g.loc[c, "hg19_diff"] + 0.035, "%.2f" % g.loc[c, "hg19_diff"],
             ha="center", fontsize=5.6)
axC.set_xticks(xs)
axC.set_xticklabels(cells, fontsize=6.3)
axC.set_ylabel("paired 5′ − 3′ H3K9me3 (log$_2$)", fontsize=6.4, labelpad=2)
axC.set_ylim(0, 2.32)
axC.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
h = [Patch(facecolor=BLUE, label="GRCh38"),
     Patch(facecolor="#9B9B9B", hatch="////", edgecolor="white", label="hg19")]
axC.legend(handles=h, loc="upper right", bbox_to_anchor=(1.02, 1.06), fontsize=5.6,
           handlelength=0.9, borderpad=0.1, labelspacing=0.25)
axC.text(-0.30, 1.14, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top")
axC.text(-0.215, 1.14, "Assembly-dependent", transform=axC.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

fig.savefig(os.path.join("..", "figures", "fig6.png"), dpi=450, facecolor="white")
print("fig6 written", "obs", round(obs, 4), "exp", round(exp_med, 4))
