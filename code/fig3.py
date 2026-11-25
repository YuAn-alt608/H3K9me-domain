# -*- coding: utf-8 -*-
"""Figure 3, redrawn from the source tables."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from data import r, style, ORANGE, BLUE, GREY, DGREY

style(plt)

rob = r("k562_census_robust.tsv")
m = rob[rob.family == "MER11A"].copy().sort_values("median_delta", ascending=False)
m["pt"] = 100 * m.median_delta
m = m.reset_index(drop=True)
MOD = ["SIN3A", "HDAC1", "HDAC2", "RCOR1", "MTA2", "MTA3", "GATAD2A", "KDM1A"]

cl = r("k562_closure_by_family.tsv").set_index("family")
BF = ["MER11A", "MER11B", "MER11C", "LTR12C", "MER41B"]

fig = plt.figure(figsize=(6.0, 2.82))
gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.845, bottom=0.215,
                      wspace=0.30, width_ratios=[1.32, 1])

# ---------------- A the 525-protein panel ----------------
axA = fig.add_subplot(gs[0, 0])
x = np.arange(len(m))
axA.axhline(0, color="#8C8C8C", lw=0.6)
axA.scatter(x, m.pt, s=3.2, color="#BFBFBF", edgecolor="none", zorder=2)
mod = m[m.target.isin(MOD)]
axA.scatter(mod.index, mod.pt, marker="v", s=16, facecolor=ORANGE, edgecolor="white",
            linewidth=0.4, zorder=4)
top = m.head(6)
block = "\n".join("%-8s %+.1f" % (row.target, row.pt) for _, row in top.iterrows())
axA.text(56, 10.55, block, fontsize=5.4, color="#3D3D3D", va="top", ha="left",
         family="DejaVu Sans", linespacing=1.45)
axA.annotate("co-repressor module\n(highest member KDM1A, +0.14)",
             xy=(mod.index.max(), 1.0), xytext=(len(m) * 0.78, -3.3),
             fontsize=5.4, color=ORANGE, ha="center", va="top",
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.5))
axA.set_xlabel("525 chromatin-associated proteins, ranked", fontsize=6.4, labelpad=2)
axA.set_ylabel("occupancy minus length-matched\nflanks (percentage points)", fontsize=6.4,
               labelpad=2)
axA.set_xlim(-14, len(m) + 14)
axA.set_ylim(-4.9, 10.9)
axA.set_yticks([-4, -2, 0, 2, 4, 6, 8, 10])
axA.set_xticks([])
axA.text(-0.075, 1.15, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.015, 1.15, "Adaptor and HP1 present, co-repressors absent",
         transform=axA.transAxes, fontsize=7.1, fontweight="bold", va="top")

# ---------------- B closure by family ----------------
axB = fig.add_subplot(gs[0, 1])
xs = np.arange(len(BF))
w = 0.27
hp1 = 100 * cl.loc[BF, "hp1_given_writer"].values
eff = 100 * cl.loc[BF, "eff_given_writer"].values
k9 = 50 * cl.loc[BF, "k9_med_writer"].values
axB.bar(xs - w, hp1, width=w, color=BLUE, lw=0)
axB.bar(xs, eff, width=w, color=ORANGE, lw=0)
axB.bar(xs + w, k9, width=w, color="#9B9B9B", lw=0, hatch="////", edgecolor="white")
for i, f in enumerate(BF):
    axB.text(i - w, hp1[i] + 1.5, "%.0f" % hp1[i], ha="center", fontsize=5.0, color=BLUE)
    axB.text(i, eff[i] + 1.5, "%.0f" % eff[i], ha="center", fontsize=5.0, color=ORANGE)
    axB.text(i + w, k9[i] + 1.5, "%.1f" % cl.loc[f, "k9_med_writer"], ha="center",
             fontsize=5.0, color="#5A5A5A")
axB.set_xticks(xs)
axB.set_xticklabels(["%s\nn=%d" % (f, cl.loc[f, "n_writer"]) for f in BF], fontsize=6.0)
axB.set_ylabel("% of adaptor-positive copies", fontsize=6.4, labelpad=2)
axB.set_ylim(0, 132)
axB.set_yticks([0, 20, 40, 60, 80, 100])
axB.tick_params(axis="x", pad=9)
h = [Patch(facecolor=BLUE, label="HP1 | adaptor$^+$"),
     Patch(facecolor=ORANGE, label="co-repressor | adaptor$^+$"),
     Patch(facecolor="#9B9B9B", hatch="////", edgecolor="white",
           label="median H3K9me3 (×50)")]
axB.legend(handles=h, loc="upper right", bbox_to_anchor=(1.01, 1.01), ncol=1,
           fontsize=5.5, handlelength=1.0, borderpad=0.1, labelspacing=0.25)
axB.text(-0.155, 1.15, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.085, 1.15, "Repeats across MER11",
         transform=axB.transAxes, fontsize=7.1, fontweight="bold", va="top")

fig.savefig(os.path.join("..", "figures", "fig3.png"), dpi=450, facecolor="white")
print("fig3 written")
