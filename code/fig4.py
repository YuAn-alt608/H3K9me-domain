# -*- coding: utf-8 -*-
"""Figure 4, redrawn from the source tables."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from data import r, style, FAM, ORANGE, BLUE, GREEN, DGREY, BM
import os

style(plt)

g = pd.read_csv(os.path.join(BM, "encoding_grid.tsv"), sep="\t")
ROWS = [("A", "accessibility", "A: DNase peak (binary)"),
        ("A", "H3K9me3", "A: H3K9me3"),
        ("B", "accessibility", "B: element/flank ratio"),
        ("B", "H3K9me3", "B: H3K9me3"),
        ("C", "accessibility (element)", "C: element DNase"),
        ("C", "accessibility (flank)", "C: neighbourhood DNase"),
        ("C", "H3K9me3", "C: H3K9me3")]

k5 = r("d8b_encoding_fit.tsv")
k5 = k5[k5.model == "C_split"].set_index("family").loc[FAM]
h2 = r("d3_hepg2_encoding_fit.tsv")
h2 = h2[h2.encoding.astype(str).str.startswith("C")].set_index("family").loc[FAM]

fig = plt.figure(figsize=(6.0, 2.97))
gs = fig.add_gridspec(1, 2, left=0.205, right=0.985, top=0.845, bottom=0.185,
                      wspace=0.24, width_ratios=[1.35, 1])

# ---------------- A three encodings ----------------
axA = fig.add_subplot(gs[0, 0])
n = len(ROWS)
ys = np.arange(n)[::-1]
axA.axvspan(0.28, 4.6, ymin=0, ymax=0.0, color="none")
axA.add_patch(plt.Rectangle((0.26, -0.48), 3.2, 2.96,
                            facecolor="#F1F1F1", edgecolor="none", zorder=0))
axA.axvline(1, color="#7A7A7A", lw=0.6, ls=(0, (3, 2.5)), zorder=1)
for i, (enc, term, lab) in enumerate(ROWS):
    y = ys[i]
    for cell, mk, col in [("K562", "o", BLUE), ("HepG2", "s", GREEN)]:
        row = g[(g.cell == cell) & (g.encoding == enc) & (g.label == term)].iloc[0]
        off = -0.13 if cell == "K562" else 0.13
        axA.plot([row.lo, row.hi], [y + off] * 2, color=col, lw=0.9,
                 solid_capstyle="butt", zorder=3)
        axA.plot([row.lo, row.lo], [y + off - 0.06, y + off + 0.06], color=col, lw=0.7)
        axA.plot([row.hi, row.hi], [y + off - 0.06, y + off + 0.06], color=col, lw=0.7)
        axA.plot([row.OR], [y + off], marker=mk, ms=3.6, color=col, zorder=4, lw=0)
axA.set_xscale("log")
axA.set_xlim(0.26, 3.45)
axA.set_ylim(-0.62, n - 0.38)
axA.set_yticks(ys)
axA.set_yticklabels([r[2] for r in ROWS], fontsize=6.2)
for lab, (enc, term, _) in zip(axA.get_yticklabels(), ROWS):
    if term == "H3K9me3":
        lab.set_color(ORANGE)
axA.set_xticks([0.3, 0.5, 1, 2, 3])
axA.set_xticklabels(["0.3", "0.5", "1", "2", "3"], fontsize=6.2)
axA.set_xlabel("odds ratio per SD (95% CI)", fontsize=6.4, labelpad=2)
axA.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axA.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["0.3", "0.5", "1", "2", "3"]))
axA.text(0.268, ys[5], "primary encoding", fontsize=5.4, color="#6E6E6E",
         rotation=90, va="center", ha="left")
h = [Line2D([0], [0], marker="o", color=BLUE, lw=0.9, ms=3.6, label="K562"),
     Line2D([0], [0], marker="s", color=GREEN, lw=0.9, ms=3.6, label="HepG2")]
axA.legend(handles=h, loc="lower right", bbox_to_anchor=(0.99, 0.01), fontsize=5.9,
           handlelength=1.3, borderpad=0.1, labelspacing=0.3)
axA.text(-0.52, 1.13, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.435, 1.13, "The ratio encoding flips the sign of accessibility",
         transform=axA.transAxes, fontsize=7.1, fontweight="bold", va="top")

# ---------------- B per-family mark coefficient ----------------
axB = fig.add_subplot(gs[0, 1])
x = np.arange(len(FAM))
w = 0.38
axB.bar(x - w / 2, k5.b_k9.values, width=w, color=BLUE, lw=0, label="K562")
axB.bar(x + w / 2, h2.b_k9.values, width=w, color=GREEN, lw=0, label="HepG2")
axB.axhline(0, color="#8C8C8C", lw=0.6)
axB.set_xticks(x)
axB.set_xticklabels(FAM, rotation=52, ha="right", fontsize=5.9)
axB.set_ylabel("H3K9me3 coefficient (log-odds per SD)", fontsize=6.4, labelpad=2)
axB.set_ylim(-1.02, 0.66)
axB.set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4])
axB.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=1, fontsize=5.9,
           handlelength=0.9, borderpad=0.1, labelspacing=0.25)
axB.text(0.99, 0.045, "K562 positive in 7/10; HepG2 positive in 4/10",
         transform=axB.transAxes, fontsize=5.7, color=DGREY, ha="right", va="bottom")
axB.text(-0.20, 1.13, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.13, 1.13, "The mark's contribution is K562-specific",
         transform=axB.transAxes, fontsize=7.1, fontweight="bold", va="top")

fig.savefig(os.path.join("..", "figures", "fig4.png"), dpi=450, facecolor="white")
print("fig4 written")
