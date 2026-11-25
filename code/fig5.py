# -*- coding: utf-8 -*-
"""Figure 5, redrawn from the source tables."""
import numpy as np
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from statsmodels.stats.proportion import proportion_confint
from data import r, style, ORANGE, BLUE, GREEN, GREY, DGREY, TAB

style(plt)

rg = r("rna_readout_groups.tsv")
ORDER = ["bound-open", "bound-closed", "unbound-open", "unbound-closed"]
LBL = {"bound-open": "marked\nopen", "bound-closed": "marked\nclosed",
       "unbound-open": "unmarked\nopen", "unbound-closed": "unmarked\nclosed"}
rg = rg.set_index("group").loc[ORDER]
if "pct" not in rg.columns:
    rg["pct"] = 100 * rg.rna_any
    lo, hi = proportion_confint(rg.rna_any * rg.n, rg.n, method="wilson")
    rg["lo"], rg["hi"] = 100 * lo, 100 * hi

pc = r("rna_positive_control.tsv")
genes = pc[pc.kind == "gene"].sort_values("rpm_per_kb", ascending=False)
rr = pc[pc.kind == "random"].rpm_per_kb

rob = r("robustness.tsv")
anti = rob[rob.result.str.startswith("D1 antisense")].iloc[0]

fig = plt.figure(figsize=(6.0, 2.72))
gs = fig.add_gridspec(1, 3, left=0.100, right=0.960, top=0.860, bottom=0.275,
                      wspace=0.58, width_ratios=[1.0, 1.12, 0.74])

# ---------------- A fraction of copies with a read ----------------
axA = fig.add_subplot(gs[0, 0])
x = np.arange(4)
cols = ["#C95B0C", "#C95B0C", "#8C8C8C", "#8C8C8C"]
hat = [None, "////", None, "////"]
for i in range(4):
    axA.bar(i, rg.pct.iloc[i], width=0.62, color=cols[i], lw=0,
            hatch=hat[i], edgecolor="white")
    axA.plot([i, i], [rg.lo.iloc[i], rg.hi.iloc[i]], color="#3A3A3A", lw=0.8, zorder=4)
axA.set_xticks(x)
axA.set_xticklabels(["%s\nn=%s" % (LBL[g].replace("\n", " "), format(int(rg.n.loc[g]), ","))
                    for g in ORDER], fontsize=5.9, rotation=32, ha="right",
                   rotation_mode="anchor")
axA.set_ylabel("copies with ≥ 1 read (%)", fontsize=6.4, labelpad=2)
axA.set_ylim(0, 24.5)
axA.set_yticks([0, 5, 10, 15, 20])
axA.set_xlim(-0.62, 3.62)
def _br(ax, x1, x2, y, h, label):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="#8A3B00", lw=0.55,
            clip_on=False, zorder=6)
    ax.text((x1 + x2) / 2.0, y + h + 0.35, label, ha="center", va="bottom",
            fontsize=5.6, color="#8A3B00", zorder=6)

_br(axA, 0.0, 1.0, 18.6, 0.7, "OR 0.91, p = 0.60")
_br(axA, 1.0, 3.0, 21.4, 0.7, "OR 1.45, p = 3.6×10$^{-5}$")
axA.tick_params(axis="x", pad=4)
axA.text(-0.30, 1.14, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.195, 1.14, "No transcription gain", transform=axA.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- B positive control ----------------
axB = fig.add_subplot(gs[0, 1])
names = list(genes.name) + ["200 random\n2 kb windows"]
vals = list(genes.rpm_per_kb) + [rr.mean()]
ys = np.arange(len(names))[::-1]
axB.barh(ys[:-1], vals[:-1], height=0.62, color=GREEN, lw=0)
for i, (n, v) in enumerate(zip(names[:-1], vals[:-1])):
    axB.text(v * 1.16, ys[i], "%.0f" % v if v >= 10 else "%.1f" % v, va="center",
             fontsize=5.6, color="#3A3A3A")
axB.plot([0.085, max(vals[:-1]) * 1.9], [ys[-1]] * 2, color="#B9B9B9", lw=0.5,
         ls=(0, (1.2, 1.6)))
axB.plot([vals[-1]] * 2, [ys[-1] - 0.16, ys[-1] + 0.16], color="#C95B0C", lw=1.6)
axB.text(0.115, ys[-1] - 0.42, "median 0, mean 0.11", fontsize=5.5, color="#8A3B00",
         va="top")
axB.set_xscale("log")
axB.set_xlim(0.085, 4200)
axB.set_ylim(-0.95, len(names) - 0.35)
axB.set_yticks(ys)
axB.set_yticklabels(names, fontsize=6.1)
axB.set_xlabel("RPM / kb", fontsize=6.4, labelpad=2)
axB.set_xticks([0.1, 1, 10, 100, 1000])
axB.set_xticklabels(["0.1", "1", "10", "100", "1000"], fontsize=6.1)
axB.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axB.text(-0.335, 1.14, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.255, 1.14, "Pipeline dynamic range", transform=axB.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- C strand test of the library ----------------
import pandas as pd
STR = pd.read_csv(os.path.join(TAB, "rna_strandedness.tsv"), sep="\t")
axC = fig.add_subplot(gs[0, 2])
xs = list(range(len(STR))) + [len(STR) + 0.75, len(STR) + 1.95]
vals = list(STR.sense_fraction) + [0.500, 0.500]
cols = ["#4A4A4A"] * len(STR) + [ORANGE, ORANGE]
for i, (v, c) in enumerate(zip(vals, cols)):
    axC.bar(i, v, width=0.62, color=c, lw=0)
axC.axhline(0.5, color="#8A8A8A", lw=0.6, ls=(0, (3, 2)))
axC.set_xticks(xs)
axC.set_xticklabels([g.replace("MALAT1", "MALAT1") for g in STR.gene] + ["elem.\nopen",
                                                                      "elem.\nclosed"],
                    rotation=50, ha="right", fontsize=5.0, rotation_mode="anchor")
axC.set_ylim(0, 1.12)
axC.set_yticks([0, 0.5, 1.0])
axC.set_xlim(-0.7, len(xs) + 0.45)
axC.set_ylabel("reads matching gene strand", fontsize=6.2, labelpad=2)
axC.tick_params(axis="x", pad=3)
axC.text(-0.42, 1.16, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top")
axC.text(-0.30, 1.16, "Library is unstranded", transform=axC.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

fig.savefig(os.path.join("..", "figures", "fig5.png"), dpi=450, facecolor="white")
print("fig5 written")
