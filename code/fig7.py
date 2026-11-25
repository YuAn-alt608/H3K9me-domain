# -*- coding: utf-8 -*-
"""Figure 7 -- the perturbation tracks, described rather than tested.

Three panels, all descriptive:
  A  the global shift of every contrast (median mutant-minus-wild-type over all 33,246 copies)
     with its 100 kb block-bootstrap interval, next to the genome-wide depth of the two libraries
  B  the same shift per family, referenced to the global shift, with bootstrap intervals
  C  TNRC18 after the SETDB1 knockout, against the unperturbed GFP-tagged control of the same
     cells: three libraries per family, so that the knockout can be read against a condition in
     which SETDB1 is intact but no perturbation was made
No P value is drawn: one library per condition.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from data import style, ORANGE, BLUE, GREEN, GREY, DGREY, TAB

style(plt)

G = pd.read_csv(os.path.join(TAB, "perturbation_global.tsv"), sep="\t")
F = pd.read_csv(os.path.join(TAB, "perturbation_descriptive.tsv"), sep="\t")
D = pd.read_csv(os.path.join(TAB, "perturbation_depth.tsv"), sep="\t").set_index("track")
X = pd.read_csv(os.path.join(TAB, "gfp_contrast.tsv"), sep="\t")
XG = pd.read_csv(os.path.join(TAB, "gfp_global.tsv"), sep="\t")

ORDER = ["HDAC2", "TRIM28", "SETDB1", "H3K9me3", "TNRC18 KO"]
COL = {"HDAC2": ORANGE, "TRIM28": BLUE, "SETDB1": GREEN, "H3K9me3": "#8A8A8A",
       "TNRC18 KO": "#7A7A7A"}

fig = plt.figure(figsize=(6.5, 2.70))
gs = fig.add_gridspec(1, 3, left=0.100, right=0.990, top=0.775, bottom=0.235,
                      wspace=0.60, width_ratios=[1.0, 1.30, 1.14])

# ---------------- A global shift ----------------
axA = fig.add_subplot(gs[0, 0])
y = np.arange(len(ORDER))[::-1]
for i, k in enumerate(ORDER):
    r = G[G.contrast == k].iloc[0]
    axA.plot([r.global_lo, r.global_hi], [y[i], y[i]], color=COL[k], lw=1.4, solid_capstyle="round")
    axA.plot([r.global_median_delta], [y[i]], marker="o", ms=3.6, color=COL[k],
             mec="white", mew=0.6, zorder=4)
axA.axvline(0, color="#4A4A4A", lw=0.7, zorder=0)
axA.set_yticks(y)
axA.set_yticklabels(ORDER, fontsize=6.2)
axA.set_ylim(-0.6, len(ORDER) - 0.4)
axA.set_xlim(-0.32, 0.16)
axA.set_xticks([-0.3, -0.15, 0, 0.1])
axA.set_xlabel("median difference, mutant − WT (log$_2$)", fontsize=6.0, labelpad=1)
axA.text(-0.32, 1.255, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.22, 1.255, "Every contrast has a global offset", transform=axA.transAxes,
         fontsize=7.0, fontweight="bold", va="top")
axA.text(0.0, -0.30, "bars = 100 kb block bootstrap", transform=axA.transAxes,
         fontsize=5.2, color=DGREY)

# ---------------- B per family, referenced ----------------
axB = fig.add_subplot(gs[0, 1])
FAMS = ["L1PA3", "L1PA4", "L1PA5", "L1PA6", "LTR12C",
        "MER11A", "MER11B", "MER11C", "MER41A", "MER41B"]
tracks = [("H3K9me3", "#8A8A8A"), ("HDAC2", ORANGE), ("TRIM28", BLUE)]
x = np.arange(len(FAMS))
for j, (key, col) in enumerate(tracks):
    vals, lo, hi = [], [], []
    for f in FAMS:
        r = F[ (F.contrast == key) & (F.family == f) ]
        if not len(r):
            vals.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
        r = r.iloc[0]
        vals.append(r.delta_referenced)
        lo.append(r.delta_referenced - r.ref_lo)
        hi.append(r.ref_hi - r.delta_referenced)
    xpos = x + (j - 1) * 0.26
    axB.bar(xpos, vals, width=0.23, color=col, lw=0, label=key)
    axB.errorbar(xpos, vals, yerr=[lo, hi], fmt="none", ecolor="#3A3A3A", elinewidth=0.6,
                 capsize=1.4, capthick=0.6, zorder=5)
axB.axhline(0, color="#4A4A4A", lw=0.7)
axB.set_xticks(x)
axB.set_xticklabels(FAMS, rotation=42, ha="right", fontsize=5.4, rotation_mode="anchor")
axB.set_ylabel("difference, referenced to the global shift", fontsize=6.0, labelpad=2)
axB.set_ylim(-0.35, 1.45)
axB.set_yticks([-0.2, 0, 0.5, 1.0])
axB.legend(loc="upper left", bbox_to_anchor=(-0.03, 1.04), fontsize=5.6, ncol=3,
           handlelength=0.8, borderpad=0.1, columnspacing=0.7, frameon=False)
axB.text(-0.185, 1.255, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.105, 1.255, "The mark moves, the effectors do not",
         transform=axB.transAxes, fontsize=7.0, fontweight="bold", va="top")

# ---------------- C knockout, against an unperturbed control ----------------
axC = fig.add_subplot(gs[0, 2])
f3 = ["LTR12C", "MER11B", "MER11A", "L1PA3"]
ypos = np.arange(len(f3))[::-1]
wt, gfp, ko = [], [], []
for f in f3:
    r = F[(F.contrast == "TNRC18 KO") & (F.family == f)].iloc[0]
    g = X[(X.contrast == "TNRC18 GFP - endo") & (X.family == f)].iloc[0]
    wt.append(float(r.med_b))
    gfp.append(float(g.med_a))
    ko.append(float(r.med_a))
gshift = float(XG[XG.contrast == "TNRC18 KO - GFP"].global_median_delta.iloc[0])
SER = [(wt, -0.24, ORANGE, "parental"),
       (gfp, 0.0, "#9B7FB0", "TNRC18-GFP"),
       (ko, 0.24, "#8C8C8C", "SETDB1 KO")]
for vals, dy, col, lab in SER:
    axC.barh(ypos + dy, vals, height=0.22, color=col, lw=0, label=lab, zorder=3)
for vals, dy, col, lab in SER:
    for yy, v in zip(ypos + dy, vals):
        axC.text(v + (0.05 if v >= 0 else -0.05), yy, "%+.2f" % v,
                 va="center", ha="left" if v >= 0 else "right", fontsize=4.7,
                 color="#3A3A3A", zorder=4)
axC.axvline(0, color="#4A4A4A", lw=0.7)
axC.set_yticks(ypos)
axC.set_yticklabels(f3, fontsize=6.0)
axC.set_ylim(-0.62, len(f3) - 0.25)
axC.set_xlim(-1.15, 1.30)
axC.set_xticks([-1.0, -0.5, 0, 0.5, 1.0])
axC.set_xlabel("TNRC18 element-over-flank (log$_2$)", fontsize=6.0, labelpad=2)
axC.legend(loc="lower center", bbox_to_anchor=(0.46, 1.005), fontsize=4.2, ncol=3,
           handlelength=0.6, borderpad=0.0, columnspacing=0.7,
           handletextpad=0.35, frameon=False)
axC.text(-0.28, 1.30, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top")
axC.text(-0.20, 1.30, "Knockout vs unperturbed control", transform=axC.transAxes,
         fontsize=6.8, fontweight="bold", va="top")
axC.text(0.985, 0.04, "KO $-$ GFP: %+.2f" % gshift, transform=axC.transAxes, fontsize=5.0,
         color=DGREY, ha="right", va="bottom")

fig.savefig(os.path.join("..", "figures", "fig7.png"), dpi=450, facecolor="white")
print("fig7 written")
