# -*- coding: utf-8 -*-
"""Figure 2, redrawn from the source tables."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from data import r, style, FAM, TRACKS, KEY, ORANGE, BLUE, GREEN, GREY, DGREY, TAB
from matplotlib.patches import Patch

style(plt)
rng = np.random.default_rng(11)

cp = pd.read_csv(r"D:/NewBeginning/BIOinformation/Age2.3/03_results/tables/d7_copies.tsv.gz", sep="\t")
cols = ["dnase", "k9"] + [KEY[t] for t in TRACKS]
sd = cp.groupby("repName")[cols].std(ddof=1).loc[FAM]
sd.columns = [{"dnase": "DNase", "k9": "H3K9me3"}.get(x, x) for x in sd.columns]

rel = r("d9_reliability.tsv").pivot(index="family", columns="track", values="reliability").loc[FAM]
rel.columns = [{"k9": "H3K9me3", "tnrc18": "TNRC18", "INPUT": "ChIP-seq input"}.get(c, c) for c in rel.columns]

rc = r("d9_rho_corrected.tsv").pivot(index="family", columns="partner", values="rho_obs").loc[FAM]
rc.columns = [{"k9": "H3K9me3", "dnase": "DNase", "INPUT": "ChIP-seq input"}.get(c, c) for c in rc.columns]

ft = r("d7_family_track.tsv").set_index("family").loc[FAM]
cc = r("tnrc18_crosscell.tsv")
k9 = cc[(cc.cell == "K562") & (cc.track == "tnrc18")].set_index("family").loc[FAM]
auc = k9.auc_k9
k9frac = 100 * ft.k9_pos / ft.n

cs = r("tnrc18_crosscell_summary.tsv")

fig = plt.figure(figsize=(6.0, 6.45))
gs = fig.add_gridspec(4, 2, left=0.155, right=0.985, top=0.940, bottom=0.078,
                      hspace=0.92, wspace=0.46, height_ratios=[1, 1, 0.92, 0.90])

# ---------------- A within-family SD ----------------
axA = fig.add_subplot(gs[0, 0])
orderA = sd.median().sort_values(ascending=False).index.tolist()
for i, t in enumerate(orderA):
    v = sd[t].values
    axA.scatter(i + rng.uniform(-0.13, 0.13, len(v)), v, s=7, facecolor="#BDBDBD",
                edgecolor="#A8A8A8", linewidth=0.25, zorder=2)
    axA.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color="#3A3A3A", lw=1.4,
             solid_capstyle="butt", zorder=3)
axA.set_xticks(range(len(orderA)))
axA.set_xticklabels([t.upper() if t == "tnrc18" else t for t in orderA],
                    rotation=48, ha="right", fontsize=6.2)
axA.get_xticklabels()[orderA.index("tnrc18")].set_color(ORANGE)
axA.set_ylabel("within-family SD (log$_2$)", fontsize=6.5, labelpad=2)
axA.set_ylim(0.45, 3.15)
axA.set_yticks([1.0, 1.5, 2.0, 2.5, 3.0])
axA.text(-0.23, 1.16, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.14, 1.16, "Dispersion within family", transform=axA.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- B split-half reliability ----------------
axB = fig.add_subplot(gs[0, 1])
orderB = rel.median().sort_values(ascending=False).index.tolist()
for i, t in enumerate(orderB):
    v = rel[t].values
    axB.scatter(i + rng.uniform(-0.13, 0.13, len(v)), v, s=7, facecolor="#BDBDBD",
                edgecolor="#A8A8A8", linewidth=0.25, zorder=2)
    axB.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color="#3A3A3A", lw=1.4,
             solid_capstyle="butt", zorder=3)
axB.axhline(0, color="#8C8C8C", lw=0.6)
axB.set_xticks(range(len(orderB)))
axB.set_xticklabels(orderB, rotation=48, ha="right", fontsize=6.2)
axB.get_xticklabels()[orderB.index("TNRC18")].set_color(ORANGE)
axB.set_ylabel("split-half reliability", fontsize=6.5, labelpad=2)
axB.set_ylim(-1.05, 1.12)
axB.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
axB.annotate("LTR12C\n(punctate)", xy=(orderB.index("TNRC18"), -0.89),
             xytext=(orderB.index("TNRC18") - 1.9, -0.80), fontsize=5.6, color=ORANGE,
             ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.6))
axB.text(-0.23, 1.16, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.14, 1.16, "Split-half reliability", transform=axB.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- C TNRC18 against every partner ----------------
axC = fig.add_subplot(gs[1, 0])
orderC = rc.median().sort_values().index.tolist()   # bottom = largest
for i, t in enumerate(orderC):
    v = rc[t].values
    axC.scatter(v, np.full(len(v), i) + rng.uniform(-0.13, 0.13, len(v)), s=7,
                facecolor="#BDBDBD", edgecolor="#A8A8A8", linewidth=0.25, zorder=2)
    axC.plot([np.median(v)] * 2, [i - 0.28, i + 0.28], color="#3A3A3A", lw=1.4,
             solid_capstyle="butt", zorder=3)
axC.axvline(0, color="#8C8C8C", lw=0.6)
axC.set_yticks(range(len(orderC)))
axC.set_yticklabels(orderC, fontsize=6.2)
for lab in axC.get_yticklabels():
    if lab.get_text() == "H3K9me3":
        lab.set_color(ORANGE)
axC.set_xlabel("within-family ρ with TNRC18", fontsize=6.5, labelpad=2)
axC.set_xlim(-0.20, 0.40)
axC.set_ylim(-0.7, len(orderC) - 0.3)
axC.text(-0.23, 1.14, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top")
axC.text(-0.14, 1.14, "Correlation with partner tracks", transform=axC.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- D AUC against marked fraction ----------------
axD = fig.add_subplot(gs[1, 1])
axD.axvline(0.5, color="#8C8C8C", lw=0.6, ls=(0, (3, 2.5)))
axD.text(0.5045, 49.5, "chance", fontsize=5.6, color="#6E6E6E", ha="left", va="top")
axD.scatter(auc.values, k9frac.values, s=16, facecolor=ORANGE, edgecolor="white",
            linewidth=0.5, zorder=3)
LAB = {"LTR12C": (-8, 9), "MER41B": (-26, 6), "MER11A": (-28, 4), "MER41A": (3, 6),
       "MER11B": (3, 6), "MER11C": (-32, -2), "L1PA5": (-26, 4), "L1PA4": (-34, 6),
       "L1PA6": (5, 2), "L1PA3": (-34, -7)}
for f in FAM:
    dx, dy = LAB[f]
    axD.annotate(f, (auc[f], k9frac[f]), xytext=(dx, dy), textcoords="offset points",
                 fontsize=5.5, color=DGREY, ha="left")
axD.set_xlabel("AUC of H3K9me3 for TNRC18 occupancy", fontsize=6.5, labelpad=2)
axD.set_ylabel("H3K9me3-positive\ncopies (%)", fontsize=6.2, labelpad=2)
axD.set_xlim(0.432, 0.608)
axD.set_ylim(-5, 53)
axD.set_yticks([0, 10, 20, 30, 40])
axD.text(-0.215, 1.14, "D", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")
axD.text(-0.105, 1.14, "Discrimination and mark prevalence", transform=axD.transAxes,
         fontsize=7.1, fontweight="bold", va="top")

# ---------------- E three cell lines ----------------
axE = fig.add_subplot(gs[2, :])
panel = [("K562\nTNRC18", cs[(cs.cell == "K562") & (cs.track == "tnrc18")].iloc[0], False),
         ("HeLa\nTNRC18", cs[(cs.cell == "HeLa") & (cs.track == "tnrc18")].iloc[0], False),
         ("HEK293\nTNRC18", cs[(cs.cell == "HEK293") & (cs.track == "tnrc18")].iloc[0], False),
         ("HEK293\nIgG", cs[(cs.cell == "HEK293") & (cs.track == "igg")].iloc[0], True)]
xs = np.arange(len(panel))
axE.bar(xs[:3], [p[1].rho_median for p in panel[:3]], width=0.52, color=ORANGE, lw=0)
axE.bar(xs[3], panel[3][1].rho_median, width=0.52, color="#9A9A9A", lw=0,
        hatch="////", edgecolor="white")
for i, (lab, row, igg) in enumerate(panel):
    axE.plot([i, i], [row.rho_min, row.rho_max], color="#3A3A3A", lw=0.8, zorder=4)
    axE.text(i, row.rho_median + 0.012, "%.3f" % row.rho_median, ha="center",
             va="bottom", fontsize=6.0, color="#222222")
axE.set_xticks(xs)
axE.set_xticklabels([p[0] for p in panel], fontsize=6.3)
axE.set_ylabel("ρ with H3K9me3", fontsize=6.5, labelpad=2)
axE.set_ylim(0, 0.34)
axE.set_yticks([0, 0.1, 0.2, 0.3])
axE.annotate("", xy=(2.74, 0.176), xytext=(2.26, 0.155),
             arrowprops=dict(arrowstyle="-", color="#6E6E6E", lw=0.6))
axE.text(1.30, 0.300, "same antibody series; TNRC18 does not exceed its own IgG",
         fontsize=5.9, color=DGREY, ha="left", va="center")
axE.set_xlim(-0.6, 3.6)
axE.text(-0.048, 1.16, "E", transform=axE.transAxes, fontsize=9, fontweight="bold", va="top")
axE.text(0.005, 1.16, "Three cell lines, matched IgG",
         transform=axE.transAxes, fontsize=7.1, fontweight="bold", va="top")

# ---------------- F genome-wide co-localisation ----------------
CL = pd.read_csv(os.path.join(TAB, "genomewide_coloc.tsv"), sep="\t")
NM = pd.read_csv(os.path.join(TAB, "coloc_normalised.tsv"), sep="\t")
axF = fig.add_subplot(gs[3, :])
PC = pd.read_csv(os.path.join(TAB, "positive_control_coloc.tsv"), sep="\t")
scales = [10, 100, 1000]
xpos = np.arange(3)


def raw_vals(track):
    return [PC[(PC.scale_kb == s) & (PC.track_a == track) & (PC.track_b == "H3K9me3")
               & (PC.normalised == "raw")].spearman.iloc[0] for s in scales]


def norm_vals(track):
    return [PC[(PC.scale_kb == s) & (PC.track_a == track) & (PC.track_b == "H3K9me3")
               & (PC.normalised == "IgG-normalised")].spearman.iloc[0] for s in scales]


groups = [("K562 TNRC18", ORANGE), ("HEK293 TNRC18", BLUE), ("HEK293 TRIM28", GREEN)]
for j, (lab, col) in enumerate(groups):
    if lab.startswith("K562"):
        vals = [CL[(CL.bin_kb == s) & (CL.track_a == "K562_TNRC18") &
                   (CL.track_b == "K562_H3K9me3")].spearman.iloc[0] for s in scales]
        nm = NM[NM.cell == "K562"].set_index("measure").spearman
        nvals = [nm["signal / control"]] * 3
    elif lab.startswith("HEK293 TNRC18"):
        vals = [CL[(CL.bin_kb == s) & (CL.track_a == "HEK293_TNRC18") &
                   (CL.track_b == "HEK293_H3K9me3")].spearman.iloc[0] for s in scales]
        nm = NM[NM.cell == "HEK293"].set_index("measure").spearman
        nvals = [nm["signal / control"]] * 3
    else:
        vals = raw_vals("TRIM28")
        nvals = norm_vals("TRIM28")
    xo = xpos + (j - 1) * 0.26
    axF.bar(xo, vals, width=0.23, color=col, lw=0,
            label=lab.replace("HEK293 ", "").replace("K562 ", "K562 "))
    for x, v in zip(xo, vals):
        axF.text(x, v / 2 if v < 0 else v + 0.012, "%.2f" % v, ha="center",
                 va="center" if v < 0 else "bottom", fontsize=5.2, color="white" if v < 0 else col,
                 fontweight="bold")
    axF.plot(xo, nvals, marker="o", ms=3.0, mfc="white", mec=col, mew=0.8, lw=0)
axF.axhline(0, color="#4A4A4A", lw=0.7)
axF.set_xticks(xpos)
axF.set_xticklabels(["10 kb", "100 kb", "1 Mb"], fontsize=6.4)
axF.set_xlim(-0.68, 2.68)
axF.set_ylim(-0.46, 0.26)
axF.set_yticks([-0.4, -0.3, -0.2, -0.1, 0, 0.2])
axF.set_ylabel("ρ (track vs H3K9me3)", fontsize=6.4, labelpad=2)
axF.set_xlabel("bin size", fontsize=6.5, labelpad=1)
from matplotlib.patches import Patch as _P
axF.legend(handles=[_P(facecolor=ORANGE, label="K562 TNRC18"),
                    _P(facecolor=BLUE, label="HEK293 TNRC18"),
                    _P(facecolor=GREEN, label="HEK293 TRIM28")],
           loc="lower left", bbox_to_anchor=(0.0, 1.00), fontsize=5.4, ncol=3,
           handlelength=0.8, borderpad=0.1, columnspacing=0.8, frameon=False)
axF.text(-0.045, 1.46, "F", transform=axF.transAxes, fontsize=9, fontweight="bold", va="top")
axF.text(0.000, 1.46, "Genome-wide co-localisation",
         transform=axF.transAxes, fontsize=7.1, fontweight="bold", va="top")
axF.text(0.985, 0.05, "bars: raw   circles: control-normalised", transform=axF.transAxes,
         fontsize=5.3, color=DGREY, ha="right", va="bottom")

fig.savefig("fig2.png", dpi=450, facecolor="white")
print("fig2 written")
