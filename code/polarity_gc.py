# -*- coding: utf-8 -*-
"""5' versus 3' polarity of H3K9me3 and of DNase accessibility, with the paired
5'-3' difference in GC content regressed out.

Input
  data/tables/l1pa_windows.tsv          7,222 full-length L1PA copies, 5' window
                                        (first 1 kb) and 3' window (last 1 kb)
  data/tables/l1pa5_copy_tracks_v2.tsv.gz
                                        GC content of the same two windows
Output
  data/tables/l1pa_polarity_gc.tsv

Repeats the numbers quoted in the manuscript: DNase 98.6% of copies 5' > 3'
observed and 24.0% after the GC difference is removed (slope 0.69 per GC unit);
H3K9me3 57.7% observed and 92.9% after, so base composition accounts for the
accessibility polarity and not for the mark's.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "data", "tables")


def paired(df):
    p = df[df.part == "L1PA_5p"].set_index("unit_id")
    q = df[df.part == "L1PA_3p"].set_index("unit_id")
    return p.join(q, lsuffix="_5", rsuffix="_3", how="inner")


W = paired(pd.read_csv(os.path.join(TAB, "l1pa_windows.tsv"), sep="\t"))
G = paired(pd.read_csv(os.path.join(TAB, "l1pa5_copy_tracks_v2.tsv.gz"), sep="\t"))
j = W.join(G[["gc_5", "gc_3"]], how="inner")
print("paired copies", len(j))

rows = []
for track in ["h1_dnase", "h1_k9"]:
    y = (j[track + "_5"] - j[track + "_3"]).to_numpy()
    x = (j["gc_5"] - j["gc_3"]).to_numpy()
    ok = np.isfinite(y) & np.isfinite(x)
    y, x = y[ok], x[ok]
    b = float(np.dot(x - x.mean(), y - y.mean()) / np.dot(x - x.mean(), x - x.mean()))
    r = y - b * x
    rows += [
        dict(track=track, design="observed", n=len(y), slope=np.nan,
             share_5_gt_3=float((y > 0).mean()), median_diff=float(np.median(y))),
        dict(track=track, design="GC difference removed", n=len(y), slope=b,
             share_5_gt_3=float((r > 0).mean()), median_diff=float(np.median(r))),
    ]
    print("%-9s n=%d  observed %.1f%%  slope %+0.3f  residual %.1f%%  median resid %+0.3f"
          % (track, len(y), 100 * (y > 0).mean(), b, 100 * (r > 0).mean(), np.median(r)))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(TAB, "l1pa_polarity_gc.tsv"), sep="\t", index=False)
print("wrote", os.path.join(TAB, "l1pa_polarity_gc.tsv"))
