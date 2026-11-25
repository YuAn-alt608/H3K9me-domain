# H3K9me3 abundance does not identify the reader at transposable elements

Code and derived tables for a reanalysis of public occupancy data at ten
transposable-element families. The question is whether the abundance of
H3K9me3 on a copy predicts which reader protein occupies it: structural
readers (the HP1 paralogues CBX1, CBX3, CBX5) track the mark within every
family, an effector-coupled reader (TNRC18) does not, and at the copy level
what predicts co-repressor occupancy is accessibility rather than the mark.

This repository contains **no primary data**. Every input is a public
ENCODE, GEO or UCSC track; the ENCODE accessions are listed in the
manuscript's Materials and Methods. What is here is the code that turns those
tracks into the numbers in the figures, and the derived tables it reads.

## Layout

```
code/            analysis and plotting scripts
data/tables/     derived tables (one row per family, per track or per copy)
figures/         output of the plotting scripts, as used in the manuscript
```

## Reproducing the figures

```
pip install -r requirements.txt
cd code
python fig1.py           # -> ../figures/fig1.png
...
python fig7.py
python fig6_full.py      # -> ../figures/fig6_composition.png
python polarity_gc.py    # -> ../data/tables/l1pa_polarity_gc.tsv
```

Each script reads only from `data/tables/` and writes only to `figures/`
(or, for `polarity_gc.py`, back into `data/tables/`). No script downloads
anything; the derived tables are frozen here so that the numbers do not move
if a public track is re-released.

## What each script does

| Script | Figure / table | Reads |
|---|---|---|
| `fig1.py` | Figure 1, Table 1 | `d7_within_rho`, `d10_continuous`, `d7_matched_or`, `d7_peak_rho` |
| `fig2.py` | Figure 2 | `d7_copies`, `d9_reliability`, `d9_rho_corrected`, `d7_family_track`, `tnrc18_crosscell`, `genomewide_coloc`, `coloc_normalised`, `positive_control_coloc` |
| `fig3.py` | Figure 3, Table 2 | `k562_census_robust`, `k562_closure_by_family` |
| `fig4.py` | Figure 4 | `encoding_grid`, `d8b_encoding_fit`, `d3_hepg2_encoding_fit` |
| `fig5.py` | Figure 5, Table 3 | `rna_readout_groups`, `rna_positive_control`, `robustness`, `rna_strandedness` |
| `fig6.py` | Figure 6A, 6C | `l1pa_windows`, `g38_vs_hg19` |
| `fig6_full.py` | Figure 6 (three panels) | `l1pa5_pool_expansion` |
| `fig7.py` | Figure 7 | `perturbation_global`, `perturbation_descriptive`, `perturbation_depth`, `gfp_contrast`, `gfp_global` |
| `polarity_gc.py` | 5′–3′ polarity, text and Figure 6A | `l1pa_windows`, `l1pa5_copy_tracks_v2` (GC columns) |

`code/data.py` holds the shared loaders, the family and track order, and the
plot style, so that every panel of every figure is drawn on the same axis
conventions.

## What this repository deliberately does not do

`RESULTS.md` lists, claim by claim, what is carried in the manuscript and what
was computed during the analysis and then dropped, with the reason. Five
results were withdrawn in the course of the work — a mark-matched expectation
that did not survive a composition control, a per-copy residual that was never
recomputed within composition strata, a twelve-subfamily gradient, a
read-through inference from a library that turned out not to be
strand-specific, and a family-level percentile computed across two different
genome builds. They are documented there rather than deleted, because the
reason each one fails is the reason the surviving claims are stated the way
they are.

## Licence

MIT for the code (see `LICENSE`). The derived tables are released under CC0.
