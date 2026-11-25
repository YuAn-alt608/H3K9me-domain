# What is carried, and what was dropped

Every number in the manuscript is produced by a script in `code/` from a table
in `data/tables/`. This file records the decisions behind that list: which
results survived the controls the paper applies, and which were computed at
some point during the analysis and then dropped. Nothing here is a
supplementary claim; it is the audit trail behind the seven figures.

---

## 1. Carried

| # | Claim | Number | Script | Table |
|---|---|---|---|---|
| 1 | Within every family, H3K9me3 abundance predicts the occupancy of CBX1, CBX3 and CBX5, and not that of TNRC18 | mean within-family ρ 0.512 (CBX3), 0.455 (CBX5), 0.417 (CBX1) against 0.060 (TNRC18); Table 1 places TNRC18 at the 0th percentile of 36 random proteins | `fig1.py` | `d7_within_rho`, `d7_family_track` |
| 2 | The contrast holds at matched positive rate and at 300 bp resolution | TNRC18 lowest in 9 of 10 families; ρ −0.05 to 0.09 against 0.20 to 0.65 for CBX3 | `fig1.py` | `d7_matched_or`, `d7_peak_rho` |
| 3 | The TNRC18 track is punctate, reproducible at about half the mark's reliability, and does not exceed its matched IgG | split-half median 0.44 (mark 0.91); HEK293 0.149 against IgG 0.173 | `fig2.py` | `d9_reliability`, `tnrc18_crosscell` |
| 4 | Every track in the deposit follows family-level accessibility, so the family-level contrast is unusable | ρ with DNase 1.00 (H3K9me3), 0.79 (TNRC18), 0.73 (POLR2A) | `fig2.py` | `d7_family_track` |
| 5 | On MER11A the adaptor and HP1 are present together and the co-repressor module is not | 187 of 854 copies adaptor-positive, HP1 on 67.4% of them, module on 15.5% | `fig3.py` | `k562_census_robust`, `k562_closure_by_family` |
| 6 | Copy-level co-repressor occupancy is predicted by accessibility at the element and in its neighbourhood, in two cell lines; the mark term is specific to K562 | neighbourhood OR 2.70 per SD, element 1.28; mark 1.31 in K562 and 0.94 in HepG2 (CI upper bound 1.08) | `fig4.py` | `encoding_grid`, `d8b_encoding_fit`, `d3_hepg2_encoding_fit` |
| 7 | Accessibility at marked copies does not buy transcription | marked-open against marked-closed OR 0.91 (P = 0.60); marked-closed against unmarked-closed 1.45 (P = 3.6 × 10⁻⁵) | `fig5.py` | `rna_readout_groups`, `robustness` |
| 8 | The library cannot arbitrate by strand | 0.475–0.515 of reads match the gene strand at five expressed genes | `fig5.py` | `rna_strandedness` |
| 9 | The L1PA 5′ UTR excess over matched windows is base composition | 1.96-fold matched on the mark alone, 1.06 (1.01–1.10) within H3K9me3 × GC strata, 0.76–1.27 across five designs on an independent 120,000-window pool | `fig6.py`, `fig6_full.py` | `l1pa_windows`, `l1pa5_pool_expansion` |
| 10 | The mark's 5′–3′ polarity survives the same GC control, the accessibility polarity does not | DNase 98.6% → 24.0% of copies 5′ > 3′; H3K9me3 57.7% → 92.9% | `polarity_gc.py` | `l1pa_windows`, `l1pa5_copy_tracks_v2` |
| 11 | The perturbation tracks are described, not tested | mark moves +0.54 to +1.26 log₂ at MER11 and MER41, co-repressors −0.11 to +0.17 | `fig7.py` | `perturbation_global`, `perturbation_descriptive` |

Two conventions that make the table readable. Every within-family statistic is
computed on the same copies with the same metric (element mean over flank mean,
ε = 0.01, identical windows), so a comparison between tracks is a comparison of
the same quantity. Every family-level statistic is reported as a description
only, because all tracks follow family-level accessibility and ten points
cannot separate the mark from openness.

## 2. Dropped

Each of these was computed, and each is absent from the manuscript for a stated
reason. They are listed in the order the analysis reached them.

1. **The mark-matched expectation, 1.96-fold (ratio of medians) and 2.29-fold
   (median of per-copy ratios), with 88.5% of copies above their own
   expectation.** Matched on H3K9me3 alone, the 5′ UTR is compared against
   windows with GC 0.405 against its own 0.545. Within H3K9me3 × GC strata the
   ratio falls to 1.06 (1.01–1.10), and the share above expectation to 34–67%.
   The comparison measured base composition. Only the composition-matched
   estimate is carried.

2. **Per-copy conservation of the mark-matched residual between H1 and H9
   (ρ = 0.733, n = 7,144).** A residual defined against the expectation that
   item 1 removes cannot be tested on a second cell line; it was never
   recomputed within composition strata, and is not reported as a result.

3. **A gradient across twelve subfamilies, 1.53- to 2.93-fold above matched
   expectation.** Same defect as item 2, and the same consequence.

4. **Read-through from the strand composition of element-internal reads**
   (median antisense fraction 0.500). The library (ENCFF001RPY) is not
   strand-specific: at five highly expressed genes the share of reads matching
   the annotated strand is 0.475–0.515. A fraction of 0.5 inside the elements
   is a property of the library. Withdrawn; the section now reports the strand
   test as a negative control rather than as evidence.

5. **A family-level percentile of the co-repressor module computed against 36
   proteins scored in a different build** (GRCh38 copies, union of IDR peaks,
   against hg19 element-over-flank positivity). The two numbers were not the
   same quantity. Recomputed with one build, one copy universe and one
   definition of occupancy: 49th percentile of the 495 panel proteins, 18th of
   4,000 random eight-protein unions. The recomputed values are carried; the
   original percentile is not.

6. **Negative-binomial coefficients from a fit without a length offset, and
   from a subset restricted to copies carrying a mappability value**
   (0.25 and 1.77). The offset is required: without it the length term loads
   onto mappability. Both parameterisations are tabulated
   (Supplementary Tables S1a and S1b) and neither carries a claim.

## 3. Boundaries

These are not dropped results but limits that apply to the carried ones, and
they are stated in the manuscript rather than here: one library per
perturbation condition, so no P value is attached to Figure 7; the K562 and
HepG2 arms differ in assembly, mark track and panel, so the mark term is
reported as cell-line specific; the occupancy panels contain few
effector-coupled readers, so a reader-class property cannot be separated from
a TNRC18-specific one in public data; and the L1PA 5′ UTR magnitudes depend on
how control windows are matched and on the reference assembly, which is why
every estimate in that section is reported against a composition-matched
expectation.
