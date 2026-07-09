# Which amplifications reach the protein: a gene-intrinsic account of copy-number-to-protein transmission, with application to antibody-drug conjugate target selection in 3q-amplified non-small-cell lung cancer

## Abstract

Somatic copy-number amplification is pervasive in cancer, but its functional
consequence is realised only if the extra gene dosage reaches the protein, and that
transmission is attenuated gene by gene. Using matched proteogenomic profiles across
six tumour types, we build a leakage-safe model that predicts, per case and per gene,
whether an amplified gene is over-expressed at the protein level, and we ask what the
model is detecting. On the highly amplified 3q cohort the model reaches an area under
the precision-recall curve of 0.7515 against a base rate of 0.4191 (lift 1.79); on the
pan-cancer cohort it reaches 0.6465 against 0.3983 (lift 1.62). Its distinctive value
is not average discrimination but resolution within the highly amplified regime, where
copy number alone is uninformative. Four independent lines of evidence then show that
the transmission the model detects is a stable, gene-intrinsic property. It transfers
across tumour lineages (mean pairwise rank concordance 0.840, Kendall W 0.867). It is
predictable from external gene properties alone, with no protein-derived feature, at a
leave-gene-out rank correlation of 0.547. It is not positional: holding out entire
chromosome arms changes that correlation by 0.001. And its mechanism reconciles with
post-transcriptional buffering, dominated by gene dosage-sensitivity and protein
biophysics. Adding a protein-protein interaction network improves prediction only by
widening gene coverage, not because interaction topology is intrinsically informative.
Applied to antibody-drug conjugate target selection in 3q-amplified non-small-cell lung
cancer, and filtered on five orthogonal axes, this framework yields three surface
targets that survive every axis (ATP13A3, ATP1B3, PLSCR1) and holds back two that are
high on transmission but carry a specific liability (TFRC, essential; ITGB5,
lineage-unstable). About two thirds of per-gene transmission variance remains
unexplained by gene properties, so the framework is a gene-intrinsic prior over the
genome, not a per-gene oracle; we make that bound explicit with a kept-in error, AP2M1.

All model metrics are reported on unique (case, gene) pairs, deduplicated, with the
random-precision base rate beside every area under the precision-recall curve as lift.
Seed 2 throughout; British spelling; numeric section references.

---

## 1. Introduction

Somatic copy-number alteration is among the most pervasive classes of genomic change in
cancer, spanning whole-chromosome aneuploidy through arm-level and focal events
(Beroukhim et al. 2010, 10.1038/nature08822; Zack et al. 2013, 10.1038/ng.2760).
Recurrent amplifications draw particular interest because they recur across tumours on
regions harbouring oncogenic drivers, and because the genes they carry are, in
principle, dosage-elevated and therefore candidate dependencies and therapeutic targets
(Sanchez-Vega et al. 2018, 10.1016/j.cell.2018.03.035). The functional consequence of
an amplification, however, is realised at the level of protein, and the path from gene
dosage to protein is neither direct nor uniform. An account that stops at the DNA, or
even at the transcript, can mislead about which genes are actually over-produced.

That path is attenuated. Gene dosage propagates towards protein through transcription
and then translation, and the signal is buffered at more than one step, so that protein
abundance is partially decoupled from copy number across many genes (Goncalves et al.
2017, 10.1016/j.cels.2017.08.013). A well-characterised contributor is the
stoichiometric control of protein complexes, in which a subunit produced in excess of
its partners is degraded rather than accumulated (Stingele et al. 2012,
10.1038/msb.2012.40; Sousa et al. 2019, 10.1074/mcp.RA118.001280). Crucially, the
degree of buffering is gene-specific: some amplifications transmit fully to protein
while others are largely buffered (Vogel and Marcotte 2012, 10.1038/nrg3185), and
recent work resolves how buffering depends on complex membership and dosage-sensitivity
(Heller et al. 2026, 10.1038/s44320-026-00187-9). The practical implication is that
copy number alone is a weak guide to which amplified genes are over-expressed at the
protein level, which is the level that determines function and therapeutic
accessibility.

Large proteogenomic efforts now quantify copy number, RNA and protein in common tumours
across several tumour types (Gillette et al. 2020, 10.1016/j.cell.2020.06.013; Li et
al. 2023, 10.1016/j.ccell.2023.06.009), making it possible to ask, gene by gene and
case by case, whether an amplification reaches the protein, and to model the
determinants of transmission from data rather than to assume them.

The translational setting is amplification of chromosome arm 3q, a recurrent event in
squamous malignancies and in non-small-cell lung cancer more broadly (McCaughan et al.
2010, 10.1164/rccm.201001-0005OC; Jeon et al. 2023, 10.1111/1759-7714.15045). The
amplified region is large and contains many co-amplified genes whose individual
transmission to protein is not established. Because so many genes are amplified
together, 3q-amplified non-small-cell lung cancer is a context in which distinguishing
the amplifications that reach surface protein from those that do not is both difficult
and directly useful. That usefulness is sharpest for antibody-drug conjugates, an
expanding modality (Nelson et al. 2023, 10.1146/annurev-med-071322-065903) whose
viability depends on the target being abundantly and selectively present as protein on
the cell surface (Gazzah et al. 2022, 10.1016/j.annonc.2021.12.012). An amplified gene
is an attractive target only if its amplification produces surface protein; selecting
on copy number or transcript alone risks committing a programme to genes whose protein
is attenuated.

This work develops an attenuation-aware framework that predicts protein-level
transmission within an amplified region and gates the result on surface localisation.
Its central claim, established here from four independent directions, is that
transmission is a stable, gene-intrinsic property, and that this is precisely what makes
a transmission-aware prior useful for target selection where copy number is blind.

---

## 2. Methods

### 2.1 Data and cohorts

The framework is built on matched proteogenomic profiles (copy number, RNA, protein)
across six tumour types (LSCC, LUAD, GBM, PDA, CCRCC, UCEC), with an independent cohort
reserved for transfer assessment. Two analysis cohorts are used throughout: a
pan-cancer cohort (ALL) and a 3q-amplified cohort (3Q). Amplification is defined as a
ploidy-adjusted copy number above one together with an integer copy number of at least
three (the R4 definition), the regime in which copy number alone ceases to discriminate
protein outcome.

### 2.2 Outcome and leakage-safe design

The outcome is protein-derived: whether an amplified gene is over-expressed at the
protein level in a given case. No protein-derived feature is ever used as a predictor.
Splits are grouped by case and fixed across all experiments for comparability; any
value propagated over a network or neighbourhood is computed from training-fold cases
only and recomputed per fold, so that validation cases never inform their own features.
This design is verified rather than asserted: a per-fold leakage self-check confirms
zero train-test feature overlap on every fold-aware feature.

### 2.3 Model, evaluation and deduplication

The model is a gradient-boosted decision ensemble (XGBoost, 600 trees, single-thread,
seed 2), with feature groups chosen by forward group selection. Five-fold stacking
produces byte-identical duplicate (case, gene) rows that inflate both the base rate and
the area under the precision-recall curve, so every metric is computed after
deduplication to unique (case, gene) pairs. The random-precision base rate is reported
beside every area under the precision-recall curve as lift, and precision-at-K per case
and per-subtype lift are used for cross-cohort reads.

### 2.4 The four gene-intrinsic experiments

Four post-hoc experiments, none of which retrains the model beyond a leave-tumour-out
refit, test whether the transmission the model detects is gene-intrinsic: a transfer
experiment across the six lineages (Section 3.2); a gene-intrinsic transmissibility
predictor built from external features alone (Section 3.3); a positional control that
holds out whole chromosome arms (Section 3.4); and a network-integration experiment
that asks whether a protein-protein interaction network adds anything beyond gene
coverage (Section 3.6). The mechanism is then reconciled with an independent
post-transcriptional buffering analysis (Section 3.5).

---

## 3. Results

### 3.1 A transmission model whose value is in the amplified regime

On the R4 amplified cohorts the model discriminates protein-level over-expression well
above the base rate. In the pan-cancer ALL cohort it reaches an area under the
precision-recall curve of 0.6465 against a base rate of 0.3983 (lift 1.62); in the
3q-amplified 3Q cohort it reaches 0.7515 against 0.4191 (lift 1.79). Per-tumour lifts
are consistent across lineages (Table 2, Table 1), ranging in ALL from 1.50 (PDA) to
1.71 (LSCC). The point of the model is not that it beats copy number on average, but
that it resolves which amplifications reach the protein inside the highly amplified
regime, where copy number is uninformative by construction. That resolution is only
worth trusting if it reflects something stable about the genes rather than an artefact
of any one cohort. The next four sections show that it does.

### 3.2 Transmission transfers across tumour lineages

Leaving out one tumour type at a time and ranking genes by predicted transmission in
the held-out lineage, the per-gene rankings agree across the six lineages at a mean
pairwise rank concordance of 0.840 (range 0.782 to 0.901) and a Kendall coefficient of
concordance of 0.867 (Figure 18a; Table 6). Against a null that holds the models fixed
and permutes the outcome label, the observed concordance stands far outside the null
(standardised distance +271.5; permutation p at the resampling floor). A gene's
transmission propensity is therefore not a lineage-specific accident; it transfers. The
median across-lineage stability ratio is 3.6, and 990 of 6,715 genes (14.7 per cent)
are strictly stable at a ratio of at most two. A genuinely context-dependent tail of
633 genes (9.4 per cent) moves substantially by lineage; these are named rather than
hidden (Table 6, Figure 4), and they matter for target selection because a
lineage-unstable gene is a poorer target than its pan-cancer rank suggests.

### 3.3 Transmission is predictable from gene properties alone

If transmission is gene-intrinsic, it should be predictable from properties of the gene
measured independently of any proteomics. Using only external features (dosage
sensitivity, protein biophysics, mRNA, evolution, complex membership, network
centrality, expression breadth) and a leave-gene-out design, the predictor recovers the
observed transmissibility at a rank correlation of 0.547 and an R-squared of 0.322
(Figure 18b; Figure 7; Table 14). No protein-derived feature enters the predictor; the
outcome is the protein-derived quantity. This is the strongest form of the
gene-intrinsic claim: a gene's coupling can be anticipated before its protein is
measured. The predictor is deliberately conservative on the one borderline feature
group (expression breadth), whose removal lowers the correlation only to 0.518.

### 3.4 The signal is not positional

Transmissibility could in principle be an artefact of genomic position, a gene inheriting
its apparent coupling from strong co-amplified neighbours. A leave-chromosome-arm-out
control, which holds out all arms in turn so no gene is predicted using a same-arm
neighbour, gives a rank correlation of 0.5462 against the leave-gene-out 0.5474, a drop
of 0.001 (Figure 18c; Table 15). The coupling is a property of the gene, not of its
address. This control is what licenses treating the predictor's output as a gene-level
prior rather than a positional summary.

### 3.5 The mechanism reconciles with post-transcriptional buffering

TreeSHAP attributes the prediction to protein biophysics and dosage-sensitivity first
(biophysics 31 per cent, dosage 24 per cent of grouped attribution; mRNA 15 per cent,
expression breadth and network 9 to 10 per cent each, evolution 4 per cent, complex 1
per cent) (Figure 18d; Figure 8; Table 16). Dosage-sensitivity dominates exactly as an
independent buffering analysis predicts: on the residual readout of post-transcriptional
buffering, dosage-sensitivity features add +0.104 in rank terms over baseline, almost
three times the +0.036 contributed by complex membership, and once dosage annotations
are present the additional value of complex membership collapses to +0.0019. The
apparent essentiality sign is a coding convention, not a contradiction: the
negative-is-essential encoding correlates -0.345 with transmissibility and the
positive-is-essential encoding +0.375, both meaning that more essential genes transmit
more, the same biology the buffering analysis reports at -0.29 to -0.33. The two
analyses, built on different readouts, converge on gene dosage-sensitivity as the
operative property.

### 3.6 The network adds coverage, not a new mechanism

A protein-protein interaction network improves prediction in both cohorts when stacked
on the incumbent model (STRING stack: ALL +0.057, 3Q +0.049 in area under the
precision-recall curve; Table 18; Figure 15), and all improvements clear a
10,000-permutation label-shuffle null with the models held fixed (Table 20). But once
the network is restricted to the genes the incumbent KEGG neighbourhood already covers,
the advantage nearly vanishes (ALL +0.004; 3Q -0.015, no evidence of improvement). The
relationship is monotonic in coverage: the sparser the network, the more the swap arms
lose, and the weakest network (OmniPath, 43 to 46 per cent coverage) loses outright
(Table 20; Figure 17; Table 23). The network helps by widening the gene set for which
any neighbourhood feature can be computed, not because interaction topology is
intrinsically more informative than pathway membership. This is itself a gene-intrinsic
result: the coupling does not need the network; the network is a coverage device. The
co-prediction structure within 3q is consistent with this reading, showing coherent
complex-level co-regulation among 3q genes broadly (CORUM enrichment, 19 of 21
complexes at FDR below 0.05 in the 3Q cohort; Table 17) rather than a tight module
among any small set of named genes.

### 3.7 The join between the prior and the transfer result

The transmissibility predictor and the transfer experiment measure two different
quantities and must not be conflated. The predictor estimates a gene's average
transmission across the pan-cancer cohort; the transfer experiment measures the
variance of that transmission across lineages. A gene can therefore be high on average
yet lineage-unstable, with no inconsistency between the two results. ITGB5 is the
worked example: the predictor ranks it high (percentile 0.71), while the transfer
experiment places it in the context-dependent tail (stability ratio 6.5). Both readings
are correct; they answer different questions, and a target-selection argument needs both.

---

## 4. Application: antibody-drug conjugate target selection in 3q-amplified non-small-cell lung cancer

The four evidence lines are not an end in themselves; they are the filter that turns a
ranked list of amplified surface proteins into a defensible shortlist. Each of the 24
hardened 3q surface candidates is scored on five orthogonal axes at once (Table 26): the
gene-intrinsic transmissibility prior, the observed model transmission, network
isolation, transfer context-dependence, and buffering or essentiality liability.

The prior is neither better nor worse on the surface proteins that matter here than it
is genome-wide: rank accuracy on the 24 candidates is 0.553, essentially identical to
the whole-atlas 0.547. That equivalence is important because it means the honest bound
carries over directly. An R-squared of 0.322 leaves roughly two thirds of per-gene
transmission variance unexplained by gene properties, and for 6 of the 24 candidates
that residual flips the call relative to the prior. The prior is a gene-intrinsic prior
over the genome, not a per-gene oracle.

AP2M1 is the clearest single illustration of that limit, and it is more convincing than
the caveat stated in the abstract because it is an error kept in the table rather than
removed from it. AP2M1 ranks high enough on the prior to be noticed (predicted
percentile 0.68) and is an essential gene, yet its observed transmission is only 0.10:
the prior over-predicts it substantially. A per-gene oracle would not do this; a
gene-intrinsic prior that captures average behaviour and misses gene-specific
regulation will, and pointing to a kept-in over-prediction is the most direct evidence
of the bound.

Applying the corroboration guard, of the five candidates the prior labels high on both
predicted and observed transmission, three survive every axis: ATP13A3, ATP1B3 and
PLSCR1 (Table 26; Figure 16). These carry none of the five flags. The word for them is
proportionate: they are the leads that survive every axis, not validated targets.
Consistent with the honest bound, these are transcript-level and prediction-level calls,
not surface-protein measurements, and about two thirds of per-gene variance remains
unexplained; each would require direct measurement of surface-protein abundance and
tumour-versus-normal selectivity before it could be called a target.

The remaining two prior-high candidates are held back deliberately, each for a specific
and different reason. TFRC has the highest predicted transmissibility of all 24
(percentile 0.78) and the highest observed transmission (0.74); it is genuinely the
strongest transmitter in the set. But its DepMap dependency effect of -0.94 marks it as
broadly essential, an on-target-toxicity liability that narrows its therapeutic window
regardless of how well it transmits. ITGB5 is corroborated high by the prior
(percentile 0.71) yet sits in the context-dependent transfer tail (stability ratio
6.5): its high pan-cancer rank is not a promise of stability in any given lineage.
Corroboration by the prior is not the same as clean, and the guard exists precisely to
stop a strong transmitter being silently re-promoted to a lead. Twelve further
candidates are network-isolated among 3q genes, which reinforces rather than weakens
their gene-intrinsic reading: their transmission is not a neighbourhood effect. The
honest translational position is three leads that survive every filter, with TFRC and
ITGB5 named and held back, not quietly dropped or quietly restored.

---

## 5. Discussion

The result of this work is a single claim with a translational corollary. The claim is
that copy-number-to-protein transmission is a stable, gene-intrinsic, positional-free
property: it transfers across lineages, it is predictable from gene biology measured
without any proteomics, it survives holding out whole chromosome arms, and it reconciles
mechanistically with post-transcriptional buffering through gene dosage-sensitivity. The
four experiments approach the claim from independent directions and converge, which is
stronger evidence than any one of them alone, because their failure modes do not
overlap: transfer could hold while prediction failed, prediction could hold on
positional confounding, and none of these would survive if the signal were a
cohort-specific artefact.

The corollary is that a transmission-aware prior is most useful exactly where copy
number is blind. In the highly amplified 3q region, where many genes are co-amplified
and copy number cannot discriminate protein outcome, a gene-intrinsic prior over
transmission adds real information, and gating it on surface localisation and on four
orthogonal liability axes yields a shortlist that is defensible on its face. Three
surface leads survive every axis; two strong transmitters are held back on named,
specific liabilities.

The bounds are as important as the claim, and stated plainly. First, the prior explains
about a third of per-gene transmission variance, so it ranks the genome well but
mispredicts individual genes, AP2M1 being the kept-in example; the leads are
prediction-level and transcript-level calls, not surface-protein measurements, and the
next step for any lead is direct protein-level and selectivity measurement. Second, the
network result is deliberately deflationary: interaction topology does not add a
mechanism beyond gene coverage, and the analysis does not pursue network biology as an
explanation. Third, the co-prediction structure within 3q reflects complex-level
co-regulation among 3q genes broadly rather than a tight module among any small named
set; the co-prediction score measures concordance of model outputs, which need not
imply functional co-regulation, and genes may share feature profiles without being
mechanistically linked. Where individual 3q genes such as ACTL6A have been discussed as
exemplars of chromatin-related co-amplification, that chromatin work is parallel and
collaborative to this thesis and is not claimed as a result of it; within the present
analysis, only the transmission and co-prediction structure is load-bearing, and only
TBL1XR1 among the frequently named 3q genes behaves as a consistent co-prediction hub,
coupling with replication and chromatin machinery rather than with a specific 3q trio.

The framework is deliberately leakage-safe and reproducible: the outcome is never a
predictor, splits are fixed and grouped by case, neighbourhood features are recomputed
per fold from training cases only and verified to have zero train-test overlap, metrics
are deduplicated, and the random baseline is reported beside every score. Within those
constraints, the contribution is an attenuation-aware account of which amplifications
reach the protein, and a principled, honestly bounded way to turn that account into
antibody-drug conjugate target hypotheses where the amplification signal alone cannot.
