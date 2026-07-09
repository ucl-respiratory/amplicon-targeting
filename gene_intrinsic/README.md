# Gene-intrinsic copy-number-to-protein transmissibility: export folder

This folder is a clean, self-contained export of the gene-intrinsic strand of the
CN_targeting project: a leakage-safe model of how somatic copy-number amplification
transmits to protein, evidence that the transmission it detects is a stable
gene-intrinsic property, and the downstream application to antibody-drug conjugate
(ADC) target selection in 3q-amplified non-small-cell lung cancer. It contains
scripts, outputs, figures, the manuscript and the provenance record only. No raw
downloaded datasets and nothing genome-scale are included; the large external source
files (CCLE proteomics and copy number, Human Protein Atlas, GTEx) are referenced by
origin and hash in the manifests but are not shipped here.

Everything a claim in the manuscript depends on is present in this folder, and every
number is provenance-verified against its source artefact. The verification is not an
assertion: see the provenance register, the verified manifest, and the two
number-verification tables described in Section 4 below, and the file-level
`SHA256SUMS.txt` at the root. The exported figures and tables have been checked to be
byte-identical (SHA256) to the values recorded in the provenance register.

## 1. What the work is

Somatic amplification is realised only if the extra gene dosage reaches the protein,
and that transmission is attenuated gene by gene. The model predicts, per case and
per gene, whether an amplified gene is over-expressed at the protein level. Its
distinctive value is resolution within the highly amplified regime, where copy number
alone is uninformative, not average discrimination across all genes.

On the 3q cohort the model reaches an area under the precision-recall curve of 0.7515
against a base rate of 0.4191 (lift 1.79); on the pan-cancer cohort, 0.6465 against
0.3983 (lift 1.62). Four independent lines of evidence then show the transmission is
gene-intrinsic: it transfers across tumour lineages (mean pairwise rank concordance
0.840, Kendall W 0.867); it is predictable from external gene properties alone, with
no protein-derived feature (leave-gene-out rank correlation 0.547); it is not
positional (holding out whole chromosome arms changes that correlation by 0.001); and
its mechanism reconciles with post-transcriptional buffering, dominated by
dosage-sensitivity and protein biophysics. Adding a protein-protein interaction
network improves prediction only by widening gene coverage, not because interaction
topology is intrinsically informative.

## 2. Two definitions that fix the outcome

- **Amplified (R4).** A gene in a case is amplified when its ploidy-adjusted copy
  number exceeds one integer step above rounded ploidy AND its integer copy number is
  at least three: `cn >= round(ploidy) + 1` AND `cn >= 3`. This is the single
  amplified definition used for the model anchor, the label, and every metric.

- **Outcome (binary, per (case, gene)).** The target is a binary protein
  over-expression call for each unique (caseid, gene) pair on the R4-filtered cohort,
  taken from the canonical pipeline. All metrics are deduplicated to unique
  (caseid, gene) pairs before computation, because fold stacking duplicates rows.
  On the pan-cancer (ALL) cohort this basis is a row-level base rate of 0.3983 over
  1,758,884 rows, reducing to 700,211 unique pairs. No protein-derived quantity is
  ever a predictor; the outcome is the only protein-derived object in the model.

## 3. The honest bounds

These bounds are stated plainly and hold throughout.

- **The predictor is a prior, not an oracle.** Predicting transmissibility from gene
  properties alone gives a leave-gene-out R-squared of 0.322 (Spearman 0.547, 10-fold,
  seed 2, 600 trees). Roughly two thirds of per-gene transmission variance is left
  unexplained by gene properties. The framework is therefore a gene-intrinsic prior
  over the genome, not a per-gene oracle. A deliberately kept-in error, AP2M1, marks
  this bound in the manuscript.

- **The ADC leads reduced after the surface check.** The upstream pipeline produced
  three leads that survived every liability axis: ATP13A3, ATP1B3 and PLSCR1. "Clean"
  there meant carrying no essentiality, buffering or lineage-risk flag; it never meant
  surface-accessible, because membrane topology was never tested upstream. A separate
  two-part surface check (Human Protein Atlas localisation plus UniProt extracellular
  topology) then reduced this set. ATP13A3 was removed on ectodomain topology: it is a
  10-pass endosomal P5B-type ATPase with zero extracellular residues, so there is no
  accessible cell-surface ectodomain for an antibody to bind, despite ATP13A3 being the
  strongest and most faithfully reproduced transmitter. That leaves ATP1B3 (Na/K-ATPase
  beta-3, the surface antigen CD298; single-pass, 223-amino-acid ectodomain) as the one
  unambiguous ADC-viable lead, and PLSCR1 as marginal and provisional (a 13-amino-acid
  extracellular tail, carried pending epitope confirmation). Two high-transmission
  genes are held back for specific liabilities: TFRC (essential, DepMap dependency
  -0.94) and ITGB5 (lineage-unstable transfer).

- **This is in-silico validation, not experimental proof.** Every result here is
  computed from public multi-omics and annotation. The out-of-distribution CCLE check
  strengthens the gene-intrinsic claim and the surface filter sharpens the lead set,
  but surface-protein density on tumour cells, ectodomain binding, and any therapeutic
  index remain bench questions. The leads are hypotheses, not validated targets. The
  therapeutic-window axis does not rescue any lead into a comfortable position: the one
  surface-viable lead (ATP1B3 / CD298) has the narrowest normal-tissue window, and the
  TCGA lung survival scan is exploratory, underpowered at three leads, and null after
  multiple-testing correction; it is orthogonal to ADC target validity and must not be
  read as evidence against the leads.

## 4. Folder layout

```
export_gene_intrinsic/
  README.md                     this file
  SHA256SUMS.txt                SHA256 of every file in this export

  manuscript/
    manuscript_gene_intrinsic.docx    publication-ready Word
    manuscript_gene_intrinsic.pdf     publication-ready PDF (figures + tables inline)
    manuscript_gene_intrinsic.md      Markdown source

  code/
    transmissibility_predictor.py     the gene-intrinsic predictor (self-contained;
                                       pip deps only: numpy, pandas, xgboost, scipy,
                                       scikit-learn). Header documents the label, the
                                       8 feature groups / 47 columns, and the
                                       leakage audit.
    feature_group_head_to_head_skill/ the leakage-safe head-to-head harness (the
      SKILL.md                        SKILL and its kernel). run_head_to_head encodes
      kernel.py                       four guarantees and refuses to report if any is
                                       skipped: bit-exact anchor gate, label-shuffle
                                       permutation null, grouped-by-case fixed CV,
                                       and dedup-to-unique-(case,gene) metrics with
                                       Random-AP reported beside AU-PR.

  figures_tables/
    consolidated_figure_table_provenance.csv   the provenance register: 44 rows, one
                                       per consolidated figure/table, each with source
                                       path, source SHA256, Science version id,
                                       consolidated filename and consolidated SHA256.
    tournament_master/                the consolidated figure and table set (44 files),
                                       named in the tournament_master convention and
                                       indexed 1:1 by the register. Figures 4-8 and 18
                                       and Tables 5-16 are the gene-intrinsic core;
                                       the remainder (purity backbone, co-prediction,
                                       network, hardened ADC) are the surrounding
                                       thesis context the manuscript cites.

  provenance/
    science_phase_manifest.csv        the verified manifest: 59 deliverables, each with
                                       local path, resolves_on_disk, SHA256 and Science
                                       version id.
    science_phase_verification.csv    per-claim source verification.
    science_phase_persistence_gate.csv  end-of-phase persistence gate.
    manuscript_number_verification.csv  40/40 manuscript numbers recomputed from source
                                       CSVs (a real source comparison, not a
                                       string-presence check).
    references_verified.csv           the verified reference list (live CrossRef +
                                       PubMed, zero retracted).

  validation/
    validation_summary.md             the independent in-silico validation write-up.
    predictor_ccle_validation.csv/.png  out-of-distribution CCLE check of the predictor.
    predictor_ccle_validation_summary.csv
    leads_validation.csv              the three leads: CCLE transmission, HPA
                                       localisation, UniProt topology, surface verdict.
    therapeutic_window.csv/.png       GTEx normal-tissue burden per lead.
    clinical_relevance.csv            TCGA lung survival scan (exploratory, null).
    validation_number_verification.csv
    validation_persistence_gate.csv
```

## 5. Leakage-safe design

Leakage safety is the project's highest priority and is enforced, not assumed. No
protein-derived feature is ever a predictor; splits are grouped by case; the fixed CV
splits are used unchanged; any patient-level value propagated over a network or
neighbourhood is computed from training-fold cases only and recomputed per fold, so a
validation case never informs its own features. The predictor's features come from
external sources only (gnomAD, DepMap, CORUM, UniProt, Ensembl, GTEx, STRING) and the
label from CPTAC proteomics; the two never share a source, so the model cannot be
circular. A data-driven audit against the feature manifest confirms that 0 of 16
feature sources reference CPTAC tumour proteomics or the label, and no feature column
is byte-identical to the label. The one borderline group (G8, GTEx normal-tissue RNA)
is reported with and without.

## 6. Reproducibility

Deterministic throughout: seed 2 everywhere, XGBoost 600 trees, single-thread where it
affects reproducibility, metrics deduplicated to unique (caseid, gene) with Random-AP
(base rate) reported beside every AU-PR. The predictor header records the 8 feature
groups and 47 encoded columns; G4 (protein turnover / half-life) is left unavailable
rather than fabricated, because no proteome-wide half-life atlas was reachable. Package
versions, data hashes, and external-resource names and hashes are recorded in the
manifests under `provenance/` and in the validation manifest referenced there.

British English; no em dashes; methods first, claims second.
