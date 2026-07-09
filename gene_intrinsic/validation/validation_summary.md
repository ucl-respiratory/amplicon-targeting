# Independent in-silico validation of the transmissibility predictor and the 3q ADC leads

## 0. What this validation does and does not do

This is an in-silico validation against data the transmissibility model never saw:
the CCLE cell-line panel (copy number and quantitative proteomics), the Human
Protein Atlas, GTEx, and TCGA lung. It raises or lowers confidence in the
existing results and, where the evidence demanded, it changed the lead set. It is
not experimental proof. Every claim below is a transcript-and-annotation-level
call from public data, not a bench measurement.

The headline: the transmission model transfers out of distribution, weakly but
genuinely, and a stricter surface-accessibility filter removed a candidate that
should not have survived. That is the validation behaving as designed, not the
lead set collapsing.

## 1. Stage 1: the predictor validates out of distribution on CCLE

The CPTAC-trained gene-level transmissibility atlas was tested against an
independent observed transmission label built in the CCLE cell-line panel. The
two sample sets are disjoint (CPTAC: 560 primary-tumour patient cases; CCLE:
1,804 immortalised cell lines; zero identifier or name overlap), so this is a
genuine out-of-distribution test on the 377 lines with matched copy number and
proteomics.

CCLE has no matched-normal reference, so a CPTAC-identical (tumour-versus-normal)
label is impossible. The closest analogue is a between-line contrast: for each
gene, the fraction of amplified lines whose protein exceeds the panel's own
unamplified lines for that gene. This measures the same underlying biology
(does amplification raise protein) against a different baseline, so a weak
correlation is interpreted as a label-construction gap, not silently as
predictor failure.

Result: a weak but genuine positive correlation, agreed by both label variants.
Between-line Spearman +0.112 (n = 3,780 genes at at least 10 amplified lines),
z-score variant +0.122; both stable at the stricter at-least-20 threshold
(+0.111 and +0.122, n = 1,952). Both far exceed a label-shuffle null
(permutation p = 0.001; largest shuffled absolute rho about 0.05 against an
observed 0.11 to 0.12). The constitutive-expression confound, the specific
failure mode of a between-line comparator, was tested and cleared: controlling
for mean panel expression retains 84 to 87 per cent of the signal (partial rho
+0.094 to +0.107), and the correlation is stronger among low-expression genes
(+0.150) than high-expression genes (+0.054), exactly as the comparator's
geometry predicts.

Reconciliation with the within-CPTAC leave-gene-out figure of 0.547: the drop to
0.112 is an expected magnitude change, not predictor fragility. The out-of-
distribution test crosses two gaps at once, sample type (primary tumour to cell
line) and label construction (tumour-versus-normal to between-line), and the
confound analysis shows the residual signal is real and concentrated where the
between-line contrast is actually informative. The predictor is a gene-intrinsic
prior that partially generalises, not a cross-context oracle.

## 2. Stage 2: transmission holds for all three leads; the surface check removes ATP13A3

On observed CCLE transmission, all three leads corroborate their CPTAC
prediction, on a reasonable number of amplified lines:

- ATP13A3: 0.694 on 49 amplified lines (CPTAC observed 0.684, a close match).
- ATP1B3: 0.774 on 31 amplified lines.
- PLSCR1: 0.733 on 30 amplified lines.

The two-part surface check, run because localisation annotation alone is not
sufficient for an antibody-drug conjugate target, separates them decisively:

- ATP13A3 fails. UniProt topology shows a 10-pass endosomal P5B-type ATPase
  (recycling, early and late endosome membrane) with zero extracellular
  topological domains; every non-cytoplasmic loop is lumenal, facing the
  endosome interior. HPA places its main signal at the nucleoli. There is no
  accessible cell-surface ectodomain for an antibody to bind. This independently
  confirms the parallel-work exclusion, and it does so despite ATP13A3 being the
  strongest and most faithfully reproduced transmitter. ATP13A3 is removed as an
  ADC target on surface-accessibility grounds.
- ATP1B3 passes. Single-pass, a 223-amino-acid extracellular domain,
  plasma-membrane localised (HPA Supported); this is Na/K-ATPase beta-3, the
  surface antigen CD298. It is the one unambiguous ADC-viable lead.
- PLSCR1 is marginal. Single-pass, but only a 13-amino-acid extracellular tail
  with the bulk of the protein (288 residues) cytoplasmic and a substantial
  nuclear pool. A 13-residue tail is not a confident ADC epitope; it is carried
  as provisional pending epitope confirmation.

This exposes a real limit of the upstream pipeline: "clean" there meant carrying
no essentiality, buffering or lineage-risk flag; it never meant surface-
accessible. Membrane topology was never tested, so a clean lead could still fail
as an ADC target for want of an ectodomain.

## 3. Stage 3: therapeutic window, with ATP1B3's normal-tissue burden the binding constraint

None of the three is tumour-enriched; HPA calls all three "low tissue
specificity, detected in all" and "low cancer specificity". Normal-tissue RNA
burden (GTEx, 36 tissues) with the surface verdict attached to each:

- ATP1B3, the one surface-viable lead, has the narrowest window: median 81.7
  nTPM across all 36 tissues, and higher in normal lung (183.5) than in most
  other tissues, with substantial heart (97.4) and kidney (103.0). As CD298 this
  is a known, broadly expressed surface antigen, so off-tumour toxicity is the
  binding constraint for exactly the lead that survived the surface filter. The
  ectodomain that makes it targetable is present on vital normal tissue too.
- PLSCR1: median 25.1 nTPM, expressed in all 36 tissues, a narrow window;
  provisional anyway on the epitope question.
- ATP13A3: median 2.7 nTPM, the widest window on paper, but moot as a surface
  target given Stage 2.

The therapeutic-window axis does not rescue any lead into a comfortable
position: the surface-viable lead has the worst window, and the widest-window
gene is not surface-accessible.

## 4. Stage 4: survival is null and exploratory, and orthogonal to target validity

In TCGA lung (LUAD 494 patients, 177 events; LUSC 477 patients, 203 events), no
lead shows an overall-survival association surviving multiple-testing
correction. The strongest expression signal is PLSCR1 in LUAD (Cox HR 1.14 per
standard deviation, raw p = 0.028, FDR 0.17); every amplification log-rank test
is FDR above 0.38, and the nominal LUSC amplification trends point, if anything,
towards better survival with amplification. HPA independently calls all three
unprognostic.

This is exploratory and orthogonal. A null survival scan neither supports nor
refutes an ADC target: antibody-drug conjugate value rests on transmission,
surface accessibility and the therapeutic window, not on the gene's prognostic
association. It is reported for completeness and must not be read as evidence
against the leads.

## 5. What is now independently supported, and what still needs the bench

Independently supported by this validation:

- The transmissibility model transfers out of distribution. An entirely
  independent cell-line panel, under a deliberately different label, recovers a
  real gene-intrinsic signal (Stage 1). This strengthens the central claim that
  transmission propensity is a gene-intrinsic property, beyond the CPTAC cohorts.
- ATP1B3 survives every axis testable in silico: it transmits (Stage 2), it is
  plasma-membrane localised with a substantial 223-amino-acid ectodomain
  (Stage 2), and it is robustly expressed in lung tumours (Stage 3). It is the
  one lead that is coherent as an ADC target on paper.

Still requires the bench, and cannot be settled by any public dataset used here:

- Surface-protein density on tumour cells, not RNA or annotation: an ADC needs
  sufficient cell-surface antigen copies per cell, which only protein-level
  surface quantification (for example flow cytometry or quantitative surface
  proteomics) can establish.
- Epitope accessibility, especially for PLSCR1: whether the 13-residue
  extracellular tail is real, exposed and bindable in the native membrane.
- Tumour-versus-normal selectivity at the protein level: the ATP1B3 normal-tissue
  burden (Stage 3) is an RNA-level warning that must be resolved with protein
  immunohistochemistry across vital tissues before ATP1B3 could be considered
  safe.
- Any therapeutic window in vivo: tolerability and the actual off-tumour
  toxicity of a CD298-directed conjugate.

## 6. The honest bottom line

The transmission model validated out of distribution, and the ADC application
behaved as designed when a stricter accessibility filter removed a candidate
(ATP13A3) that should not have survived. Of the original three clean leads, one
(ATP1B3) survives every in-silico axis but carries a serious normal-tissue
liability, one (PLSCR1) is marginal on epitope grounds, and one (ATP13A3) is out
as a surface target despite being the best transmitter. This validation raised
confidence in the model, lowered it for two of three leads on the surface and
window axes, and changed the lead set where the evidence demanded. That is the
purpose of an independent validation; it is not a substitute for the bench.
