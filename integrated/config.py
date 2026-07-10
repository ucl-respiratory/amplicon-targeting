# =============================================================================
# integrated/config.py  --  Configuration for the INTEGRATED paper pipeline
# -----------------------------------------------------------------------------
# This pipeline builds "From amplicon to antigen": a single, integrated account
# that (1) quantifies where copy-number dosage is transmitted to protein and
# where it is gated, (2) predicts transmissibility from gene properties with no
# protein data, (3) nominates surface ADC co-target sets from BOTH measured and
# predicted transmissibility (tagging which is which), and (4) demonstrates
# AND-gate selectivity in single-cell data.
#
# It reads the SAME data_download/ outputs the analysis/ pipeline consumes
# (one data contract, no duplication of the download) and writes to its own
# integrated/{figures,tables,reports}. It does not modify analysis/ or
# gene_intrinsic/.
#
# INPUT LOCATION
#   export CNT_DATA=/path/to/hackathon/data_download/from_source
#   (defaults to the sibling data_download/from_source run if unset)
#
# DETERMINISM
#   All stochastic steps use SEED below. Given identical data_download tables
#   the pipeline is reproducible; a fresh data_download run drifts within the
#   tolerances in data_download/DOWNLOAD_PIPELINE.md.
# =============================================================================
import os
from pathlib import Path

# ---- Seed (single global; matches the two source pipelines' conventions) ----
SEED              = 2       # gene-intrinsic predictor + empirical-Bayes resamples
SEED_JAMES_SPLIT  = 42      # xgboost train_test_split (matches analysis/)
SEED_PATIENT_CV   = 2022    # patient-level CV splits (matches analysis/)

# ---- Root resolution (identical contract to analysis/00_config.py) ----------
def _root() -> Path:
    env = os.environ.get("CNT_DATA")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve().parent
    for cand in (here.parent / "data_download" / "from_source",
                 here.parent / "data_download" / "from_source_run"):
        if cand.exists():
            return cand
    return (here.parent / "data_download" / "from_source").resolve()

DATA_ROOT = _root()
# The data_download R pipeline writes its assembled tables to <root>/out/tables
# (DIR_TAB = OUT/tables in 00_config.R). Prefer that; fall back to <root>/tables
# for older layouts.
def _tab_default():
    for cand in (DATA_ROOT / "out" / "tables", DATA_ROOT / "tables"):
        if (cand / "str_omic_table_rebuilt.csv").exists():
            return cand
    return DATA_ROOT / "out" / "tables"
TAB    = Path(os.environ.get("CNT_TABLES", _tab_default()))

# CELLxGENE malignant single-cell slices are written by the census extraction
# (data_download stage 16) to <root>/tables, independent of the pipeline's
# out/tables. Resolve them wherever they landed.
def _tab_census():
    for cand in (DATA_ROOT / "tables", DATA_ROOT / "out" / "tables", TAB):
        if (cand / "cellxgene_LUAD_malignant.parquet").exists():
            return cand
    return DATA_ROOT / "tables"
_TAB_CENSUS = _tab_census()
RAW    = Path(os.environ.get("CNT_RAW",    DATA_ROOT / "data" / "raw"))
DEPMAP = Path(os.environ.get("CNT_DEPMAP", DATA_ROOT / "data" / "raw" / "depmap"))
PARSED = Path(os.environ.get("CNT_PARSED", DATA_ROOT / "data" / "parsed"))

# ---- Output dirs (under integrated/ by default) -----------------------------
OUT     = Path(os.environ.get("CNT_INTEGRATED_OUT",
                              Path(__file__).resolve().parent))
DIR_FIG = OUT / "figures"
DIR_TAB = OUT / "tables"
DIR_REP = OUT / "reports"
for d in (DIR_FIG, DIR_TAB, DIR_REP):
    d.mkdir(parents=True, exist_ok=True)

# ---- Input file paths (the CONTRACT with data_download) ---------------------
# Same keys/paths as analysis/00_config.py so the two pipelines stay in lockstep.
PATHS = {
    "str_omic":       TAB / "str_omic_table_rebuilt.csv",
    "protein_core":   TAB / "omic_table_protein_core.parquet",
    "annotated":      TAB / "omic_table_annotated.parquet",
    "describeprot":   TAB / "describeprot_gene_features.parquet",
    "atac":           PARSED / "atac_gene_promoter_accessibility.parquet",
    "topology":       TAB / "surface_topology_uniprot.parquet",
    "ensembl_coords": TAB / "ensembl_gene_coords.parquet",
    "gtex_bulk":      TAB / "gtex_v8_median_tpm.parquet",
    "corum":          PARSED / "corum_humanComplexes.txt",
    "dep_prob":       DEPMAP / "CRISPRGeneDependency.csv",
    "common_ess":     DEPMAP / "AchillesCommonEssentialControls.csv",
    "nonessential":   DEPMAP / "AchillesNonessentialControls.csv",
    "cxg_luad":       _TAB_CENSUS / "cellxgene_LUAD_malignant.parquet",
    "cxg_lscc":       _TAB_CENSUS / "cellxgene_LSCC_malignant.parquet",
    "cxg_gbm":        _TAB_CENSUS / "cellxgene_GBM_malignant.parquet",
    "hpa_normal":     RAW / "normal_tissue.tsv",
    "hpa_singlecell": RAW / "rna_single_cell_type.tsv",
}

# ---- Amplification / co-elevation parameters (CN>=1.4 primary, per decision)--
TUMOR_CODES     = ["CCRCC", "GBM", "LSCC", "LUAD", "PDA", "UCEC"]
AMP_THRESHOLD   = 1.4    # cn_adjusted >= 1.4 == amplified (ploidy-adjusted, V9 basis)
RECURRENCE_FREQ = 0.20   # amplicon recurrent when amplified in >=20% of a type
RECURRENCE_MIN_N= 8
FDR_ALPHA       = 0.10   # BH-FDR for co-elevation nomination
COELEV_MIN_FRAC = 0.50   # elevated in >=50% of amplified tumours
REL_TISSUE_HI   = 0.80   # prot.rel.tissue > 0.8 == co-elevated
DEP_ESSENTIAL   = 0.50   # DepMap dependency >= 0.5 == essential
SC_DETECT_NTPM  = 10     # single-cell detection floor
SC_BINDING_NTPM = 25     # binding-relevant threshold
SC_DEPTH_BINS   = 10     # per-cell depth (nnz) deciles for the depth-stratified co-detection null
N_BOOTSTRAP     = 2000
N_PERMUTATION   = 1000

# ---- Predictor: gene-property feature groups (NO protein-derived feature) ----
# Mirrors the gene_intrinsic transmissibility predictor. The outcome is the
# protein-derived transmissibility; predictors are gene properties only.
PREDICTOR_GROUPS = {
    "dosage":     ["gnomad_LOEUF", "gnomad_pLI", "gnomad_mis_z",
                   "dep_mean_effect", "dep_frac_dependent"],
    "complex":    ["in_complex", "n_complexes", "complex_size", "has_complex"],
    "biophysics": ["length", "mol_weight", "isoelectric_point", "gravy",
                   "aggregation_propensity", "tm_domain_count", "signal_peptide",
                   "vsl2_disorder", "psipred_helix", "psipred_strand",
                   "psipred_coil", "asaquick_buried"],
    "mrna":       ["transcript_length", "gc_content", "n_isoforms",
                   "utr5_length", "utr3_length", "codon_optimality"],
    "evolution":  ["dn_ds", "phylop_mean", "gene_age_proxy"],
    "function":   ["is_tf", "is_kinase", "is_receptor", "is_enzyme"],
    "breadth":    ["n_tissues_expressed", "tau"],
    "network":    ["degree", "weighted_degree", "betweenness"],
}

# ---- Empirical-Bayes combine parameters -------------------------------------
EB_REFERENCE_MIN_CASES = 150   # genes with >= this many amplified cases = reference truth
EB_RESAMPLES           = 40    # resamples per cohort-size grid point
EB_COHORT_GRID         = [0, 5, 8, 12, 18, 25, 35, 50, 75, 110, 160, 232]
EB_HOLDOUT_LINEAGE     = "LSCC"  # lineage held out for the extrapolation demo

# ---- Confidence tiers for nominated surface targets -------------------------
# Every nominated antigen is tagged by the evidence supporting it.
CONF_TIERS = {
    "measured_high":  "co-elevation measured (FDR<0.1) AND surface-accessible",
    "measured_pred":  "measured + prior-corroborated (both high)",
    "predicted_only": "prior-nominated (no proteome in this context); ranked, lower confidence",
}

# ---- Surface-topology gate --------------------------------------------------
MIN_ECTODOMAIN_AA = 20   # UniProt extracellular residues required for an accessible epitope

# ---- Reused pre-computed tables (gene_intrinsic export; committed to repo) ---
# The transmissibility atlas (measured observed + predicted OOF prior) and the
# gene feature table are committed under gene_intrinsic/. Until a fresh
# data_download+predictor re-fit lands, these provide the predictor outputs.
GI_EXPORT = Path(__file__).resolve().parent.parent / "gene_intrinsic" / "figures_tables" / "tournament_master"
GI_PATHS = {
    "atlas":         GI_EXPORT / "tournament_master_tab12_transmissibility_atlas.csv",
    "feature_table": GI_EXPORT / "tournament_master_tab09_gene_feature_table.csv",
    "transfer":      GI_EXPORT / "tournament_master_tab06_transfer_stability.csv",
    "shap":          GI_EXPORT / "tournament_master_tab16_transmissibility_shap.csv",
    "leads_valid":   Path(__file__).resolve().parent.parent / "gene_intrinsic" / "validation" / "leads_validation.csv",
}

def require(*keys):
    """Assert named data_download inputs exist; point to data_download if not."""
    missing = [(k, str(PATHS[k])) for k in keys if not PATHS[k].exists()]
    if missing:
        lines = "\n".join(f"    - {k}: {p}" for k, p in missing)
        raise FileNotFoundError(
            "Required data_download outputs are missing:\n" + lines +
            f"\n  CNT_DATA resolved to: {DATA_ROOT}\n"
            "  Run: Rscript data_download/run_all.R\n")

def have(*keys) -> bool:
    """True if all named inputs exist (for runnable-now vs needs-data gating)."""
    return all(PATHS[k].exists() for k in keys)
