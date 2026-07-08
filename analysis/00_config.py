# =============================================================================
# 00_config.py  --  Configuration for the CN-targeting ANALYSIS pipeline
# -----------------------------------------------------------------------------
# The analysis pipeline assumes `data_download/` has already been run. It reads
# that pipeline's OUT/ tables (+ a few RAW/ auxiliary files) and reproduces every
# figure and in-text value of the manuscript. Nothing here downloads data.
#
# INPUT LOCATION
#   Set CNT_DATA to the data_download output root (the dir that contains
#   tables/, and whose sibling data/raw + data/depmap hold auxiliary files):
#       export CNT_DATA=/path/to/data_download/from_source
#   Individual paths can be overridden with CNT_<NAME> env vars (see PATHS).
#
# DETERMINISM
#   All stochastic steps use the seeds pinned here. Given identical input tables
#   the pipeline is bit-for-bit reproducible; a *fresh* data_download run drifts
#   within the tolerances documented in data_download/DOWNLOAD_PIPELINE.md.
# =============================================================================
import os
from pathlib import Path

# ---- Seeds (pinned; do not change without re-verifying values) --------------
RANDOM_SEED       = 42     # global default (numpy / sklearn splits)
SEED_JAMES_SPLIT  = 42     # xgboost train_test_split
SEED_BOOTSTRAP    = 0      # same-cell donor-block bootstrap + permutation
SEED_PATIENT_CV   = 2022   # patient-level model splits (James' convention)

# ---- Root resolution --------------------------------------------------------
def _root() -> Path:
    env = os.environ.get("CNT_DATA")
    if env:
        return Path(env).expanduser().resolve()
    # fall back to a sibling data_download/from_source run
    here = Path(__file__).resolve().parent
    for cand in (here.parent / "data_download" / "from_source",
                 here.parent / "data_download" / "from_source_run"):
        if cand.exists():
            return cand
    # last resort: current dir
    return Path(os.environ.get("CNT_DATA", here.parent / "data_download" / "from_source")).resolve()

DATA_ROOT = _root()
TAB   = Path(os.environ.get("CNT_TABLES", DATA_ROOT / "tables"))
RAW   = Path(os.environ.get("CNT_RAW",    DATA_ROOT / "data" / "raw"))
DEPMAP= Path(os.environ.get("CNT_DEPMAP", DATA_ROOT / "data" / "raw" / "depmap"))
PARSED= Path(os.environ.get("CNT_PARSED", DATA_ROOT / "data" / "parsed"))

# ---- Output dirs (under analysis/out by default) ----------------------------
OUT      = Path(os.environ.get("CNT_ANALYSIS_OUT",
                               Path(__file__).resolve().parent / "out"))
DIR_FIG  = OUT / "figures"
DIR_TAB  = OUT / "tables"
DIR_REP  = OUT / "reports"
for d in (DIR_FIG, DIR_TAB, DIR_REP):
    d.mkdir(parents=True, exist_ok=True)

# ---- Input file paths (the CONTRACT with data_download) ---------------------
PATHS = {
    # core modelling feature matrix (stage 13)
    "str_omic":      TAB / "str_omic_table_rebuilt.csv",
    # per-layer / annotation parquet (stages 11-13)
    "protein_core":  TAB / "omic_table_protein_core.parquet",
    "describeprot":  TAB / "describeprot_gene_features.parquet",
    # chromatin accessibility (stage 14)
    "atac":          TAB / "atac_gene_promoter_accessibility.parquet",
    # membrane topology (stage 15)
    "topology":      TAB / "surface_topology_uniprot.parquet",
    # Ensembl gene coordinates (stage 17) -> gene length / 1Mb density features
    "ensembl_coords": TAB / "ensembl_gene_coords.parquet",
    # GTEx v8 bulk median TPM, gene x tissue (stage 18) -> off-target safety
    "gtex_bulk":     TAB / "gtex_v8_median_tpm.parquet",
    # CORUM complexes (stage 05)
    "corum":         PARSED / "corum_humanComplexes.txt",
    # DepMap 23Q4 (stage 05)
    "dep_prob":      DEPMAP / "CRISPRGeneDependency.csv",
    "common_ess":    DEPMAP / "AchillesCommonEssentialControls.csv",
    "nonessential":  DEPMAP / "AchillesNonessentialControls.csv",
    "depmap_cn":     DEPMAP / "Copy_Number_Public_23Q4.csv",
    "depmap_expr":   DEPMAP / "Expression_Public_23Q4.csv",
    "depmap_model":  DEPMAP / "Model.csv",
    "depmap_prot":   DEPMAP / "Proteomics.csv",
    # HPA (stage 04)
    "hpa_normal":    RAW / "normal_tissue.tsv",
    "hpa_singlecell":RAW / "rna_single_cell_type.tsv",
    # CELLxGENE malignant-cell slices (stage 16)
    "cxg_luad":      TAB / "cellxgene_LUAD_malignant.parquet",
    "cxg_lscc":      TAB / "cellxgene_LSCC_malignant.parquet",
    # ADC atlas nominated genes (stage 04 raw -> adapter)
    "adc_genes":     TAB / "adc_distinct_genes.csv",
    # CPTAC matched tumour/normal proteomics (stage 03 CPTAC ingestion) -> off-target safety
    "cptac_matched_protein": TAB / "cptac_matched_protein.csv",
    "cptac_matched_pheno":   TAB / "cptac_matched_pheno.csv",
}

# ---- Analysis parameters (pinned to the manuscript) -------------------------
TUMOR_CODES        = ["CCRCC", "GBM", "LSCC", "LUAD", "PDA", "UCEC"]
AMP_THRESHOLD      = 1.4    # cn_adjusted >= 1.4 == amplified (>=40% gain over ploidy)
RECURRENCE_FREQ    = 0.20   # amplicon recurrent when amplified in >=20% of a type
RECURRENCE_MIN_N   = 8      # and in >=8 samples
FDR_ALPHA          = 0.10   # BH-FDR threshold for co-elevation nomination
COELEV_MIN_FRAC    = 0.50   # elevated in >=50% of amplified tumours
REL_TISSUE_HI      = 0.80   # prot.rel.tissue > 0.8 == co-elevated (James' target)
DEP_ESSENTIAL      = 0.50   # DepMap dependency >= 0.5 == essential
SC_DETECT_NTPM     = 10     # HPA single-cell detection floor
SC_BINDING_NTPM    = 25     # binding-relevant threshold
N_BOOTSTRAP        = 2000   # donor-block bootstrap resamples (Fig 6 CIs)
N_PERMUTATION      = 1000   # marginal-preserving permutation test

# GTEx tissue -> cancer-type map (mirror of data_download)
GTEX_TISSUE_MAP = {"CCRCC":"Kidney - Cortex","GBM":"Brain - Cortex","LSCC":"Lung",
                   "LUAD":"Lung","PDA":"Pancreas","UCEC":"Uterus"}

def require(*keys):
    """Assert the named input files exist; raise a clear pointer to data_download."""
    missing = [(k, str(PATHS[k])) for k in keys if not PATHS[k].exists()]
    if missing:
        lines = "\n".join(f"    - {k}: {p}" for k, p in missing)
        raise FileNotFoundError(
            "Required data_download outputs are missing:\n" + lines +
            f"\n  CNT_DATA resolved to: {DATA_ROOT}\n"
            "  Run `Rscript data_download/run_all.R` first, or set CNT_DATA to its output root.")

if __name__ == "__main__":
    print("CNT_DATA  =", DATA_ROOT)
    print("OUT       =", OUT)
    present = {k: PATHS[k].exists() for k in PATHS}
    for k, ok in present.items():
        print(f"  [{'x' if ok else ' '}] {k:14s} {PATHS[k]}")
