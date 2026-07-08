# =============================================================================
# cnt_io.py  --  shared loaders + value-recording for the analysis pipeline
# -----------------------------------------------------------------------------
# Every numbered stage imports from here. Loaders are cached (lru_cache) so the
# big feature table is read once per process. `record_values(stage, d)` writes a
# per-stage JSON into OUT/reports/values/, which stage 90 aggregates into the
# paper_values.json manifest and the verification report.
# =============================================================================
import json, functools
from pathlib import Path
import pandas as pd, numpy as np

import importlib.util, sys
_cfg_path = Path(__file__).resolve().parent / "00_config.py"
_spec = importlib.util.spec_from_file_location("cnt_config", _cfg_path)
cfg = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(cfg)

PATHS, TUMOR_CODES = cfg.PATHS, cfg.TUMOR_CODES
DIR_FIG, DIR_TAB, DIR_REP = cfg.DIR_FIG, cfg.DIR_TAB, cfg.DIR_REP
VALDIR = DIR_REP / "values"; VALDIR.mkdir(parents=True, exist_ok=True)

# ---- loaders ----------------------------------------------------------------
@functools.lru_cache(maxsize=1)
def load_str_omic():
    """The wide James-compatible feature matrix (gene x caseid rows)."""
    cfg.require("str_omic")
    return pd.read_csv(PATHS["str_omic"], low_memory=False)

@functools.lru_cache(maxsize=1)
def load_atac():
    cfg.require("atac"); return pd.read_parquet(PATHS["atac"])

@functools.lru_cache(maxsize=1)
def load_topology():
    cfg.require("topology"); return pd.read_parquet(PATHS["topology"])

@functools.lru_cache(maxsize=1)
def load_corum():
    """CORUM human complexes -> dict complex_id -> set(genes)."""
    cfg.require("corum")
    df = pd.read_csv(PATHS["corum"], sep="\t", dtype=str)
    gcol = "subunits_gene_name" if "subunits_gene_name" in df.columns else \
           [c for c in df.columns if "gene_name" in c.lower()][0]
    icol = "complex_id" if "complex_id" in df.columns else df.columns[0]
    out = {}
    for _, r in df.iterrows():
        genes = [g for g in str(r[gcol]).split(";") if g and g != "nan"]
        if genes: out[r[icol]] = set(genes)
    return out

@functools.lru_cache(maxsize=1)
def load_dep_prob():
    """DepMap CRISPRGeneDependency: models x genes (0-1). Strips ' (Entrez)'."""
    cfg.require("dep_prob")
    df = pd.read_csv(PATHS["dep_prob"], index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df

@functools.lru_cache(maxsize=1)
def load_common_essential():
    cfg.require("common_ess")
    df = pd.read_csv(PATHS["common_ess"])
    col = df.columns[0]
    return set(str(g).split(" (")[0] for g in df[col].dropna())

@functools.lru_cache(maxsize=1)
def load_nonessential():
    cfg.require("nonessential")
    df = pd.read_csv(PATHS["nonessential"])
    col = df.columns[0]
    return set(str(g).split(" (")[0] for g in df[col].dropna())

@functools.lru_cache(maxsize=1)
def load_depmap_prot():
    cfg.require("depmap_prot")
    df = pd.read_csv(PATHS["depmap_prot"], index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df.loc[:, ~df.columns.duplicated()].copy()

@functools.lru_cache(maxsize=1)
def load_depmap_model():
    cfg.require("depmap_model"); return pd.read_csv(PATHS["depmap_model"])

@functools.lru_cache(maxsize=1)
def load_depmap_cn():
    cfg.require("depmap_cn")
    df = pd.read_csv(PATHS["depmap_cn"], index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]; return df

def load_hpa_singlecell():
    cfg.require("hpa_singlecell")
    return pd.read_csv(PATHS["hpa_singlecell"], sep="\t")

def load_hpa_normal():
    cfg.require("hpa_normal")
    return pd.read_csv(PATHS["hpa_normal"], sep="\t")

def load_cellxgene(code):
    key = {"LUAD":"cxg_luad","LSCC":"cxg_lscc"}[code]
    cfg.require(key); return pd.read_parquet(PATHS[key])

# ---- value recording --------------------------------------------------------
def _json_safe(x):
    if isinstance(x, (np.integer,)):  return int(x)
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.ndarray,)):  return x.tolist()
    if isinstance(x, dict):  return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [_json_safe(v) for v in x]
    return x

def record_values(stage: str, values: dict):
    """Write OUT/reports/values/<stage>.json (one per stage)."""
    p = VALDIR / f"{stage}.json"
    json.dump(_json_safe(values), open(p, "w"), indent=1)
    return p

def savefig(fig, name):
    """Save a figure to OUT/figures/<name> at 150 dpi (tight)."""
    out = DIR_FIG / name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    return out


def load_cellxgene(code):
    """Load a data_download stage-16 malignant-cell detection slice for a cancer
    type ('LUAD' or 'LSCC'). Columns: donor_id, nnz, then one 0/1 column per gene."""
    import pandas as pd
    key = {"LUAD": "cxg_luad", "LSCC": "cxg_lscc"}[code]
    return pd.read_parquet(PATHS[key])
