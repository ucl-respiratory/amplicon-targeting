# =============================================================================
# _depmap.py  --  DepMap 23Q4 essentiality + cross-cohort validation (shared)
# -----------------------------------------------------------------------------
# gene_dependency_summary(): per-gene mean CRISPR dependency, fraction of lines
#   dependent (>0.5), and common-essential / nonessential control flags. Built
#   from CRISPRGeneDependency.csv + Achilles{Common,Non}essentialControls.csv.
# depmap_validation(): cross-cohort/cross-platform check that amplicon-high
#   cell lines of the matched lineage over-express each nominated surface target
#   (one-sided Mann-Whitney on CCLE proteomics), BH-FDR corrected.
# Used by stages 21-24 and Fig 3. Paper-exact.
# =============================================================================
import functools
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_io, _amplicon
from cnt_io import cfg

LINEAGE_MAP = {"LUAD":"Lung","LSCC":"Lung","CCRCC":"Kidney","UCEC":"Uterus",
               "PDA":"Pancreas","GBM":"CNS/Brain"}

def _bh_fdr(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    ranked = np.empty(n); cummin = 1.0
    for i in range(n-1, -1, -1):
        idx = order[i]; val = p[idx]*n/(i+1); cummin = min(cummin, val); ranked[idx] = cummin
    return np.clip(ranked, 0, 1)

@functools.lru_cache(maxsize=1)
def gene_dependency_summary():
    dep = cnt_io.load_dep_prob()                 # models x genes, cols de-Entrez'd
    common_essential = cnt_io.load_common_essential()
    nonessential = cnt_io.load_nonessential()
    gd = pd.DataFrame({"mean_dependency": dep.mean(axis=0),
                       "frac_lines_dependent": (dep > cfg.DEP_ESSENTIAL).mean(axis=0),
                       "n_lines": dep.notna().sum(axis=0)})
    gd.index.name = "gene"; gd = gd.reset_index()
    gd["common_essential"] = gd.gene.isin(common_essential).astype(int)
    gd["nonessential"] = gd.gene.isin(nonessential).astype(int)
    return gd

@functools.lru_cache(maxsize=1)
def _good_surf():
    """Reliably co-elevated, surface, lift>0.15, elevated in >=50% amplified."""
    co = _amplicon.amplicon_coelevation()
    gs = co[(co.is_surface == 1) & (co.fdr < 0.1) & (co.lift > 0.15) &
            (co.p_high_amp >= cfg.COELEV_MIN_FRAC)].copy()
    gs["arm"] = gs.band.str.extract(r"^(\d+[pq])")
    return gs

@functools.lru_cache(maxsize=1)
def depmap_validation():
    model = cnt_io.load_depmap_model()
    prot  = cnt_io.load_depmap_prot()
    cn    = cnt_io.load_depmap_cn()
    df    = cnt_io.load_str_omic()

    def amp_score_lines(genes_arm, ids):
        gg = [g for g in genes_arm if g in cn.columns]
        return cn.loc[[i for i in ids if i in cn.index], gg].mean(axis=1)

    def validate_target(gene, tumor_code, arm):
        lin = LINEAGE_MAP.get(tumor_code)
        ids = [m for m in model[model.OncotreeLineage == lin].ModelID if m in prot.index and m in cn.index]
        if len(ids) < 12 or gene not in prot.columns: return None
        arm_genes = df[df.arm == arm].gene.unique()
        sc = amp_score_lines(arm_genes, ids).dropna()
        if len(sc) < 12: return None
        amp = sc[sc >= sc.quantile(0.66)].index
        noamp = sc[sc < sc.quantile(0.5)].index
        hi = prot.loc[[l for l in amp if l in prot.index], gene].dropna()
        lo = prot.loc[[l for l in noamp if l in prot.index], gene].dropna()
        if len(hi) < 4 or len(lo) < 4: return None
        u, p = mannwhitneyu(hi, lo, alternative="greater")
        return {"gene":gene,"tumor_code":tumor_code,"arm":arm,"lineage":lin,
                "amp_median":float(hi.median()),"noamp_median":float(lo.median()),
                "delta":float(hi.median()-lo.median()),"n_amp":len(hi),"n_noamp":len(lo),"p":float(p)}

    targets = _good_surf()[["gene","tumor_code","arm"]].drop_duplicates()
    val = [validate_target(r.gene, r.tumor_code, r.arm) for r in targets.itertuples()]
    val = pd.DataFrame([v for v in val if v]).sort_values("p")
    val["bh_fdr"] = _bh_fdr(val.p.values)
    return val

@functools.lru_cache(maxsize=1)
def cotarget_dependency_annotated():
    gs = _good_surf().merge(gene_dependency_summary(), on="gene", how="left")
    gs["dependency_class"] = np.where(
        (gs.mean_dependency >= cfg.DEP_ESSENTIAL) | (gs.common_essential == 1), "essential_driver-like",
        np.where(gs.mean_dependency.notna(), "dispensable_passenger", "unknown"))
    val = depmap_validation()
    valu = val.drop_duplicates(["gene","tumor_code","arm"])[["gene","tumor_code","arm","delta","p","bh_fdr"]].rename(
        columns={"p":"depmap_p","bh_fdr":"depmap_fdr"})
    return gs.merge(valu, on=["gene","tumor_code","arm"], how="left")
