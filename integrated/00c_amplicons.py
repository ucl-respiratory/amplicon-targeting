# =============================================================================
# 00c_amplicons.py -- pan-cancer recurrent-amplicon landscape + co-elevation
# -----------------------------------------------------------------------------
# Detects recurrent cytoband amplicons across ALL CPTAC cohorts in str_omic,
# then runs the paper-exact per-(tumor, band, gene) one-sided Fisher co-elevation
# test (high tissue-relative protein enriched in band-amplified vs non-amplified
# tumours). This is the candidate universe for surface-antigen nomination and is
# pan-cancer by construction -- not restricted to any tumour type.
#
# Amplified gene   : cn_adjusted >= AMP_THRESHOLD (1.4)
# Amplified band   : >= COELEV_MIN_FRAC (50%) of the band's genes amplified in a tumour
# Recurrent band   : amplified in >= RECURRENCE_FREQ (20%) AND >= RECURRENCE_MIN_N (8) samples
# Co-elevated gene : Fisher(prot.rel.tissue>0.8 | band-amplified vs not), BH-FDR < FDR_ALPHA
#
# Emits: tables/recurrent_amplicons.csv, tables/amplicon_coelevation.csv
#        figures/fig_amplicon_landscape.png
# =============================================================================
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cnt_shared as ish
import config as cfg

def _bh_fdr(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    ranked = np.empty(n); cummin = 1.0
    for i in range(n-1, -1, -1):
        idx = order[i]; val = p[idx]*n/(i+1); cummin = min(cummin, val); ranked[idx] = cummin
    return np.clip(ranked, 0, 1)

def recurrent_amplicons(d):
    band_samp = d.groupby(["tumor_code","cytogenetic_location","caseid"]).agg(
        frac_genes_amp=("amplified","mean")).reset_index()
    band_samp["band_amplified"] = (band_samp.frac_genes_amp >= cfg.COELEV_MIN_FRAC).astype(int)
    br = band_samp.groupby(["tumor_code","cytogenetic_location"]).agg(
        n_amp_samples=("band_amplified","sum"), n_samples=("band_amplified","size")).reset_index()
    br["recurrence"] = br.n_amp_samples / br.n_samples
    br["recurrent"] = (br.recurrence >= cfg.RECURRENCE_FREQ) & (br.n_amp_samples >= cfg.RECURRENCE_MIN_N)
    return br, band_samp

def coelevation(d, band_samp, recurrent_keys):
    d2 = d.merge(band_samp[["tumor_code","cytogenetic_location","caseid","band_amplified"]],
                 on=["tumor_code","cytogenetic_location","caseid"], how="left")
    d2["is_rec"] = [(t,c) in recurrent_keys for t,c in zip(d2.tumor_code, d2.cytogenetic_location)]
    dd = d2[d2.is_rec].copy()
    recs = []
    for (tc, band, gene), g in dd.groupby(["tumor_code","cytogenetic_location","gene"]):
        amp = g[g.band_amplified == 1]; noamp = g[g.band_amplified == 0]
        if len(amp) < 8 or len(noamp) < 5: continue
        a = int(amp.prot_high.sum()); b = len(amp)-a
        c = int(noamp.prot_high.sum()); e = len(noamp)-c
        try: _, pval = fisher_exact([[a,b],[c,e]], alternative="greater")
        except Exception: pval = np.nan
        recs.append(dict(tumor_code=tc, band=band, gene=gene,
                         p_high_amp=amp.prot_high.mean(), p_high_noamp=noamp.prot_high.mean(),
                         lift=amp.prot_high.mean()-noamp.prot_high.mean(),
                         n_amp=len(amp), n_noamp=len(noamp), fisher_p=pval,
                         is_surface=int(g.is_surface.iloc[0]) if "is_surface" in g else 0,
                         is_secreted=int(g.is_secreted.iloc[0]) if "is_secreted" in g else 0))
    co = pd.DataFrame(recs).dropna(subset=["fisher_p"]).reset_index(drop=True)
    co["fdr"] = _bh_fdr(co.fisher_p.values)
    return co

def main():
    ish.apply_style(sizes=(9,8,7))
    df = ish.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","cytogenetic_location","prot.rel.tissue"]).copy()
    d["amplified"] = (d.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)
    d["prot_high"] = (d["prot.rel.tissue"] > cfg.REL_TISSUE_HI).astype(int)

    br, band_samp = recurrent_amplicons(d)
    rec = br[br.recurrent].copy()
    rec_keys = set(zip(rec.tumor_code, rec.cytogenetic_location))
    print(f"recurrent cytoband amplicons: {len(rec)} across {rec.tumor_code.nunique()} cohorts")
    print(rec.groupby("tumor_code").size().to_dict())

    co = coelevation(d, band_samp, rec_keys)
    co_sig = co[co.fdr < cfg.FDR_ALPHA]
    print(f"co-elevation tests: {len(co)}; FDR<{cfg.FDR_ALPHA}: {len(co_sig)} "
          f"({co_sig.is_surface.sum()} surface, {co_sig.is_secreted.sum()} secreted)")

    rec.to_csv(cfg.DIR_TAB / "recurrent_amplicons.csv", index=False)
    co.to_csv(cfg.DIR_TAB / "amplicon_coelevation.csv", index=False)

    # ---- figure: recurrent amplicons per cohort ----
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    cohort_counts = rec.groupby("tumor_code").size().sort_values(ascending=True)
    coelev_by_cohort = co_sig.groupby("tumor_code").gene.nunique().reindex(cohort_counts.index).fillna(0)
    y = np.arange(len(cohort_counts))
    ax.barh(y, cohort_counts.values, color="#4C72B0", label="recurrent amplicons")
    for i,(n_amp,n_ce) in enumerate(zip(cohort_counts.values, coelev_by_cohort.values)):
        ax.text(n_amp+1, i, f"{int(n_ce)} co-elevated genes (FDR<0.1)", va="center", fontsize=7, color="#333")
    ax.set_yticks(y); ax.set_yticklabels(cohort_counts.index)
    ax.set_xlabel("recurrent cytoband amplicons"); ax.set_title("Pan-cancer recurrent-amplicon landscape")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(cfg.DIR_FIG / "fig_amplicon_landscape.png", dpi=200)

    ish.record("00c_amplicons", {
        "n_recurrent_amplicons": int(len(rec)),
        "n_cohorts": int(rec.tumor_code.nunique()),
        "per_cohort": rec.groupby("tumor_code").size().to_dict(),
        "n_coelev_tests": int(len(co)),
        "n_coelev_sig": int(len(co_sig)),
        "note_surface": "surface accessibility is annotated downstream via live UniProt topology "
                        "(str_omic is_surface is empty because the TCSA surfaceome file was unreachable)",
        "cohorts_no_coelev": "GBM (no tissue-relative proteome), PDA (no recurrent cytoband amplicon meets 20%/8-sample bar)",
        "fdr_alpha": cfg.FDR_ALPHA, "amp_threshold": cfg.AMP_THRESHOLD,
    })

if __name__ == "__main__":
    main()
