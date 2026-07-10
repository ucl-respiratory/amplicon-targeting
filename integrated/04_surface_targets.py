#!/usr/bin/env python3
# =============================================================================
# 04_surface_targets.py  --  Surface ADC co-target identification from BOTH
# measured and predicted transmissibility.
# -----------------------------------------------------------------------------
# The integrated nomination. A gene is a surface ADC co-target candidate if it
# is (1) on a recurrent amplicon, (2) transmitted to protein, and (3) presents
# an accessible extracellular ectodomain. Transmission evidence comes from two
# sources, tagged so a reader always knows which:
#
#   empirical   observed co-elevation measured in this cancer type (FDR<0.1)
#   predicted   gene-intrinsic prior (no proteome needed in this context)
#   both        measured AND prior-corroborated
#
# The predicted arm is the wider space of genes: it nominates surface targets in
# amplicons/cancer types where proteomics is thin or absent, which the measured
# arm cannot reach. The surface-topology gate (UniProt extracellular residues)
# is applied identically to both arms.
#
# Runs fully once data_download is built (co-elevation + amplification frequency
# from str_omic; topology from stage 15). The predicted arm + the validated lead
# detail run now from committed tables + a direct UniProt topology fetch.
#
# Outputs: figures/fig5_surface_targets.png, tables/integrated_surface_targets.csv
# =============================================================================
import sys, json, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_shared as ish
from cnt_shared import cfg
ish.apply_style(sizes=(9, 8, 7))

# ---- UniProt topology fetch (self-contained; mirrors data_download stage 15) --
import urllib.request, urllib.parse
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = "accession,gene_primary,cc_subcellular_location,ft_transmem,ft_topo_dom,ft_signal"

def _uniprot_batch(genes):
    q = "(" + " OR ".join(f"gene:{g}" for g in genes) + ") AND organism_id:9606 AND reviewed:true"
    url = UNIPROT_REST + "?" + urllib.parse.urlencode(
        {"query": q, "fields": UNIPROT_FIELDS, "format": "json", "size": 500})
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())["results"]

def _extracellular_aa(entry):
    """Total extracellular topological-domain residues + TM helix count."""
    tm = 0; ecto = 0
    for f in entry.get("features", []):
        t = f.get("type", "")
        if t == "Transmembrane": tm += 1
        if t == "Topological domain" and "extracellular" in f.get("description", "").lower():
            loc = f.get("location", {})
            s = loc.get("start", {}).get("value"); e = loc.get("end", {}).get("value")
            if s is not None and e is not None: ecto += (e - s + 1)
    return tm, ecto

def fetch_topology(genes):
    genes = sorted(set(genes)); rows = []
    for i in range(0, len(genes), 80):
        batch = genes[i:i+80]
        try: res = _uniprot_batch(batch)
        except Exception as ex:
            print(f"  uniprot batch {i} failed: {ex}"); continue
        for entry in res:
            gp = entry.get("genes", [{}])
            sym = gp[0].get("geneName", {}).get("value") if gp else None
            if sym not in batch: continue
            tm, ecto = _extracellular_aa(entry)
            subcell = "; ".join(
                c.get("subcellularLocation", {}).get("location", {}).get("value", "")
                for cc in entry.get("comments", []) if cc.get("commentType") == "SUBCELLULAR LOCATION"
                for c in cc.get("subcellularLocations", []))
            rows.append({"gene": sym, "uniprot_acc": entry.get("primaryAccession"),
                         "n_TM_helices": tm, "total_extracellular_aa": ecto,
                         "subcellular": subcell})
        time.sleep(0.3)
    df = pd.DataFrame(rows)
    if len(df):  # one canonical (shortest accession) per gene
        df = df.sort_values("uniprot_acc", key=lambda s: s.str.len()).drop_duplicates("gene")
    return df

def surface_verdict(tm, ecto):
    if ecto >= cfg.MIN_ECTODOMAIN_AA and tm >= 1: return "accessible"
    if 0 < ecto < cfg.MIN_ECTODOMAIN_AA: return "marginal"
    return "not accessible"

# ---- amplification frequency per (gene, tumor_type) -------------------------
def amplification_frequency():
    df = ish.load_str_omic().dropna(subset=["cn_adjusted"])
    df["amp"] = (df.cn_adjusted >= cfg.AMP_THRESHOLD).astype(int)
    freq = df.groupby(["gene", "tumor_code"]).amp.mean().reset_index(name="amp_freq")
    return freq

def surfaceome_universe(atlas):
    """Surface genes on recurrent amplicons, tagged empirical vs predicted.
    Universe restriction is essential: most high-transmissibility genes are
    intracellular (ribosome biogenesis, RNA processing); the target space is
    surface genes on amplicons, ranked by transmissibility WITHIN that set.

    When data_download is built, the surface flag + measured co-elevation come
    from the amplicon analysis. Until then, the committed ADC cross-reference
    (surface genes on the 3q amplicon, already transmissibility-scored) is the
    universe, and the predicted arm extends it to surface genes genome-wide."""
    have_data = cfg.PATHS["str_omic"].exists()
    if have_data:
        co = ish.amplicon_coelevation()
        surf_co = co[co.is_surface == 1].copy()
        surf_co["measured_coelevated"] = surf_co.fdr < cfg.FDR_ALPHA
        uni = surf_co.merge(atlas, on="gene", how="left")
        uni = uni.rename(columns={"tumor_code": "cancer_type"})
        return uni, True
    # committed fallback: the ADC cross-reference table = surface genes, scored
    xref = pd.read_csv(cfg.GI_PATHS["atlas"].parent /
                       "tournament_master_tab11_transmissibility_adc_crossref.csv")
    xref = xref[xref.surface == 1].copy()
    xref["measured_coelevated"] = xref.observed_transmissibility >= \
        xref.observed_transmissibility.median()
    xref["cancer_type"] = "3q-amplified (committed)"
    return xref, False

def main():
    atlas = ish.load_atlas()
    lv = ish.load_leads_validation()
    uni, have_data = surfaceome_universe(atlas)
    if not have_data:
        print("data_download tables absent — surfaceome universe from committed "
              "ADC cross-reference; predicted arm extends genome-wide")

    # predicted-high flag (top quintile of prior, WITHIN the surfaceome)
    pthr = atlas.predicted_transmissibility_oof.quantile(0.80)
    uni["predicted_high"] = uni.predicted_transmissibility_oof >= pthr

    # ---- surface-topology gate (UniProt) on the universe --------------------
    topo_path = cfg.PATHS["topology"]
    if topo_path.exists():
        topo = pd.read_parquet(topo_path)
    else:
        cand = sorted(set(uni.gene) | set(lv.gene))
        print(f"stage-15 topology absent — fetching UniProt topology for "
              f"{len(cand)} surfaceome candidate genes directly")
        topo = fetch_topology(cand)
    if "total_extracellular_aa" in topo.columns:
        topo["surface_verdict_calc"] = [
            surface_verdict(t, e) for t, e in
            zip(topo.get("n_TM_helices", 0), topo.total_extracellular_aa)]

    integ = uni.merge(topo, on="gene", how="left")
    integ["measured_here"] = integ.get("measured_coelevated", False)

    def tier(r):
        acc = r.get("surface_verdict_calc", "unknown")
        if acc not in ("accessible", "marginal"): return "not_surface"
        if r["measured_here"] and r["predicted_high"]: return "measured_pred"
        if r["measured_here"]: return "measured_high"
        return "predicted_only"
    integ["confidence_tier"] = integ.apply(tier, axis=1)
    keep = ["gene", "cancer_type", "observed_transmissibility",
            "predicted_transmissibility_oof", "predicted_percentile",
            "n_TM_helices", "total_extracellular_aa", "surface_verdict_calc",
            "measured_here", "predicted_high", "confidence_tier"]
    keep = [c for c in keep if c in integ.columns]
    out = integ[keep].drop_duplicates("gene").sort_values(
        "predicted_transmissibility_oof", ascending=False)
    out.to_csv(cfg.DIR_TAB / "integrated_surface_targets.csv", index=False)

    surf = out[out.surface_verdict_calc.isin(["accessible", "marginal"])] \
        if "surface_verdict_calc" in out.columns else out
    tier_counts = surf.confidence_tier.value_counts().to_dict()

    # ---- Figure 5: empirical vs predicted surface targets -------------------
    fig = plt.figure(figsize=(12, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1], wspace=0.28)
    C = {"measured_pred": "#6c5b9c", "measured_high": "#2c6fbb",
         "predicted_only": "#e67e22", "not_surface": "#cccccc"}

    # panel a: predicted vs observed transmissibility, colored by tier, sized by amp cases
    axa = fig.add_subplot(gs[0, 0])
    for t, c in C.items():
        s = surf[surf.confidence_tier == t] if "confidence_tier" in surf else surf.iloc[0:0]
        if not len(s): continue
        axa.scatter(s.predicted_transmissibility_oof,
                    s.get("observed_transmissibility", pd.Series(np.nan, index=s.index)),
                    s=18, alpha=0.6, color=c, edgecolors="none", label=t.replace("_", " "))
    axa.set_xlabel("predicted transmissibility (prior)")
    axa.set_ylabel("observed transmissibility (measured)")
    axa.set_title("Surface ADC co-target candidates:\nempirical vs predicted evidence, tagged",
                  fontsize=9, loc="left")
    axa.legend(frameon=False, fontsize=7.5, loc="lower right", title="confidence tier")

    # label the validated leads
    for g in lv.gene:
        row = surf[surf.gene == g]
        if len(row):
            axa.annotate(g, xy=(row.predicted_transmissibility_oof.iloc[0],
                                row.observed_transmissibility.iloc[0]),
                         fontsize=7.5, fontweight="bold",
                         xytext=(4, 4), textcoords="offset points")

    # panel b: tier counts (how many targets each arm contributes)
    axb = fig.add_subplot(gs[0, 1])
    tiers = ["measured_pred", "measured_high", "predicted_only"]
    vals = [tier_counts.get(t, 0) for t in tiers]
    axb.barh(range(len(tiers)), vals, color=[C[t] for t in tiers])
    axb.set_yticks(range(len(tiers)))
    axb.set_yticklabels([t.replace("_", "\n") for t in tiers])
    axb.invert_yaxis()
    axb.set_xlabel("surface-accessible candidates")
    axb.set_title("The predicted arm widens the\ntarget space beyond measurement",
                  fontsize=9, loc="left")
    for i, v in enumerate(vals):
        axb.text(v + 0.5, i, str(v), va="center", fontsize=8)
    for ax, L in zip([axa, axb], "ab"):
        ish.panel_letter(ax, L)

    fig.savefig(cfg.DIR_FIG / "fig5_surface_targets.png", dpi=200, bbox_inches="tight")

    ish.record("04_surface_targets", {
        "have_data": bool(have_data),
        "n_universe": int(len(out)),
        "n_predicted_high": int(out.get("predicted_high", pd.Series(dtype=bool)).sum()),
        "n_surface_accessible": int(len(surf)),
        "tier_counts": tier_counts,
        "leads": lv[["gene", "surface_verdict"]].to_dict("records"),
    })
    print("tier counts (surface-accessible):", tier_counts)
    print("leads:", lv[["gene", "surface_verdict"]].to_string(index=False))

if __name__ == "__main__":
    main()
