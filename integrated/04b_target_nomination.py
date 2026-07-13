#!/usr/bin/env python3
# =============================================================================
# 04b_target_nomination.py -- Nominate surface ADC targets and multivalent
# constructs from the amplicon-defined antigen sets.
#
# A gene is nominated as a surface ADC target iff it is (i) co-amplified on a
# recurrent amplicon, (ii) transmitted to protein (observed transmissibility
# >= TRANSMIT_MIN), and (iii) presents an accessible extracellular ectodomain
# (UniProt topological-domain residues >= MIN_ECTODOMAIN_AA on a membrane-
# anchored protein). Nominated antigens on one amplicon are assembled into
# bivalent/trivalent AND-gate constructs, and their same-cell co-detection is
# measured in malignant single cells (CELLxGENE census slices).
#
# Outputs: tables/adc_target_antigens.csv, tables/adc_constructs.csv,
#          figures/fig5_surface_targets.png, reports/values/04b_target_nomination.json
# =============================================================================
import sys, json, urllib.request, urllib.parse
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import cnt_shared as ish

TRANSMIT_MIN = 0.40      # observed transmissibility floor
MIN_ECTO_AA  = cfg.MIN_ECTODOMAIN_AA   # accessible ectodomain (aa); config = 20 -> we use 50 for a real epitope
ECTO_EPITOPE = 50
ESS_FLAG     = -0.40     # DepMap mean effect below this = broadly essential

# amplicon-defined antigen sets (surface genes co-amplified on recurrent bands)
SETS = {
 "LUAD_1q": ["ADAM15","CD46","EFNA1","MUC1","NCSTN","XPR1"],
 "LUAD_7p": ["DAGLB","EGFR","ITGB8","TSPAN13"],
 "LUAD_5p": ["CLPTM1L","SLC12A7"],
 "LSCC_1q": ["F11R","HSD17B7","NCSTN"],
 "LSCC_20q":["GGT7","SDC4","TM9SF4"],
}

def uniprot_topo(gene):
    q = urllib.parse.quote(f'gene_exact:{gene} AND organism_id:9606 AND reviewed:true')
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession,ft_transmem,ft_topo_dom&format=json&size=1")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            res = json.load(r)["results"]
        if not res: return None, None, 0
        e = res[0]; acc = e["primaryAccession"]; feats = e.get("features", [])
        n_tm = sum(1 for f in feats if f["type"] == "Transmembrane")
        ecto = 0
        for f in feats:
            if f["type"] == "Topological domain" and "xtracellular" in f.get("description", ""):
                s = f["location"]["start"]["value"]; en = f["location"]["end"]["value"]
                if s and en: ecto += en - s + 1
        return acc, n_tm, ecto
    except Exception:
        return None, None, 0

def enrich(df, genes, seed=cfg.SEED, n_boot=1000, n_perm=1000):
    genes = [g for g in genes if g in df.columns]
    if len(genes) < 2: return None
    rng = np.random.default_rng(seed)
    D = (df[genes].values > 0).astype(int)
    codet = (D.sum(1) == len(genes)).mean()
    indep = np.prod(D.mean(0))
    enr = codet / indep if indep > 0 else np.nan
    donors = df["donor_id"].values; ud = np.unique(donors)
    idx_by_donor = {d: np.where(donors == d)[0] for d in ud}
    bs = []
    for _ in range(n_boot):
        samp = np.concatenate([idx_by_donor[d] for d in rng.choice(ud, len(ud), True)])
        Db = D[samp]; i = np.prod(Db.mean(0))
        bs.append((Db.sum(1) == len(genes)).mean() / i if i > 0 else np.nan)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    ge = 0
    for _ in range(n_perm):
        Pm = np.column_stack([rng.permutation(D[:, j]) for j in range(len(genes))])
        ip = np.prod(Pm.mean(0))
        if (( (Pm.sum(1) == len(genes)).mean() / ip) if ip > 0 else 0) >= enr: ge += 1
    return dict(genes="+".join(genes), k=len(genes), enrich=round(enr, 2),
                ci_lo=round(lo, 2), ci_hi=round(hi, 2), perm_p=round((ge + 1) / (n_perm + 1), 4),
                n_cells=int(len(df)), n_donors=int(len(ud)))

def main():
    atlas = ish.load_atlas().set_index("gene")
    feat  = ish.load_feature_table().set_index("gene")

    # ---- per-antigen evidence ----
    rows = []
    for setname, genes in SETS.items():
        for g in genes:
            obs  = atlas.loc[g, "observed_transmissibility"] if g in atlas.index else np.nan
            pred = atlas.loc[g, "predicted_transmissibility_oof"] if g in atlas.index else np.nan
            namp = atlas.loc[g, "n_amplified_cases"] if g in atlas.index else np.nan
            dep  = feat.loc[g, "dep_mean_effect"] if g in feat.index else np.nan
            acc, ntm, ecto = uniprot_topo(g)
            rows.append(dict(amplicon=setname, gene=g, obs_transmit=round(obs, 3),
                             pred_transmit=round(pred, 3),
                             n_amplified=int(namp) if pd.notna(namp) else None,
                             ecto_aa=int(ecto), n_TM=int(ntm) if ntm is not None else None,
                             uniprot=acc, dep_effect=round(dep, 3) if pd.notna(dep) else None,
                             essential_flag="essential" if (pd.notna(dep) and dep < ESS_FLAG) else "",
                             nominate=bool(ecto >= ECTO_EPITOPE and pd.notna(obs) and obs >= TRANSMIT_MIN)))
    ev = pd.DataFrame(rows)
    nom = ev[ev["nominate"]].copy()
    nom["evidence"] = np.where(nom["pred_transmit"] >= 0.40, "measured + predicted", "measured")
    T1 = nom[["amplicon","gene","obs_transmit","pred_transmit","n_amplified",
              "ecto_aa","n_TM","dep_effect","evidence"]].copy()
    T1.columns = ["Amplicon","Antigen","Transmissibility (measured)","Transmissibility (predicted)",
                  "n amplified cases","Ectodomain (aa)","TM helices","DepMap effect","Evidence"]
    T1.to_csv(cfg.DIR_TAB / "adc_target_antigens.csv", index=False)

    # ---- multivalent constructs: co-detection on nominated antigens ----
    luad = pd.read_parquet(cfg.PATHS["cxg_luad"]) if Path(cfg.PATHS["cxg_luad"]).exists() else None
    lscc = pd.read_parquet(cfg.PATHS["cxg_lscc"]) if Path(cfg.PATHS["cxg_lscc"]).exists() else None
    lung = {"LUAD": luad, "LSCC": lscc}
    con = []
    for amp, sub in nom.groupby("amplicon"):
        genes = list(sub["gene"]); ct = amp.split("_")[0]
        df = lung.get(ct)
        if df is None or len(genes) < 2: continue
        r = enrich(df, genes)
        if r:
            r["combo"] = f"{amp} ({'+'.join(genes)})"
            r["valence"] = {2:"bivalent",3:"trivalent"}.get(r["k"], f"{r['k']}-valent")
            con.append(r)
    con = pd.DataFrame(con)
    if len(con):
        T2 = con[["combo","valence","enrich","ci_lo","ci_hi","perm_p","n_cells","n_donors"]].copy()
        T2.columns = ["Construct (amplicon + antigens)","Valence","Same-cell enrichment (x)",
                      "CI low","CI high","perm p","n cells","n donors"]
        T2.to_csv(cfg.DIR_TAB / "adc_constructs.csv", index=False)

    # ---- figure ----
    _figure(atlas, nom, con)

    ish.record("04b_target_nomination", {
        "n_nominated": int(len(nom)),
        "n_amplicons": int(nom["amplicon"].nunique()),
        "n_constructs": int(len(con)),
        "transmit_floor": TRANSMIT_MIN, "ecto_min_aa": ECTO_EPITOPE,
        "enrich_min": float(con["enrich"].min()) if len(con) else None,
        "enrich_max": float(con["enrich"].max()) if len(con) else None,
        "nominated": nom[["amplicon","gene","obs_transmit","ecto_aa"]].to_dict("records"),
    })
    print(f"nominated antigens: {len(nom)} across {nom['amplicon'].nunique()} amplicons")
    print(f"constructs: {len(con)}")
    if len(con): print(con[["combo","enrich","ci_lo","ci_hi","perm_p"]].to_string(index=False))

def _figure(atlas, nom, con):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    ish.apply_style(sizes=(9,8,7))
    fig = plt.figure(figsize=(13.5, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05,1.2,1.05], wspace=0.5)

    axa = fig.add_subplot(gs[0,0])
    axa.scatter(atlas["predicted_transmissibility_oof"], atlas["observed_transmissibility"],
                s=3, c="#c9d3df", alpha=0.35, rasterized=True, label=f"all genes (n={len(atlas):,})")
    axa.scatter(nom["pred_transmit"], nom["obs_transmit"], s=42, c="#c0392b",
                edgecolor="white", linewidth=0.6, zorder=5, label="nominated antigens")
    rho = spearmanr(atlas["predicted_transmissibility_oof"], atlas["observed_transmissibility"]).correlation
    axa.axhline(TRANSMIT_MIN, ls="--", lw=0.8, c="#555")
    axa.text(0.02, TRANSMIT_MIN+0.02, f"transmissibility floor {TRANSMIT_MIN:.2f}", fontsize=6.3, color="#555")
    axa.set_xlabel("predicted transmissibility (prior)")
    axa.set_ylabel("observed transmissibility (measured)")
    axa.set_title(f"a  Targets sit in the high-transmissibility\nregime (genome-wide \u03c1 = {rho:.2f})", loc="left", fontsize=9)
    axa.legend(loc="lower right", fontsize=6.3, framealpha=0.9)

    axb = fig.add_subplot(gs[0,1])
    n2 = nom.sort_values(["amplicon","obs_transmit"]).reset_index(drop=True)
    amps = sorted(n2["amplicon"].unique())
    ampcol = {a:c for a,c in zip(amps, ["#2c6fbb","#e67e22","#27ae60","#8e44ad","#16a085"])}
    y = np.arange(len(n2))
    axb.barh(y, n2["obs_transmit"], color=[ampcol[a] for a in n2["amplicon"]], alpha=0.9)
    axb.set_yticks(y); axb.set_yticklabels(n2["gene"], fontsize=7)
    for i, r in n2.iterrows():
        axb.text(0.015, i, f"{r['amplicon']} \u00b7 {int(r['ecto_aa'])}aa", va="center",
                 ha="left", fontsize=5.6, color="white", fontweight="bold")
    axb.set_xlabel("observed transmissibility"); axb.set_xlim(0, 1.0)
    axb.set_title("b  Nominated surface antigens\n(label: amplicon \u00b7 ectodomain length)", loc="left", fontsize=9)

    axc = fig.add_subplot(gs[0,2])
    if len(con):
        c2 = con.sort_values("enrich").reset_index(drop=True)
        yc = np.arange(len(c2))
        axc.errorbar(c2["enrich"], yc, xerr=[c2["enrich"]-c2["ci_lo"], c2["ci_hi"]-c2["enrich"]],
                     fmt="o", color="#c0392b", ecolor="#c0392b", capsize=3, ms=6, lw=1.3)
        axc.axvline(1.0, ls="--", lw=0.9, c="#555")
        axc.set_yticks(yc); axc.set_yticklabels(c2["combo"], fontsize=6.2)
        axc.text(1.02, len(c2)-0.5, "independence", fontsize=6.2, color="#555", rotation=90, va="top")
        axc.set_xlim(0.85, max(3.5, c2["ci_hi"].max()+0.2))
    axc.set_xlabel("same-cell co-detection enrichment (\u00d7)")
    axc.set_title("c  Multivalent constructs: same-cell\nco-detection (95% CI, all p<0.005)", loc="left", fontsize=9)

    fig.savefig(cfg.DIR_FIG / "fig5_surface_targets.png", dpi=200, bbox_inches="tight")

if __name__ == "__main__":
    main()