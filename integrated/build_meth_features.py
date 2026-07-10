#!/usr/bin/env python3
# Assemble downloaded EPIC beta files into per-gene promoter (TSS200/1500) and
# gene-body mean-beta matrices using the manifest probe->gene maps, then derive
# the two transmission-gate features: mean_promoter_meth and meth_rna_corr.
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); import config as cfg
CACHE = HERE.parent / "data_download/from_source/data/meth_cache"
man = json.load(open(HERE/"reports/meth_manifest.json"))

tss_map  = pd.read_csv(HERE/"reports/meth_probe_tss.csv")    # probe, gene
body_map = pd.read_csv(HERE/"reports/meth_probe_body.csv")
tss_probes  = tss_map.groupby("probe").gene.apply(list).to_dict()
body_probes = body_map.groupby("probe").gene.apply(list).to_dict()
keep = set(tss_probes) | set(body_probes)
print(f"[build] probes of interest: {len(keep)} (tss {len(tss_probes)}, body {len(body_probes)})", flush=True)

tss_cols, body_cols, cases = {}, {}, []
n_ok = 0
for cid, rec in man.items():
    f = CACHE / rec["file_name"]
    if not f.exists() or f.stat().st_size < 1000: continue
    b = pd.read_csv(f, sep="\t", header=None, names=["probe","beta"], na_values=["NA"])
    b = b[b.probe.isin(keep)].dropna()
    beta = dict(zip(b.probe, b.beta))
    # per-gene mean over that gene's TSS / body probes for this case
    tg, bg = {}, {}
    for p, v in beta.items():
        for g in tss_probes.get(p, ()):  tg.setdefault(g, []).append(v)
        for g in body_probes.get(p, ()): bg.setdefault(g, []).append(v)
    tss_cols[cid]  = {g: np.mean(v) for g, v in tg.items()}
    body_cols[cid] = {g: np.mean(v) for g, v in bg.items()}
    cases.append(cid); n_ok += 1
    if n_ok % 50 == 0: print(f"[build] assembled {n_ok} cases", flush=True)

tss  = pd.DataFrame(tss_cols)     # genes x cases
body = pd.DataFrame(body_cols)
print(f"[build] tss matrix {tss.shape}, body {body.shape}", flush=True)
tss.to_parquet(cfg.TAB / "methdata_tss.parquet")
body.to_parquet(cfg.TAB / "methdata_genes.parquet")

# features: mean promoter methylation, and per-gene corr(promoter meth, rna) across cases
so = pd.read_csv(cfg.PATHS["str_omic"], usecols=["gene","caseid","rna"], low_memory=False)
rna = so.pivot_table(index="gene", columns="caseid", values="rna", aggfunc="mean")
shared = [c for c in tss.columns if c in rna.columns]
print(f"[build] cases shared with RNA: {len(shared)}", flush=True)
feats = []
for g in tss.index:
    mp = tss.loc[g, shared].astype(float)
    val = mp.dropna()
    row = {"gene": g, "mean_promoter_meth": float(val.mean()) if len(val) else np.nan}
    if g in rna.index and len(val) >= 20:
        rr = rna.loc[g, val.index].astype(float)
        ok = rr.notna() & val.notna()
        if ok.sum() >= 20 and val[ok].std() > 0 and rr[ok].std() > 0:
            row["meth_rna_corr"] = float(np.corrcoef(val[ok], rr[ok])[0,1])
    feats.append(row)
fdf = pd.DataFrame(feats)
fdf.to_csv(cfg.DIR_TAB / "meth_gene_features.csv", index=False)
print(f"[build] meth_gene_features: {len(fdf)} genes; "
      f"mean_promoter_meth non-null {fdf.mean_promoter_meth.notna().sum()}, "
      f"meth_rna_corr non-null {fdf.meth_rna_corr.notna().sum()}", flush=True)
print(f"[build] median promoter meth {fdf.mean_promoter_meth.median():.3f}; "
      f"median meth_rna_corr {fdf.meth_rna_corr.median():.3f}", flush=True)
