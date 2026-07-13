#!/usr/bin/env python3
"""00d_gene_features.py -- assemble the gene-property predictor feature table from
data_download/from_source, replacing the committed gene_intrinsic feature table.

40 predictor features in 8 groups (all gene-intrinsic; no CPTAC tumour data):
  G1 dosage      gnomAD LOEUF/pLI/mis_z (raw/features/gnomad...) + DepMap dep_*
  G2 complex     CORUM (data/parsed/corum_humanComplexes.txt)
  G3 biophysics  UniProt sequence -> ProtParam (length/MW/pI/GRAVY/aggregation);
                 UniProt features -> tm_domain_count/signal_peptide;
                 DescribePROT disorder/SS -> pinned snapshot
  G5 mRNA        Ensembl biomart transcript length/GC/isoforms (+UTR/codon: NaN if absent)
  G6 evolution   dN/dS + gene_age_proxy -> pinned snapshot; phylop_mean -> NaN stub
  G7 function    UniProt keywords/EC/GO -> is_tf/kinase/receptor/enzyme + go_mf_category
  G8 breadth     GTEx v8 median TPM -> n_tissues_expressed, tau
  G9 network     STRING v12 -> degree/weighted_degree/betweenness (score>=700)

Emits integrated/tables/gene_feature_table.csv keyed on the transmissibility atlas genes.
"""
import sys, gzip, re
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

FEAT = cfg.RAW / "features"
STRING_SCORE_MIN = 700
KD = dict(zip("ACDEFGHIKLMNPQRSTVWY",
    [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
MW = dict(zip("ACDEFGHIKLMNPQRSTVWY",
    [71.08,103.14,115.09,129.12,147.18,57.05,137.14,113.16,128.17,113.16,131.19,114.10,97.12,128.13,156.19,87.08,101.10,99.13,186.21,163.18]))


def _protparam(seq):
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
    n = len(seq)
    if n == 0:
        return dict(length=0, mol_weight=np.nan, isoelectric_point=np.nan, gravy=np.nan, aggregation_propensity=np.nan)
    mw = sum(MW.get(a, 110) for a in seq) + 18.02
    gravy = sum(KD.get(a, 0) for a in seq) / n
    # simple pI via charge-balance bisection
    pk = {"Cterm":3.55,"Nterm":7.5,"D":4.05,"E":4.45,"C":9.0,"Y":10.0,"H":5.98,"K":10.0,"R":12.0}
    cnt = {a: seq.count(a) for a in "DECYHKR"}
    def charge(pH):
        pos = 1/(1+10**(pH-pk["Nterm"])) + cnt["K"]/(1+10**(pH-pk["K"])) + cnt["R"]/(1+10**(pH-pk["R"])) + cnt["H"]/(1+10**(pH-pk["H"]))
        neg = 1/(1+10**(pk["Cterm"]-pH)) + cnt["D"]/(1+10**(pk["D"]-pH)) + cnt["E"]/(1+10**(pk["E"]-pH)) + cnt["C"]/(1+10**(pk["C"]-pH)) + cnt["Y"]/(1+10**(pk["Y"]-pH))
        return pos - neg
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo+hi)/2
        if charge(mid) > 0: lo = mid
        else: hi = mid
    pI = (lo+hi)/2
    # aggregation propensity proxy: mean hydrophobicity of runs (GRAVY-like, positive scale)
    agg = sum(max(KD.get(a,0),0) for a in seq)/n * 10
    return dict(length=n, mol_weight=round(mw,0), isoelectric_point=round(pI,4),
                gravy=round(gravy,6), aggregation_propensity=round(agg,6))


def build_biophysics_function():
    up = pd.read_csv(FEAT/"uniprot_human_reviewed.tsv.gz", sep="\t")
    up = up.rename(columns={"Gene Names (primary)":"gene","Sequence":"seq",
        "Keywords":"kw","Transmembrane":"tm","Signal peptide":"sig",
        "Gene Ontology (molecular function)":"gomf","EC number":"ec","Protein names":"pname"})
    up = up.dropna(subset=["gene","seq"]).drop_duplicates("gene")
    pp = up["seq"].apply(_protparam).apply(pd.Series)
    up = pd.concat([up.reset_index(drop=True), pp.reset_index(drop=True)], axis=1)
    up["tm_domain_count"] = up["tm"].fillna("").str.count("TRANSMEM")
    up["signal_peptide"] = up["sig"].fillna("").str.contains("SIGNAL")
    kw = up["kw"].fillna("").str.lower()
    up["is_kinase"] = kw.str.contains("kinase")
    up["is_receptor"] = kw.str.contains("receptor")
    up["is_tf"] = kw.str.contains("transcription") & kw.str.contains("dna-binding|activator|repressor|regulation")
    up["is_enzyme"] = up["ec"].notna()
    up["go_mf_category"] = _gomf_category(up["gomf"].fillna(""))
    cols = ["gene","length","mol_weight","isoelectric_point","gravy","aggregation_propensity",
            "tm_domain_count","signal_peptide","is_tf","is_kinase","is_receptor","is_enzyme","go_mf_category"]
    return up[cols]


def _gomf_category(gomf_series):
    def cat(s):
        s = s.lower()
        if "transcription" in s and ("dna-binding" in s or "regulat" in s): return "transcription_regulator_activity"
        if "transporter" in s or "channel" in s: return "transporter_activity"
        if "receptor" in s or "transducer" in s: return "molecular_transducer_activity"
        if "catalytic" in s or "kinase" in s or "hydrolase" in s or "transferase" in s or "oxidoreductase" in s: return "catalytic_activity"
        if "structural" in s: return "structural_molecule_activity"
        if "binding" in s: return "binding"
        if "regulator" in s: return "molecular_function_regulator"
        return "other_molecular_function"
    return gomf_series.apply(cat)


def build_dosage():
    gn = pd.read_csv(FEAT/"gnomad.v4.1.constraint_metrics.tsv", sep="\t",
        usecols=["gene","mane_select","canonical","lof.oe_ci.upper","lof.pLI","mis.z_score"])
    gn = gn[(gn.mane_select==True)|(gn.canonical==True)]
    gn = gn.sort_values("mane_select", ascending=False).drop_duplicates("gene")
    gn = gn.rename(columns={"lof.oe_ci.upper":"gnomad_LOEUF","lof.pLI":"gnomad_pLI","mis.z_score":"gnomad_mis_z"})
    dos = gn[["gene","gnomad_LOEUF","gnomad_pLI","gnomad_mis_z"]]
    # DepMap dependency
    dp = pd.read_csv(cfg.PATHS["dep_prob"], index_col=0)
    genes = [c.split(" (")[0] for c in dp.columns]
    dep_mean = 1 - dp.mean(axis=0).values      # dependency prob mean -> effect proxy
    frac_dep = (dp > 0.5).mean(axis=0).values
    depdf = pd.DataFrame({"gene":genes,"dep_frac_dependent":frac_dep})
    # dep_mean_effect: use CRISPRGeneEffect if present, else derive from probability
    ge = cfg.DEPMAP/"CRISPRGeneEffect.csv"
    if ge.exists():
        g2 = pd.read_csv(ge, index_col=0)
        gg = [c.split(" (")[0] for c in g2.columns]
        depdf = depdf.merge(pd.DataFrame({"gene":gg,"dep_mean_effect":g2.mean(axis=0).values}), on="gene", how="outer")
    else:
        depdf["dep_mean_effect"] = -depdf["dep_frac_dependent"]
    return dos.merge(depdf, on="gene", how="outer")


def build_complex():
    cx = pd.read_csv(cfg.PARSED/"corum_humanComplexes.txt", sep="\t")
    subcol = [c for c in cx.columns if "subunits" in c.lower() and "gene" in c.lower()]
    subcol = subcol[0] if subcol else [c for c in cx.columns if "gene name" in c.lower()][0]
    from collections import defaultdict
    ncx = defaultdict(int); sizes = defaultdict(list)
    for _, r in cx.iterrows():
        subs = [s.strip() for s in re.split(r"[;,]", str(r[subcol])) if s.strip()]
        for g in set(subs):
            ncx[g] += 1; sizes[g].append(len(subs))
    rows = [{"gene":g,"n_complexes":ncx[g],"complex_size":np.mean(sizes[g]),
             "in_complex":True,"has_complex":1} for g in ncx]
    return pd.DataFrame(rows)


def build_mrna():
    mf = pd.read_csv(FEAT/"ensembl_mrna_features.tsv", sep="\t")
    mf.columns = ["gene","transcript_length","gc_content","biotype"]
    mf = mf[mf.biotype=="protein_coding"]
    tl = mf.groupby("gene")["transcript_length"].max().reset_index()  # longest isoform proxy
    gc = mf.groupby("gene")["gc_content"].first().reset_index()
    tc = pd.read_csv(FEAT/"ensembl_transcript_counts.tsv", sep="\t")
    tc.columns = ["gene","tx"]
    niso = tc.groupby("gene")["tx"].nunique().reset_index().rename(columns={"tx":"n_isoforms"})
    out = tl.merge(gc, on="gene", how="outer").merge(niso, on="gene", how="outer")
    # UTR lengths + codon optimality: Ensembl-derived but the bulk BioMart UTR/CDS
    # query is unreliable at whole-transcriptome scale; carried in the pinned snapshot
    # (external_predicted_descriptors_snapshot.csv) alongside dn_ds/gene_age.
    return out


def build_breadth():
    import gzip as _gz
    gct = cfg.RAW/"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    with _gz.open(gct, "rt") as f:
        f.readline(); f.readline()
        gx = pd.read_csv(f, sep="\t")
    gx["gene"] = gx["Description"]
    tis = gx.drop(columns=[c for c in ["Name","Description","gene"] if c in gx.columns])
    mat = tis.values
    n_expr = (mat > 1).sum(axis=1)
    mx = mat.max(axis=1, keepdims=True); mx[mx==0] = 1
    tau = (1 - mat/mx).sum(axis=1) / (mat.shape[1]-1)
    return pd.DataFrame({"gene":gx["gene"].values,"n_tissues_expressed":n_expr,"tau":np.round(tau,6)}).drop_duplicates("gene")


def build_network():
    return pd.read_parquet(FEAT/"string_network_centrality.parquet")


def build_snapshot():
    return pd.read_csv(FEAT/"external_predicted_descriptors_snapshot.csv")


def main():
    atlas = pd.read_csv(cfg.DIR_TAB/"transmissibility_atlas.csv")
    base = atlas[["gene","observed_transmissibility","n_amplified_cases"]].rename(
        columns={"observed_transmissibility":"transmissibility"})
    parts = [build_dosage(), build_complex(), build_biophysics_function(),
             build_mrna(), build_breadth(), build_network(), build_snapshot()]
    ft = base
    for p in parts:
        ft = ft.merge(p, on="gene", how="left")
    ft["phylop_mean"] = np.nan  # documented unavailable stub (all-null in GI too)
    for c in ["in_complex","n_complexes","complex_size","has_complex"]:
        if c in ft.columns:
            ft[c] = ft[c].fillna({"in_complex":False,"n_complexes":0,"has_complex":0,"complex_size":np.nan}.get(c))
    out = cfg.DIR_TAB/"gene_feature_table.csv"
    ft.to_csv(out, index=False)
    print(f"gene_feature_table.csv: {ft.shape[0]} genes x {ft.shape[1]} cols")
    return ft


if __name__ == "__main__":
    main()
