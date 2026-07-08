# =============================================================================
# 12_mechanism_decomposition.py  --  Fig 9: Shapley variance decomposition of
# between-gene protein responsiveness by mechanism (dosage transmission,
# regulatory/chromatin, complex buffering, protein structure, GTEx).
#
# Headline result: dosage transmission (CN->RNA) explains ~89% of the explainable
# between-gene variance in protein responsiveness.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from itertools import permutations
from sklearn.linear_model import LinearRegression
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cnt_io, _features, _mechanism
from cnt_io import DIR_TAB, savefig, record_values

STRUCT = ['VSL2_disorder','PSIPRED_helix','PSIPRED_strand','PSIPRED_coil',
 'MMseq2_low_conservation','MMseq2_high_conservation','MMseq2_median','DFLpred_linker',
 'ASAquick_buried','DisoRDPbind_RNA','DisoRDPbind_DNA','DisoRDPbind_PRO','MoRFchibi_morf',
 'DRNApred_RNA','DRNApred_DNA','SCRIBER_PRO']

GROUPS = {
 "Dosage transmission (CN\u2192RNA)": ["cn_rna_corr"],
 "Regulatory/chromatin": ['mean_promoter_meth','meth_rna_corr','log_gene_length','gene_density_1mb'],
 "Complex buffering (CORUM)": ['is_complex_subunit','corum_n_complexes','corum_max_complex_size'],
 "Protein structure/turnover": STRUCT,
 "Tissue expression (GTEx)": ["gtex"],
}

def build_gene_table():
    df = cnt_io.load_str_omic()
    gcols = ['gene'] + STRUCT + ['gtex']
    gl = df[gcols].drop_duplicates('gene').set_index('gene')
    # CORUM features (0 for non-subunits)
    cx = _features.corum_gene_features().set_index('gene')
    gl = gl.join(cx, how='left')
    gl['is_complex_subunit'] = gl['is_complex_subunit'].fillna(0)
    gl['corum_n_complexes'] = gl['corum_n_complexes'].fillna(0)
    gl['corum_max_complex_size'] = gl['corum_max_complex_size'].fillna(0)
    # attenuation (transmission + responsiveness target)
    att = _mechanism.per_gene_attenuation().set_index('gene')
    gl = gl.join(att[['cn_rna_corr','cn_prot_corr','attenuation']], how='inner')
    # regulatory
    reg = _features.regulatory_gene_features().set_index('gene')
    gl = gl.join(reg[['mean_promoter_meth','meth_rna_corr','log_gene_length','gene_density_1mb']], how='left')
    allf = STRUCT + ['gtex','is_complex_subunit','corum_n_complexes','corum_max_complex_size',
                     'cn_rna_corr','mean_promoter_meth','meth_rna_corr','log_gene_length','gene_density_1mb']
    for c in allf:
        gl[c] = gl[c].fillna(gl[c].median())
    return gl.dropna(subset=['cn_prot_corr'])

def r2_of(sub, groups, y, g):
    f = [c for gn in sub for c in groups[gn]]
    if not f: return 0.0
    return LinearRegression().fit(g[f].values, y).score(g[f].values, y)

def shapley(gl):
    y = gl['cn_prot_corr'].values
    gn = list(GROUPS); perms = list(permutations(gn))
    shap = {x: 0.0 for x in gn}
    for perm in perms:
        seen = []
        for x in perm:
            shap[x] += r2_of(seen + [x], GROUPS, y, gl) - r2_of(seen, GROUPS, y, gl)
            seen.append(x)
    for x in shap: shap[x] /= len(perms)
    full = r2_of(gn, GROUPS, y, gl)
    return shap, full

def main():
    gl = build_gene_table()
    shap, full = shapley(gl)
    dec = pd.DataFrame([{"mechanism": k, "shapley_r2": v, "pct": 100*v/full, "total_r2": full}
                        for k, v in shap.items()])
    dec.to_csv(DIR_TAB / "mechanism_variance_decomposition.csv", index=False)

    # figure
    plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False})
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    d = dec.sort_values("shapley_r2", ascending=True)
    colors = {"Dosage transmission (CN\u2192RNA)":"#4a6fa5"}
    bars = ax.barh(d.mechanism, d.shapley_r2,
                   color=[colors.get(m, "#9aa5b1") for m in d.mechanism])
    for b, p in zip(bars, d.pct):
        ax.text(b.get_width()+0.004, b.get_y()+b.get_height()/2, f"{p:.0f}%",
                va="center", fontsize=7.5)
    ax.set_xlabel("Shapley R\u00b2 (share of explained between-gene variance)")
    ax.set_title(f"Dosage transmission dominates protein-responsiveness variance "
                 f"(total R\u00b2 = {full:.2f})", loc="left", fontsize=8.5)
    fig.tight_layout()
    savefig(fig, "mechanism_decomposition.png"); plt.close(fig)

    top = dec.set_index("mechanism")
    record_values("12_mechanism_decomposition", {
        "total_r2": round(full, 4),
        "transmission_pct": round(top.loc["Dosage transmission (CN\u2192RNA)","pct"], 2),
        "transmission_shapley_r2": round(top.loc["Dosage transmission (CN\u2192RNA)","shapley_r2"], 4),
        "n_genes": int(len(gl)),
        "pct_by_mechanism": {k: round(v, 2) for k, v in top["pct"].items()},
    })
    print(f"total_r2={full:.4f}  transmission_pct={top.loc['Dosage transmission (CN\u2192RNA)','pct']:.2f}  n_genes={len(gl)}")

if __name__ == "__main__":
    main()
