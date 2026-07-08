# =============================================================================
# 14_regulatory_decomposition.py  --  Fig 11: (A) regulatory-augmented Shapley
# decomposition of protein responsiveness; (B) determinants of dosage
# transmission itself (methylation contributes, but transmission is mostly
# gene-intrinsic; total R^2 ~ 0.10).
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from itertools import permutations
from sklearn.linear_model import LinearRegression
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import importlib.util
from cnt_io import DIR_TAB, savefig, record_values

# reuse the gene table from stage 12
spec = importlib.util.spec_from_file_location("s12", Path(__file__).parent / "12_mechanism_decomposition.py")
s12 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s12)

STRUCT = s12.STRUCT

def r2_of(sub, groups, y, g):
    f = [c for gn in sub for c in groups[gn]]
    if not f: return 0.0
    return LinearRegression().fit(g[f].values, y).score(g[f].values, y)

def shapley(groups, y, gl):
    gn = list(groups); perms = list(permutations(gn)); shap = {x:0.0 for x in gn}
    for perm in perms:
        seen = []
        for x in perm:
            shap[x] += r2_of(seen+[x], groups, y, gl) - r2_of(seen, groups, y, gl); seen.append(x)
    for x in shap: shap[x] /= len(perms)
    return shap, r2_of(gn, groups, y, gl)

def main():
    gl = s12.build_gene_table()

    groups1 = s12.GROUPS
    shap1, full1 = shapley(groups1, gl["cn_prot_corr"].values, gl)
    d1 = pd.DataFrame([{"mechanism":k,"shapley_r2":v,"pct":100*v/full1,"total_r2":full1}
                       for k,v in shap1.items()]).sort_values("shapley_r2")

    groups2 = {
     "Promoter methylation": ['mean_promoter_meth','meth_rna_corr'],
     "Genomic context (length, density)": ['log_gene_length','gene_density_1mb'],
     "Expression baseline (GTEx)": ['gtex'],
     "Complex membership": ['is_complex_subunit','corum_n_complexes','corum_max_complex_size'],
     "Protein structure": ['VSL2_disorder','PSIPRED_helix','PSIPRED_strand','MMseq2_median','ASAquick_buried'],
    }
    shap2, full2 = shapley(groups2, gl["cn_rna_corr"].values, gl)
    d3 = pd.DataFrame([{"determinant":k,"shapley_r2":v,"pct":100*v/full2,"total_r2":full2}
                       for k,v in shap2.items()]).sort_values("shapley_r2")
    d3.to_csv(DIR_TAB / "transmission_determinants.csv", index=False)
    d1.to_csv(DIR_TAB / "decomposition_with_regulatory.csv", index=False)

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"xtick.labelsize":7,
        "ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    cmap = {"Dosage transmission (CN\u2192RNA)":"#4a6fa5","Regulatory/chromatin":"#8e5aa5",
            "Protein structure/turnover":"#5a9367","Complex buffering (CORUM)":"#c07a2d",
            "Tissue expression (GTEx)":"#8a8a8a"}
    ax = axes[0]
    ax.barh(range(len(d1)), d1.shapley_r2, color=[cmap[m] for m in d1.mechanism])
    ax.set_yticks(range(len(d1))); ax.set_yticklabels([m.replace(" (","\n(") for m in d1.mechanism], fontsize=6.6)
    for i, r in enumerate(d1.itertuples()): ax.text(r.shapley_r2+0.006, i, f"{r.pct:.0f}%", va="center", fontsize=6.5)
    ax.set_xlabel("Shapley R\u00b2"); ax.set_xlim(0, 0.6)
    ax.set_title(f"A  Protein responsiveness (+regulatory)\ntotal R\u00b2={full1:.2f}")
    cmap2 = {"Promoter methylation":"#8e5aa5","Complex membership":"#c07a2d","Protein structure":"#5a9367",
             "Expression baseline (GTEx)":"#8a8a8a","Genomic context (length, density)":"#c74a4a"}
    ax = axes[1]
    ax.barh(range(len(d3)), d3.shapley_r2, color=[cmap2[m] for m in d3.determinant])
    ax.set_yticks(range(len(d3))); ax.set_yticklabels([m.replace(" (","\n(") for m in d3.determinant], fontsize=6.6)
    for i, r in enumerate(d3.itertuples()): ax.text(r.shapley_r2+0.0004, i, f"{r.pct:.0f}%", va="center", fontsize=6.5)
    ax.set_xlabel("Shapley R\u00b2 (of dosage transmission)"); ax.set_xlim(0, max(0.045, d3.shapley_r2.max()*1.3))
    ax.set_title(f"B  What explains dosage transmission itself\ntotal R\u00b2={full2:.2f} (largely gene-intrinsic)")
    fig.suptitle("Regulatory-augmented decomposition: methylation contributes to transmission, "
                 "but transmission is mostly gene-intrinsic", fontsize=8.8, y=1.02)
    fig.tight_layout(); savefig(fig, "regulatory_decomposition.png"); plt.close(fig)

    record_values("14_regulatory_decomposition", {
        "responsiveness_total_r2": round(float(full1),4),
        "transmission_total_r2": round(float(full2),4),
        "responsiveness_pct": {k: round(100*v/full1,2) for k,v in shap1.items()},
        "transmission_determinant_pct": {k: round(100*v/full2,2) for k,v in shap2.items()},
    })
    print(f"responsiveness_total_r2={full1:.4f}  transmission_total_r2={full2:.4f}")

if __name__ == "__main__":
    main()
