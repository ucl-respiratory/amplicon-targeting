# =============================================================================
# 11_per_tissue_attenuation.py  --  Fig 8: attenuation holds in every tumour
# type (protein tracks CN ~half as strongly as mRNA). Uses cn_adjusted here
# (per-tumour ploidy-adjusted), matching the manuscript figure.
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io
from cnt_io import DIR_TAB, savefig, record_values

def main():
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","rna","prot"]).copy()
    rows = []
    for tc, g in d.groupby("tumor_code"):
        rna_c, prot_c = [], []
        for gene, gg in g.groupby("gene"):
            if len(gg) < 20: continue
            if gg.cn_adjusted.std() > 0 and gg.rna.std() > 0 and gg.prot.std() > 0:
                rna_c.append(pearsonr(gg.cn_adjusted, gg.rna)[0])
                prot_c.append(pearsonr(gg.cn_adjusted, gg.prot)[0])
        rows.append({"tumor_code":tc,"n_genes":len(rna_c),
                     "mean_cn_rna":np.nanmean(rna_c),"mean_cn_prot":np.nanmean(prot_c)})
    pt = pd.DataFrame(rows); pt["attenuation"] = pt.mean_cn_rna - pt.mean_cn_prot
    pt = pt.sort_values("mean_cn_rna")
    pt.to_csv(DIR_TAB / "per_tissue_cn_rna_prot.csv", index=False)

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    y = np.arange(len(pt)); h = 0.38
    ax.barh(y+h/2, pt.mean_cn_rna, h, color="#4a6fa5", label="CN\u2192mRNA")
    ax.barh(y-h/2, pt.mean_cn_prot, h, color="#5a9367", label="CN\u2192protein")
    for i, r in enumerate(pt.itertuples()):
        ax.text(r.mean_cn_rna+0.005, i+h/2, f"{r.mean_cn_rna:.2f}", va="center", fontsize=6.3)
        ax.text(r.mean_cn_prot+0.005, i-h/2, f"{r.mean_cn_prot:.2f}", va="center", fontsize=6.3)
    ax.set_yticks(y); ax.set_yticklabels(pt.tumor_code)
    ax.set_xlabel("Mean per-gene correlation with copy number")
    ax.set_title("Attenuation holds in every tumour type\n(protein tracks CN ~half as strongly as mRNA)")
    ax.legend(frameon=False, fontsize=7, loc="lower right"); ax.set_xlim(0, 0.48)
    fig.tight_layout(); savefig(fig, "per_tissue_attenuation.png"); plt.close(fig)

    record_values("11_per_tissue_attenuation", {
        "per_type": {r.tumor_code: {"mean_cn_rna": round(r.mean_cn_rna,4),
                     "mean_cn_prot": round(r.mean_cn_prot,4),
                     "attenuation": round(r.attenuation,4), "n_genes": int(r.n_genes)}
                     for r in pt.itertuples()}})
    print(pt.to_string(index=False))

if __name__ == "__main__":
    main()
