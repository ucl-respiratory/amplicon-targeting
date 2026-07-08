# =============================================================================
# 23_multispecific_safety.py  --  Fig 4: multispecific single-molecule reframe.
# (A) UniProt topology filter (31/62 accessible); (B) OR-gate (cocktail) vs
# AND-gate (multispecific) normal-organ Medium+ burden from HPA IHC; (C) the
# AND-gate preserves the tumour signal (all antigens elevated | amplicon).
# Writes multispecific_andgate_sets.csv (the multi-antigen accessible sets).
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io, _depmap
from cnt_io import DIR_TAB, savefig, record_values, cfg

VULN = {"Heart":["heart muscle"],"Liver":["liver"],"Lung":["lung","bronchus"],"Kidney":["kidney"],
 "Brain":["cerebral cortex","cerebellum","hippocampus","caudate","hypothalamus"],
 "GI_tract":["colon","duodenum","small intestine","stomach","stomach 1","stomach 2","rectum","esophagus","oral mucosa"],
 "Bone_marrow_blood":["bone marrow","lymph node","spleen"],"Pancreas":["pancreas"],
 "Skin":["skin","skin 1","skin 2"]}
LVL = {"Not detected":0,"Low":1,"Medium":2,"High":3}
SET_ORDER = ["LUAD_1q","LUAD_7p","LSCC_1q","PDA_1q","LSCC_20q","UCEC_1q","LUAD_5p"]

def main():
    tp = cnt_io.load_topology().set_index("gene")
    co = _depmap.cotarget_dependency_annotated()
    acc = set(tp[tp.surface_accessible].index)
    co = co.copy(); co["accessible"] = co.gene.isin(acc)
    elig = co[(co.fdr<0.1)&(co.dependency_class=="dispensable_passenger")&(co.accessible)]
    arm = elig.groupby(["tumor_code","arm"]).agg(
        n=("gene","nunique"), genes=("gene", lambda s: ";".join(sorted(set(s))))).reset_index()
    sets_multi = arm[arm.n>=2].sort_values("n", ascending=False)

    main_t = cnt_io.load_str_omic()[["gene","caseid","tumor_code","cn_adjusted","prot.rel.tissue"]]
    ihc = cnt_io.load_hpa_normal(); ihc["lvl"] = ihc.Level.map(LVL)

    def organ_level(gene, tissues):
        s = ihc[(ihc["Gene name"]==gene)&(ihc.Tissue.isin(tissues))]["lvl"].dropna()
        return s.max() if len(s) else np.nan

    def tumor_gate(genes, tc):
        d = main_t[(main_t.tumor_code==tc)&(main_t.gene.isin(genes))]
        pp = d.pivot_table(index="caseid",columns="gene",values="prot.rel.tissue",aggfunc="first")
        pc = d.pivot_table(index="caseid",columns="gene",values="cn_adjusted",aggfunc="first")
        gp = [g for g in genes if g in pp.columns and g in pc.columns]
        common = pp[gp].dropna().index.intersection(pc[gp].dropna().index)
        pp = pp.loc[common,gp]; pc = pc.loc[common,gp]
        amp = pc.mean(axis=1)>=cfg.AMP_THRESHOLD; high = pp>cfg.REL_TISSUE_HI
        return {"n_amp":int(amp.sum()),
                "AND_amp":float(high.all(axis=1)[amp].mean()) if amp.sum() else np.nan,
                "AND_noamp":float(high.all(axis=1)[~amp].mean()) if (~amp).sum() else np.nan,
                "OR_amp":float(high.any(axis=1)[amp].mean()) if amp.sum() else np.nan,"genes":gp}

    sets_final = {f"{r.tumor_code}_{r.arm}": r.genes.split(";") for r in sets_multi.itertuples()}
    rec = []
    for sn, genes in sets_final.items():
        tc = sn.split("_")[0]; tg = tumor_gate(genes, tc)
        ormed=andmed=orhigh=andhigh=0
        for organ, tis in VULN.items():
            levels=[organ_level(g,tis) for g in genes]; levels=[l for l in levels if not np.isnan(l)]
            if not levels: continue
            ormed+=int(max(levels)>=2); andmed+=int(min(levels)>=2); orhigh+=int(max(levels)>=3); andhigh+=int(min(levels)>=3)
        rec.append({"set":sn,"genes":";".join(tg["genes"]),"n_members":len(genes),
                    "OR_med_organs":ormed,"AND_med_organs":andmed,"OR_high_organs":orhigh,"AND_high_organs":andhigh,
                    "n_amp":tg["n_amp"],
                    "tumour_AND_amp":tg["AND_amp"],"tumour_AND_noamp":tg["AND_noamp"],"tumour_OR_amp":tg["OR_amp"]})
    GA = pd.DataFrame(rec)
    GA.to_csv(DIR_TAB / "multispecific_andgate_sets.csv", index=False)
    order = [s for s in SET_ORDER if s in GA.set.values]
    GAi = GA.set_index("set").loc[order]

    plt.rcParams.update({"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":6.8,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1,3,figsize=(15,4.7))
    ax=axes[0]
    keys=["surface_TM_ectodomain","surface_GPI","membrane_associated_cytoplasmic_face","intracellular_membrane","intracellular"]
    labels=["surface\nTM+ectodomain","surface\nGPI","membrane assoc.\n(cytoplasmic face)","intracellular\nmembrane","intracellular"]
    vc=tp.surface_class.value_counts(); vals=[int(vc.get(k,0)) for k in keys]; cols=["#2c7a3f","#7cae4a","#d08a3a","#b04a3a","#8a8a8a"]
    ax.barh(range(len(keys)),vals,color=cols); ax.set_yticks(range(len(keys))); ax.set_yticklabels(labels); ax.invert_yaxis()
    for i,v in enumerate(vals): ax.text(v+0.3,i,str(v),va="center",fontsize=8)
    n_acc=int(tp.surface_accessible.sum()); n_tot=len(tp)
    ax.set_xlabel(f"Nominated antigens (of {n_tot})")
    ax.set_title(f"A  Topology-aware filter (UniProt):\n{n_acc}/{n_tot} have accessible ectodomain"); ax.axhline(1.5,color="k",ls="--",lw=0.8)
    ax=axes[1]; x=np.arange(len(order)); w=0.38
    ax.bar(x-w/2,GAi.OR_med_organs,w,label="OR gate (ADC cocktail)",color="#b04a3a")
    ax.bar(x+w/2,GAi.AND_med_organs,w,label="AND gate (multispecific)",color="#2c7a3f")
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_","\n") for s in order]); ax.set_ylabel("Vulnerable organs, Medium+ burden")
    ax.set_title("B  Multispecific (AND) shrinks\nnormal-tissue burden vs cocktail (OR)"); ax.legend(frameon=False,loc="upper right")
    for i in range(len(order)):
        ax.text(i-w/2,GAi.OR_med_organs.iloc[i]+0.1,int(GAi.OR_med_organs.iloc[i]),ha="center",fontsize=6.5)
        ax.text(i+w/2,GAi.AND_med_organs.iloc[i]+0.1,int(GAi.AND_med_organs.iloc[i]),ha="center",fontsize=6.5)
    ax=axes[2]
    ax.bar(x-w/2,GAi.tumour_AND_amp,w,label="amplified tumours",color="#2c6fb0")
    ax.bar(x+w/2,GAi.tumour_AND_noamp,w,label="non-amplified",color="#b9c9de")
    ax.set_xticks(x); ax.set_xticklabels([s.replace("_","\n") for s in order]); ax.set_ylabel("Tumours co-expressing ALL members")
    ax.set_title("C  AND-gate preserves tumour signal\n(all antigens elevated | amplicon)"); ax.legend(frameon=False,loc="upper right"); ax.set_ylim(0,1)
    for i in range(len(order)):
        ax.text(i-w/2,GAi.tumour_AND_amp.iloc[i]+0.02,f"{GAi.tumour_AND_amp.iloc[i]:.2f}",ha="center",fontsize=6.3)
        ax.text(i+w/2,GAi.tumour_AND_noamp.iloc[i]+0.02,f"{GAi.tumour_AND_noamp.iloc[i]:.2f}",ha="center",fontsize=6.3)
    fig.suptitle("Multispecific single-molecule reframing: topology-filtered accessible antigen sets and AND-gate tumour selectivity",fontsize=9.5,y=1.02)
    fig.tight_layout(); savefig(fig, "multispecific_andgate.png"); plt.close(fig)

    record_values("23_multispecific_safety", {
        "n_multi_antigen_sets": int(len(sets_multi)),
        "n_surface_accessible": int(tp.surface_accessible.sum()),
        "sets": {r.set: {"n_members":int(r.n_members),"OR_med_organs":int(r.OR_med_organs),
                 "AND_med_organs":int(r.AND_med_organs),"AND_high_organs":int(r.AND_high_organs),
                 "tumour_AND_amp":round(float(r.tumour_AND_amp),4) if pd.notna(r.tumour_AND_amp) else None,
                 "tumour_AND_noamp":round(float(r.tumour_AND_noamp),4) if pd.notna(r.tumour_AND_noamp) else None}
                 for r in GAi.reset_index().itertuples()},
    })
    print(GAi[["n_members","OR_med_organs","AND_med_organs","n_amp","tumour_AND_amp","tumour_AND_noamp"]].to_string())

if __name__ == "__main__":
    main()
