# =============================================================================
# 16_statistical_robustness.py  --  Fig 13: robustness of the transmission
# result. (A) arm-grouped vs random CV (neighbour-gene leakage test);
# (B) split-half reliability (~0.74) + disattenuated R^2; (C) clinical
# discrimination AUC (weak ranker -> co-targets come from direct CPTAC+DepMap).
# =============================================================================
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold, KFold
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.ensemble import GradientBoostingClassifier
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_io, _features, _amplicon
from cnt_io import DIR_TAB, savefig, record_values, cfg

TYPES = ["LUAD","LSCC","CCRCC","GBM","UCEC"]
CPTAC_FEATS = ['mean_promoter_meth','meth_rna_corr','log_gene_length','gene_density_1mb',
               'is_complex_subunit','corum_n_complexes','corum_max_complex_size']

def _pg_r(sub, xcol, ycol, min_n=20):
    s = sub[["gene",xcol,ycol]].copy()
    s["xy"]=s[xcol]*s[ycol]; s["x2"]=s[xcol]**2; s["y2"]=s[ycol]**2
    g = s.groupby("gene").agg(mx=(xcol,"mean"),my=(ycol,"mean"),mx2=("x2","mean"),
        my2=("y2","mean"),mxy=("xy","mean"),n=(xcol,"size"))
    vx=g.mx2-g.mx**2; vy=g.my2-g.my**2; cov=g.mxy-g.mx*g.my
    return (cov/np.sqrt(vx*vy)).where(g.n>=min_n).replace([np.inf,-np.inf],np.nan).dropna()

def _cv_r2(X, y, groups=None):
    it = (KFold(5, shuffle=True, random_state=cfg.SEED_BOOTSTRAP).split(X) if groups is None
          else GroupKFold(5).split(X, y, groups))
    ps = np.zeros(len(y))
    for tr, te in it: ps[te] = LinearRegression().fit(X[tr], y[tr]).predict(X[te])
    return r2_score(y, ps)

def main():
    df = cnt_io.load_str_omic()
    d = df.dropna(subset=["cn_adjusted","rna"]).copy()
    reg = _features.regulatory_gene_features()
    corum = _features.corum_gene_features()
    accdf = cnt_io.load_atac()
    if accdf.index.name is None: accdf.index.name = "gene"
    gene_arm = df[["gene","arm"]].drop_duplicates("gene").set_index("gene")["arm"]

    def build(ct):
        tr = _pg_r(d[d.tumor_code==ct],"cn_adjusted","rna").rename("transmission").reset_index()
        tr = tr.merge(reg[["gene","mean_promoter_meth","meth_rna_corr","log_gene_length","gene_density_1mb"]],on="gene",how="left")
        tr = tr.merge(corum,on="gene",how="left")
        for c in ["is_complex_subunit","corum_n_complexes","corum_max_complex_size"]: tr[c]=tr[c].fillna(0)
        tr = tr.merge(accdf[ct].rename("atac").reset_index().rename(columns={"index":"gene"}),on="gene",how="left")
        return tr

    # split-half reliability
    rng = np.random.default_rng(cfg.SEED_BOOTSTRAP)
    rel_rows = []
    for ct in TYPES:
        g = d[d.tumor_code==ct]; cases = g.caseid.unique().copy(); rng.shuffle(cases); h=len(cases)//2
        tA=_pg_r(g[g.caseid.isin(cases[:h])],"cn_adjusted","rna")
        tB=_pg_r(g[g.caseid.isin(cases[h:])],"cn_adjusted","rna")
        common=tA.index.intersection(tB.index)
        r_half=np.corrcoef(tA[common],tB[common])[0,1]; r_full=2*r_half/(1+r_half)
        rel_rows.append({"tumor_code":ct,"n_cases":len(cases),"reliability_fulln":r_full,"r_halfsplit":r_half})
    rel = pd.DataFrame(rel_rows)

    rows = []
    for ct in TYPES:
        t = build(ct); t["arm"]=t.gene.map(gene_arm)
        t = t.dropna(subset=["transmission"]+CPTAC_FEATS+["atac","arm"])
        y=t.transmission.values; X=t[CPTAC_FEATS+["atac"]].values
        r2_rand=_cv_r2(X,y); r2_arm=_cv_r2(X,y,t.arm.values)
        relf=rel.set_index("tumor_code").loc[ct,"reliability_fulln"]
        rows.append({"tumor_code":ct,"n_genes":len(t),"n_arms":t.arm.nunique(),
                     "R2_random_CV":r2_rand,"R2_armgrouped_CV":r2_arm,"R2_disattenuated":r2_rand/relf})
    res = pd.DataFrame(rows)
    res.to_csv(DIR_TAB / "transmission_robustness.csv", index=False)

    co = _amplicon.amplicon_coelevation()
    co["reliable"]=((co.fdr<0.1)&(co.p_high_amp>=0.5)).astype(int)
    clin_rows=[]
    for ct in ["LUAD","LSCC","CCRCC","UCEC"]:
        c=co[co.tumor_code==ct].copy()
        c=c.merge(reg[["gene","mean_promoter_meth","log_gene_length","gene_density_1mb"]],on="gene",how="left")
        c=c.merge(corum,on="gene",how="left")
        for cc in ["is_complex_subunit","corum_max_complex_size"]: c[cc]=c[cc].fillna(0)
        c=c.merge(accdf[ct].rename("atac").reset_index().rename(columns={"index":"gene"}),on="gene",how="left")
        c["arm"]=c.gene.map(gene_arm)
        feats=["atac","mean_promoter_meth","log_gene_length","gene_density_1mb","is_complex_subunit","corum_max_complex_size"]
        c=c.dropna(subset=feats+["reliable","arm"])
        if c.reliable.nunique()<2 or len(c)<50: continue
        y=c.reliable.values; X=c[feats].values; nsp=min(5,c.arm.nunique()); oof=np.zeros(len(y))
        for tr,te in GroupKFold(nsp).split(X,y,c.arm.values):
            oof[te]=GradientBoostingClassifier(n_estimators=150,max_depth=3,random_state=cfg.SEED_BOOTSTRAP).fit(X[tr],y[tr]).predict_proba(X[te])[:,1]
        auc=roc_auc_score(y,oof); order=np.argsort(-oof); topk=max(10,len(y)//10)
        clin_rows.append({"tumor_code":ct,"n":len(c),"n_arms":c.arm.nunique(),"base_rate":y.mean(),
                          "AUC":auc,"precision_top_decile":y[order[:topk]].mean(),"lift":y[order[:topk]].mean()/y.mean()})
    clin = pd.DataFrame(clin_rows)

    plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"legend.fontsize":7,
        "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,
        "axes.titlelocation":"left"})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ct_o=res.tumor_code.tolist(); x=np.arange(len(ct_o)); w=0.38
    ax=axes[0]
    ax.bar(x-w/2,res.R2_random_CV,w,color="#8aa1b4",label="random CV")
    ax.bar(x+w/2,res.R2_armgrouped_CV,w,color="#4a6fa5",label="arm-grouped CV")
    ax.set_xticks(x); ax.set_xticklabels(ct_o); ax.set_ylabel("R\u00b2 (transmission)")
    ax.set_title("A  Neighbour-gene leakage test\narm-grouped \u2248 random \u2192 no inflation"); ax.legend(frameon=False,fontsize=6.8)
    ax=axes[1]
    ax.bar(x-w/2,res.R2_random_CV,w,color="#8aa1b4",label="observed R\u00b2")
    ax.bar(x+w/2,res.R2_disattenuated,w,color="#c07a2d",label="disattenuated (\u00f7reliability)")
    ax.set_xticks(x); ax.set_xticklabels(ct_o); ax.set_ylabel("R\u00b2 (transmission)")
    ax.set_title("B  Noise from ~100 tumours/type\ndeflates observed R\u00b2 (reliability\u22480.74)"); ax.legend(frameon=False,fontsize=6.8)
    ax=axes[2]
    xc=np.arange(len(clin))
    ax.bar(xc,clin.AUC,0.6,color=["#c9a97a" if a>0.55 else "#c7ccd1" for a in clin.AUC])
    ax.axhline(0.5,color="#888",lw=0.8,ls="--")
    ax.set_xticks(xc); ax.set_xticklabels(clin.tumor_code); ax.set_ylim(0.45,0.65)
    ax.set_ylabel("Arm-grouped AUC (predict reliable co-elevation)")
    ax.set_title("C  Feature model is a weak ranker\n(co-targets come from direct CPTAC+DepMap)")
    fig.suptitle("Statistical robustness of the transmission result: no leakage inflation, "
                 "noise-deflated R\u00b2, and why low R\u00b2 does not limit co-target nomination", fontsize=8.4, y=1.03)
    fig.tight_layout(); savefig(fig, "transmission_statistical_robustness.png"); plt.close(fig)

    record_values("16_statistical_robustness", {
        "mean_reliability": round(float(rel.reliability_fulln.mean()),4),
        "mean_R2_random": round(float(res.R2_random_CV.mean()),4),
        "mean_R2_armgrouped": round(float(res.R2_armgrouped_CV.mean()),4),
        "mean_R2_disattenuated": round(float(res.R2_disattenuated.mean()),4),
        "clinical_AUC": {r.tumor_code: round(r.AUC,4) for r in clin.itertuples()},
    })
    print("reliability mean:", round(rel.reliability_fulln.mean(),4))
    print(res[["tumor_code","R2_random_CV","R2_armgrouped_CV","R2_disattenuated"]].to_string(index=False))
    print("AUC:", {r.tumor_code: round(r.AUC,3) for r in clin.itertuples()})

if __name__ == "__main__":
    main()
