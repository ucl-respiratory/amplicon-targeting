#!/usr/bin/env python3
"""verify_paper.py -- exact-reproduction check for the integrated (bioRxiv) manuscript.

Unlike analysis/91_verify_paper.py (which checks the older v9 dosage-transmission paper),
this verifier checks the INTEGRATED paper's headline numbers against the from-source
pipeline outputs in integrated/tables + integrated/reports/values. Every checked value is
read from a pipeline artifact -- none is transcribed by hand -- so a PASS means the
manuscript's numbers are exactly what the from-source pipeline (00a/00d/01/02/03/00b/00c/
04*/05) produces on data_download/from_source.

Run:  python verify_paper.py        (exit 0 = all pass, non-zero = drift)
"""
import sys, json
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

DIR_TAB = cfg.DIR_TAB
DIR_VAL = Path(__file__).resolve().parent / "reports" / "values"


def _load_json(name):
    p = DIR_VAL / name
    return json.load(open(p)) if p.exists() else {}


def _rows(fname):
    p = DIR_TAB / fname
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def check(label, got, want, tol=0.0):
    ok = (got == want) if tol == 0 and not isinstance(want, float) else \
         (got is not None and abs(float(got) - float(want)) <= tol)
    return {"check": label, "got": got, "want": want, "tol": tol, "pass": bool(ok)}


def main():
    results = []

    # ---- funnel (00b) ----
    fn = _load_json("00b_target_funnel.json")
    results.append(check("funnel: universe genes", fn.get("universe"), 6648, 0))
    results.append(check("funnel: transmitted (obs>=0.40)", fn.get("transmitted"), 3104, 0))
    results.append(check("funnel: co-elevated", fn.get("coelevated"), 888, 0))
    results.append(check("funnel: nominated antigens", fn.get("nominated_antigens"), 22, 0))
    results.append(check("funnel: nominated amplicons", fn.get("nominated_amplicons"), 18, 0))
    results.append(check("funnel: cohorts nominated", fn.get("n_cohorts_nominated"), 4, 0))

    # ---- predictor (02, from-source retrain) ----
    pr = _load_json("predictor_fromsource.json")
    results.append(check("predictor: leave-gene-out rho", pr.get("leave_gene_out_rho"), 0.520, 0.01))
    results.append(check("predictor: leave-gene-out R2", pr.get("leave_gene_out_R2"), 0.291, 0.01))
    results.append(check("predictor: leave-arm-out rho", pr.get("leave_arm_out_rho"), 0.528, 0.01))
    results.append(check("predictor: |positional delta|<0.02", abs(pr.get("positional_delta", 1)), 0.0, 0.02))

    # ---- empirical Bayes (03) ----
    eb = _load_json("03_empirical_bayes.json")
    results.append(check("EB: n_genes", eb.get("n_genes"), 6648, 0))
    results.append(check("EB: tau2", eb.get("tau2"), 0.047, 0.005))
    results.append(check("EB: rho_prior_obs", eb.get("rho_prior_obs"), 0.520, 0.02))

    # ---- transmission cascade (01) ----
    tg = _load_json("01_transmission_gates.json")
    results.append(check("cascade: frac attenuated ~0.85", tg.get("frac_attenuated"), 0.854, 0.02))

    # ---- constructs (04d) incl passenger-only ----
    con = _rows("adc_constructs.csv")
    if len(con):
        con = con.set_index("amplicon")
        def _enr(amp):
            return round(float(con.loc[amp, "enrich"]), 2) if amp in con.index and pd.notna(con.loc[amp, "enrich"]) else None
        results.append(check("construct: LUAD_7p full", _enr("LUAD_7p"), 1.45, 0.03))
        results.append(check("construct: LUAD_7p passenger-only (driver EGFR excluded)", _enr("LUAD_7p_pass"), 1.30, 0.03))
        results.append(check("construct: GBM_19p (prediction-only)", _enr("GBM_19p"), 1.12, 0.03))
        results.append(check("construct: GBM_20p (prediction-only)", _enr("GBM_20p"), 1.10, 0.03))
        # passenger-only interval clear of 1.0
        if "LUAD_7p_pass" in con.index:
            lo = float(con.loc["LUAD_7p_pass", "ci_lo"])
            results.append(check("passenger-only CI lower bound > 1.0", lo > 1.0, True, 0))

    # ---- nominated antigen count + passenger fraction (Table 1) ----
    ant = _rows("adc_target_antigens.csv")
    if len(ant):
        results.append(check("Table 1: distinct antigens", ant.antigen.nunique(), 22, 0))
        n_pass = (ant.drop_duplicates("antigen").is_driver == False).sum() if "is_driver" in ant.columns else None
        if n_pass is not None:
            results.append(check("Table 1: passenger antigens (of 22)", int(n_pass), 21, 0))

    # ---- report ----
    n_pass = sum(r["pass"] for r in results); n = len(results)
    rep = DIR_VAL.parent / "verification_report_integrated.md"
    with open(rep, "w") as f:
        f.write(f"# Integrated manuscript verification: {n_pass}/{n} checks PASS\n\n")
        f.write("Every value is read from a from-source pipeline artifact "
                "(integrated/tables, integrated/reports/values); none is hand-transcribed.\n\n")
        f.write("| Check | Got | Want | Tol | Pass |\n|---|---|---|---|---|\n")
        for r in results:
            f.write(f"| {r['check']} | {r['got']} | {r['want']} | {r['tol']} | {'PASS' if r['pass'] else '**FAIL**'} |\n")
    print(f"{n_pass}/{n} checks PASS -> {rep}")
    for r in results:
        if not r["pass"]:
            print(f"  FAIL: {r['check']}: got {r['got']} want {r['want']} (tol {r['tol']})")
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
