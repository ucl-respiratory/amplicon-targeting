#!/usr/bin/env python3
# =============================================================================
# 03_empirical_bayes.py  --  Combining measured and predicted transmissibility.
# -----------------------------------------------------------------------------
# The measured co-elevation catalogue is precise where proteomics is deep and
# noisy where it is thin; the gene-property predictor is available everywhere but
# explains only ~1/3 of per-gene variance. An empirical-Bayes posterior blends
# them: it shrinks the noisy per-gene measurement toward the predictor prior in
# inverse proportion to the measurement's precision (binomial, set by the number
# of amplified cases). Where proteomics is rich the data dominate; where it is
# thin or absent the prior carries the nomination.
#
#   posterior_g = (mu_g / tau^2 + y_g / sigma_g^2) / (1/tau^2 + 1/sigma_g^2)
#     mu_g   = predicted (prior)              tau^2   = between-gene prior variance
#     y_g    = observed (measured)            sigma_g^2 = p(1-p)/n  (binomial)
#
# A thin-cohort recovery test quantifies the gain: genes with >=150 amplified
# cases are treated as reference truth; binomial noise is injected at a sweep of
# cohort sizes; raw-measurement, prior-only and posterior are scored against the
# reference by Spearman rho and RMSE. The prior-only score is the ceiling for a
# zero-proteomics context.
#
# Inputs: transmissibility atlas (measured + predicted). Runs now.
# Outputs: figures/fig4_empirical_bayes.png, tables/posterior_transmissibility.csv,
#          reports/eb_key_numbers.json
# =============================================================================
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cnt_shared as ish
from cnt_shared import cfg
ish.apply_style(sizes=(9, 8, 7))

def binom_var(p, n):
    p = np.clip(p, 1e-3, 1 - 1e-3)
    return p * (1 - p) / np.maximum(n, 1)

def posterior(mu, y, n, tau2):
    """Precision-weighted posterior with binomial measurement variance."""
    s2 = binom_var(y, n)
    w_prior = 1.0 / tau2
    w_meas = 1.0 / s2
    post = (mu * w_prior + y * w_meas) / (w_prior + w_meas)
    shrink = w_prior / (w_prior + w_meas)   # weight on the prior (0=data,1=prior)
    return post, shrink

def main():
    a = ish.load_atlas().dropna(
        subset=["observed_transmissibility", "predicted_transmissibility_oof",
                "n_amplified_cases"]).copy()
    mu = a.predicted_transmissibility_oof.values     # prior mean
    y = a.observed_transmissibility.values           # measurement
    n = a.n_amplified_cases.values.astype(float)

    # tau^2 = between-gene residual variance of measurement about prior, net of
    # mean measurement noise (method-of-moments empirical-Bayes).
    mean_meas_var = float(np.mean(binom_var(y, n)))
    tau2 = max(float(np.var(y - mu) - mean_meas_var), 1e-4)

    post, shrink = posterior(mu, y, n, tau2)
    a["posterior_transmissibility"] = post
    a["prior_weight"] = shrink
    a.sort_values("posterior_transmissibility", ascending=False)[
        ["gene", "observed_transmissibility", "predicted_transmissibility_oof",
         "n_amplified_cases", "posterior_transmissibility", "prior_weight"]
    ].to_csv(cfg.DIR_TAB / "posterior_transmissibility.csv", index=False)

    rho_prior_obs = float(spearmanr(y, mu).correlation)

    # ---- thin-cohort recovery test -----------------------------------------
    rng = np.random.default_rng(cfg.SEED)
    ref = a[a.n_amplified_cases >= cfg.EB_REFERENCE_MIN_CASES].copy()
    truth = ref.observed_transmissibility.values
    prior_ref = ref.predicted_transmissibility_oof.values
    rho_prior_floor = float(spearmanr(truth, prior_ref).correlation)

    grid = cfg.EB_COHORT_GRID
    curve = []
    for nc in grid:
        if nc == 0:
            curve.append({"n": 0, "rho_raw": np.nan, "rho_prior": rho_prior_floor,
                          "rho_post": rho_prior_floor, "rmse_raw": np.nan,
                          "rmse_prior": float(np.sqrt(np.mean((prior_ref - truth)**2))),
                          "rmse_post": float(np.sqrt(np.mean((prior_ref - truth)**2)))}); continue
        rr, rpo, er, ep = [], [], [], []
        for _ in range(cfg.EB_RESAMPLES):
            noisy = np.clip(truth + rng.normal(0, np.sqrt(binom_var(truth, nc)), len(truth)), 0, 1)
            pst, _ = posterior(prior_ref, noisy, np.full(len(truth), nc), tau2)
            rr.append(spearmanr(noisy, truth).correlation)
            rpo.append(spearmanr(pst, truth).correlation)
            er.append(np.sqrt(np.mean((noisy - truth)**2)))
            ep.append(np.sqrt(np.mean((pst - truth)**2)))
        curve.append({"n": nc, "rho_raw": float(np.mean(rr)),
                      "rho_prior": rho_prior_floor, "rho_post": float(np.mean(rpo)),
                      "rmse_raw": float(np.mean(er)),
                      "rmse_prior": float(np.sqrt(np.mean((prior_ref - truth)**2))),
                      "rmse_post": float(np.mean(ep))})
    cv = pd.DataFrame(curve)
    cv.to_csv(cfg.DIR_TAB / "eb_recovery_curve.csv", index=False)

    # crossover: smallest n where raw measurement beats the prior floor
    beat = cv[(cv.n > 0) & (cv.rho_raw > cv.rho_prior)]
    crossover_n = int(beat.n.min()) if len(beat) else None
    posterior_dominates = bool((cv[cv.n > 0].rmse_post <=
                                cv[cv.n > 0][["rmse_raw", "rmse_prior"]].min(axis=1) + 1e-9).all())

    km = {"n_genes": int(len(a)), "tau2": tau2, "mean_meas_var": mean_meas_var,
          "rho_prior_obs": rho_prior_obs, "rho_prior_floor": rho_prior_floor,
          "ref_genes": int(len(ref)), "crossover_n": crossover_n,
          "posterior_dominates_rmse": posterior_dominates,
          "shrink_median": float(np.median(shrink)), "shrink_min": float(shrink.min()),
          "shrink_max": float(shrink.max())}
    (cfg.DIR_REP / "eb_key_numbers.json").write_text(json.dumps(km, indent=2))

    # ---- Figure 4: three panels --------------------------------------------
    fig = plt.figure(figsize=(13, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 1], wspace=0.32)
    C_RAW, C_PRIOR, C_POST = "#c0392b", "#2c6fbb", "#6c5b9c"

    # panel a: shrinkage toward prior vs number of amplified cases
    axa = fig.add_subplot(gs[0, 0])
    axa.scatter(a.n_amplified_cases, shrink, s=6, alpha=0.3, color="#666",
                rasterized=True, edgecolors="none")
    axa.set_xlabel("amplified cases  (measurement precision)")
    axa.set_ylabel("weight on prior  (1 = prior, 0 = data)")
    axa.set_title("Fewer cases → lean on the prior;\nmore cases → trust the measurement",
                  fontsize=9, loc="left")

    # panel b: recovery curve (rho vs cohort size)
    axb = fig.add_subplot(gs[0, 1])
    d = cv[cv.n > 0]
    axb.plot(d.n, d.rho_raw, "o-", color=C_RAW, ms=4, lw=1.4)
    axb.plot(d.n, d.rho_post, "s-", color=C_POST, ms=4, lw=1.4)
    axb.axhline(rho_prior_floor, ls="--", color=C_PRIOR, lw=1.4)
    axb.set_xscale("log")
    axb.set_xlabel("cohort size  (amplified cases)")
    axb.set_ylabel("rank recovery of truth  (ρ)")
    axb.set_title("Posterior never below measurement or prior", fontsize=9, loc="left")
    # direct labels at right end
    xr = d.n.max()
    axb.annotate("measured only", xy=(xr, d.rho_raw.iloc[-1]), xytext=(6, -2),
                 textcoords="offset points", color=C_RAW, fontsize=7.5, va="center")
    axb.annotate("posterior", xy=(xr, d.rho_post.iloc[-1]), xytext=(6, 4),
                 textcoords="offset points", color=C_POST, fontsize=7.5, va="center")
    axb.annotate(f"prior only (ρ={rho_prior_floor:.2f})", xy=(d.n.iloc[0], rho_prior_floor),
                 xytext=(0, -12), textcoords="offset points", color=C_PRIOR, fontsize=7.5)

    # panel c: RMSE reduction of posterior vs measurement
    axc = fig.add_subplot(gs[0, 2])
    red = 100 * (d.rmse_raw - d.rmse_post) / d.rmse_raw
    axc.plot(d.n, red, "o-", color=C_POST, ms=4, lw=1.4)
    axc.axhline(0, color="#999", lw=0.8)
    axc.set_xscale("log")
    axc.set_xlabel("cohort size  (amplified cases)")
    axc.set_ylabel("RMSE reduction vs measurement  (%)")
    axc.set_title("Largest gain in the thin-cohort regime", fontsize=9, loc="left")
    for ax, L in zip([axa, axb, axc], "abc"):
        ish.panel_letter(ax, L)

    fig.savefig(cfg.DIR_FIG / "fig4_empirical_bayes.png", dpi=200, bbox_inches="tight")

    ish.record("03_empirical_bayes", km)
    print(json.dumps(km, indent=2))
    top = a.sort_values("posterior_transmissibility", ascending=False).head(6)
    print(top[["gene", "posterior_transmissibility"]].to_string(index=False))

if __name__ == "__main__":
    main()
