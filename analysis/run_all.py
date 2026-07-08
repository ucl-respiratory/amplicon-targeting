#!/usr/bin/env python3
# =============================================================================
# run_all.py  --  Master runner for the CN-targeting ANALYSIS pipeline
# -----------------------------------------------------------------------------
# Runs every numbered analysis stage in order, each of which regenerates one or
# more manuscript figures and writes a per-stage values JSON. Stage 90 then
# aggregates the values manifest and stage 91 verifies against the paper.
#
# USAGE
#   export CNT_DATA=/path/to/data_download/from_source   # data_download output
#   python run_all.py                 # all stages
#   python run_all.py --stages 10:40  # a contiguous numeric subset
#   python run_all.py --skip-cellxgene  # skip the census stage (Fig 6)
#
# Each stage is an importable module with a main(); failures are logged and (for
# the optional census stage) tolerated.
# =============================================================================
import sys, time, importlib.util, argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent

STAGES = [
    "10_transmission_reconciliation.py",
    "11_per_tissue_attenuation.py",
    "12_mechanism_decomposition.py",
    "13_coelevation_by_transmission.py",
    "14_regulatory_decomposition.py",
    "15_atac_transmission.py",
    "16_statistical_robustness.py",
    "20_tissue_amplicons.py",
    "21_dependency_validation.py",
    "22_target_funnel.py",
    "23_multispecific_safety.py",
    "24_offtarget_safety.py",
    "30_normal_singlecell_threshold.py",
    "31_tumour_samecell.py",
    "90_values_manifest.py",
    "91_verify_paper.py",
]
CENSUS_STAGE = "31_tumour_samecell.py"

def load(mod_file):
    spec = importlib.util.spec_from_file_location(mod_file[:-3], HERE / mod_file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default=None, help="numeric range lo:hi (inclusive)")
    ap.add_argument("--skip-cellxgene", action="store_true")
    args = ap.parse_args()

    sel = STAGES
    if args.stages:
        lo, hi = (int(x) for x in args.stages.split(":"))
        sel = [s for s in STAGES if lo <= int(s[:2]) <= hi]
    if args.skip_cellxgene:
        sel = [s for s in sel if s != CENSUS_STAGE]

    print("="*64); print(" CN-targeting analysis pipeline"); print(" stages:",
          ", ".join(s[:-3] for s in sel)); print("="*64)
    log = []
    for s in sel:
        print(f"\n>>> {s}")
        t = time.time()
        try:
            m = load(s)
            (m.main() if hasattr(m, "main") else None)
            status = "ok"
        except Exception as e:
            status = f"error: {e}"
            print(f"  !! {status}")
            if s != CENSUS_STAGE:   # census stage tolerated
                log.append((s, status)); break
        log.append((s, status)); print(f"  [{status}] {time.time()-t:.1f}s")
    print("\n"+"="*64)
    for s, st in log: print(f"  {s:36s} {st}")
    print("="*64)
    if any(st.startswith("error") for _, st in log if _ != CENSUS_STAGE):
        sys.exit(1)

if __name__ == "__main__":
    main()
