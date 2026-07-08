# =============================================================================
# 90_values_manifest.py  --  aggregate every recorded value into one manifest.
# Each analysis stage writes OUT/reports/values/<stage>.json via cnt_io.record_values.
# This flattens all of them into values_manifest.json (nested) and
# values_manifest.csv (one row per leaf value: stage, key path, value), the
# single source mapping each in-text number / Table 1 entry to its producing stage.
# =============================================================================
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from cnt_io import DIR_REP

def flatten(prefix, obj, rows):
    if isinstance(obj, dict):
        for k, v in obj.items(): flatten(f"{prefix}.{k}" if prefix else k, v, rows)
    elif isinstance(obj, list):
        rows.append((prefix, json.dumps(obj)))
    else:
        rows.append((prefix, obj))

def main():
    vdir = DIR_REP / "values"
    manifest = {}
    rows = []
    for jf in sorted(vdir.glob("*.json")):
        stage = jf.stem
        d = json.load(open(jf))
        manifest[stage] = d
        st_rows = []
        flatten("", d, st_rows)
        for key, val in st_rows:
            rows.append({"stage": stage, "key": key, "value": val})
    (DIR_REP / "values_manifest.json").write_text(json.dumps(manifest, indent=2))
    df = pd.DataFrame(rows)
    df.to_csv(DIR_REP / "values_manifest.csv", index=False)
    print(f"stages={len(manifest)}  leaf_values={len(df)}")
    print(f"written: values_manifest.json + values_manifest.csv")

if __name__ == "__main__":
    main()
