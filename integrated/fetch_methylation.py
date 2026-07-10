#!/usr/bin/env python3
# Download CPTAC-3 EPIC-v1 methylation betas for matched cases via GDC's BULK
# /data endpoint (POST multiple UUIDs -> one tar.gz), in resumable batches.
# Much faster than per-file GETs. No ChAMP/minfi.
import json, io, tarfile, time, urllib.request
from pathlib import Path
HERE = Path(__file__).resolve().parent
CACHE = HERE.parent / "data_download/from_source/data/meth_cache"
CACHE.mkdir(parents=True, exist_ok=True)
man = json.load(open(HERE/"reports/meth_manifest.json"))

# which are already present (resumable)
def present(fname): 
    p = CACHE/fname
    return p.exists() and p.stat().st_size > 1000
todo = [(cid, r) for cid, r in man.items() if not present(r["file_name"])]
print(f"[meth] {len(todo)}/{len(man)} still to fetch (bulk POST)", flush=True)

def bulk(batch):
    ids = [r["file_id"] for _, r in batch]
    payload = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/data",
        data=payload, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        data = r.read()
    # single-file responses are raw text; multi-file are tar.gz
    if len(ids) == 1:
        (CACHE/batch[0][1]["file_name"]).write_bytes(data); return 1
    n = 0
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile(): continue
            base = Path(m.name).name
            if base == "MANIFEST.txt": continue
            fo = tf.extractfile(m)
            if fo is None: continue
            (CACHE/base).write_bytes(fo.read()); n += 1
    return n

B = 20
done = 0
for i in range(0, len(todo), B):
    batch = todo[i:i+B]
    for attempt in range(4):
        try:
            got = bulk(batch); done += got
            print(f"[meth] batch {i//B+1}: +{got} (total new {done}/{len(todo)})", flush=True)
            break
        except Exception as e:
            if attempt == 3:
                print(f"[meth] batch {i//B+1} FAILED: {str(e)[:120]}", flush=True)
            else:
                time.sleep(5*(attempt+1))
print(f"[meth] bulk download complete: {done} new files", flush=True)
