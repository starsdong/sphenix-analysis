#!/usr/bin/env python3
"""Make per-job DST lists for (a fraction of) the run3pp DST_TRKR_TRACKS production.

usage: make_sample_lists.py <outdir> [step=10] [nper=10] [--complement] [--busy RUN1,...] [--nper-busy N] [--mid RUN1,...] [--nper-mid N]
  default     : segments with seg % step == 0           (step=10 -> the 10% sample)
  --complement: segments with seg % step != 0           (step=10 -> the remaining 90%)
  Segments are grouped per run into chunks of <= nper files (one Fun4All job per chunk; a chunk never mixes
  runs because run number / CDB timestamp come from the first file). Runs given with --busy (high pile-up,
  slow) use --nper-busy files per job instead.
writes <outdir>/lists/<run>_<chunk>.txt, <outdir>/joblists.txt (one list path per line), <outdir>/inventory.txt"""
import argparse, os, re
from collections import defaultdict

D = "/sphenix/lustre01/sphnxpro/production/run3pp/physics/ana581_2026p005_v001/DST_TRKR_TRACKS"
ap = argparse.ArgumentParser()
ap.add_argument("outdir"); ap.add_argument("step", nargs="?", type=int, default=10); ap.add_argument("nper", nargs="?", type=int, default=10)
ap.add_argument("--complement", action="store_true"); ap.add_argument("--busy", default=""); ap.add_argument("--nper-busy", type=int, default=4)
ap.add_argument("--mid", default=""); ap.add_argument("--nper-mid", type=int, default=10)
ap.add_argument("--all", action="store_true", help="select every segment (no seg %% step selection)")
ap.add_argument("--exclude", default="", help="file with DST paths already processed (one per line) to skip")
ap.add_argument("--min-age-h", type=float, default=0, help="skip DST files modified less than this many hours ago")
a = ap.parse_args()
out = os.path.abspath(a.outdir)
busy = {int(r) for r in a.busy.split(",") if r}
mid = {int(r) for r in a.mid.split(",") if r}
import time
excl = set(l.strip() for l in open(a.exclude)) if a.exclude else set()
tmin = time.time() - 3600 * a.min_age_h
nexcl = nyoung = 0
pat = re.compile(r"DST_TRKR_TRACKS_run3pp_ana581_2026p005_v001-(\d{8})-(\d{5})\.root$")
runs = defaultdict(list); ntot = 0
for rr in sorted(os.listdir(D)):
    d = os.path.join(D, rr)
    if not os.path.isdir(d): continue
    for f in os.listdir(d):
        m = pat.match(f)
        if not m: continue
        ntot += 1
        run, seg = int(m.group(1)), int(m.group(2))
        p = os.path.join(d, f)
        if p in excl: nexcl += 1; continue
        if not a.all and (seg % a.step != 0) != a.complement: continue
        runs[run].append((seg, p))
# age check (stat is slow on lustre: only the 50 highest segments of each run, where new files appear)
if a.min_age_h > 0:
    for run in runs:
        lst = sorted(runs[run]); keep = lst[:-50]
        for seg, p in lst[-50:]:
            if os.stat(p).st_mtime > tmin: nyoung += 1
            else: keep.append((seg, p))
        runs[run] = keep
os.makedirs(os.path.join(out, "lists"), exist_ok=True)
jl = []; nsel = 0; perrun = []
for run in sorted(runs):
    files = [p for s, p in sorted(runs[run])]
    nsel += len(files)
    nper = a.nper_busy if run in busy else (a.nper_mid if run in mid else a.nper)
    nj = 0
    for i in range(0, len(files), nper):
        lp = os.path.join(out, "lists", f"{run:08d}_{i//nper:04d}.txt")
        open(lp, "w").write("\n".join(files[i:i+nper]) + "\n"); jl.append(lp); nj += 1
    perrun.append(f"  {run} files {len(files)} files/job {nper} jobs {nj}")
open(os.path.join(out, "joblists.txt"), "w").write("\n".join(jl) + "\n")
sel = "all" if a.all else f"seg % {a.step} {'!=' if a.complement else '=='} 0"
sel += f"; excluded (already processed) {nexcl}; skipped younger than {a.min_age_h} h {nyoung}"
inv = (f"total DST files {ntot}\nruns {len(runs)}\nselected files ({sel}) {nsel}\njobs {len(jl)} "
       f"(<= {a.nper} files/job, mid runs <= {a.nper_mid}, busy runs <= {a.nper_busy})\n" + "\n".join(perrun) + "\n")
open(os.path.join(out, "inventory.txt"), "w").write(inv); print(inv)
