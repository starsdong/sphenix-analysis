#!/usr/bin/env python3
"""done_files.py <outdir> : print the KFP_D0_*.root names of finished OK jobs (condor_<list>.out starts with OK),
mapping each job list to its output name via the first DST of the list (KFP_D0_<run>_<5-digit segment>_00000.root)."""
import glob, os, re, sys

out = sys.argv[1]
jl = [l.strip() for l in open(os.path.join(out, 'joblists.txt')) if l.strip()]
bykey = {os.path.splitext(os.path.basename(p))[0]: p for p in jl}  # condor $Fn(lst) drops .txt
for o in glob.glob(os.path.join(out, 'logs', 'condor_*.out')):
    try:
        if not open(o).readline().startswith('OK'):
            continue
    except OSError:
        continue
    key = os.path.basename(o)[len('condor_'):-len('.out')]
    lst = bykey.get(key)
    if not lst:
        continue
    first = open(lst).readline().strip()
    m = re.search(r'-(\d{8})-(\d{5})\.root$', first)
    if not m:
        continue
    fn = 'KFP_D0_%s_%s_00000.root' % (m.group(1), m.group(2))
    if os.path.exists(os.path.join(out, fn)):
        print(fn)
