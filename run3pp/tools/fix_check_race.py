#!/usr/bin/env python3
"""fix_check_race.py : in condor/check_prod.sh, query the queued jobs BEFORE listing the job logs, so a job that
finishes between the two steps is not counted as FAILED (plain log in the listing, already gone from the queue)."""
import sys

p = sys.argv[1]
s = open(p).read()
lines = s.split('\n')
iq = [i for i, l in enumerate(lines) if 'queued.txt' in l and 'condor_q' in l]
ils = [i for i, l in enumerate(lines) if l.startswith('ls -f $out/logs')]
assert len(iq) == 1 and len(ils) == 1, (iq, ils)
iq, ils = iq[0], ils[0]
if iq < ils:
    print('already ordered')
    sys.exit(0)
qline = lines[iq]
com = lines[iq - 1] if lines[iq - 1].startswith('#') else None
del lines[iq]
if com is not None:
    del lines[iq - 1]
ins = ['# jobs still in the queue (their logs are plain while running); queried BEFORE listing the logs so that a job',
       '# finishing in between is not flagged as failed', qline]
lines[ils:ils] = ins
open(p, 'w').write('\n'.join(lines))
print('reordered')
