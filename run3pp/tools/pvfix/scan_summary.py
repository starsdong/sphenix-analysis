#!/usr/bin/env python3
"""scan_summary.py <txt> [<txt> ...] [--m PDG] : for each selection of a search output (pvfix_d0_scan*.C or charm_search.C)
print the window z, and the spread of the window-scan z away from the signal (|centre - m| >= 40 MeV): r.m.s., and the
z calibrated by that r.m.s. (z / rms_scan). Scan points overlap (40 MeV window, 10 MeV step), so they are correlated."""
import re, sys, math
args = sys.argv[1:]
m0 = 1.86484
if '--m' in args:
    i = args.index('--m'); m0 = float(args[i + 1]); del args[i:i + 2]
for fn in args:
    s = open(fn).read()
    print(f'== {fn} (signal at {m0})')
    zwin = {}
    for line in s.splitlines():
        mm = re.match(r'^(\w+)\s.*excess\s+(-?\d+) \+-\s+(\d+) \((-?[\d.]+|-?inf) sigma\), narrow z\s+(-?[\d.]+)', line)
        if mm and not line.startswith('ALT'):
            zwin[mm.group(1)] = (int(mm.group(2)), int(mm.group(3)), float(mm.group(4)), float(mm.group(5)))
    for line in s.splitlines():
        mm = re.match(r'^SCAN (\w+) largest .*?;(?: centre:z/narrow_z)?(.*)$', line)
        if not mm:
            continue
        tag = mm.group(1)
        pts = []
        for tok in mm.group(2).split():
            c, z = tok.split(':')
            pts.append((float(c), float(z.split('/')[0])))
        away = [z for c, z in pts if abs(c - m0) >= 0.0399]
        rms = math.sqrt(sum(z * z for z in away) / len(away))
        nbig = sum(1 for z in away if abs(z) >= abs(zwin[tag][2]))
        ex, eex, z, nz = zwin[tag]
        near = ' '.join(f'{c:.3f}:{zz:+.1f}' for c, zz in pts if abs(c - m0) < 0.031)
        print(f'  {tag:7s} excess {ex:6d} +- {eex:4d}  z {z:+5.1f}  narrow {nz:+5.2f} | scan away: n {len(away):2d} rms {rms:4.2f}, |z|>=|z_win| {nbig:2d} | z/rms {z / rms:+5.2f} | near: {near}')
