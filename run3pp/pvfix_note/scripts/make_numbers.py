#!/usr/bin/env python3
"""make_numbers.py <note_root> : turn the logs of the pvfix analysis note into LaTeX macros (tables/numbers.tex) and
generated tables (tables/*.tex). Every number quoted in the note comes from here; PASS/FAIL checks are emitted as
macros. Macro names are letters only."""
import re
import statistics
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
assert (root / 'main.tex').exists() and (root / 'logs').is_dir(), root
L = root / 'logs'
T = root / 'tables'
T.mkdir(exist_ok=True)
M = {}


def mac(k, v):
    assert re.fullmatch(r'[A-Za-z]+', k), k
    assert k not in M, 'duplicate macro ' + k
    M[k] = str(v)


def rd(name):
    return (L / name).read_text()


def big(n):
    return f'{int(n):,}'.replace(',', '{,}')


def mil(n, d=1):
    return f'{n / 1e6:.{d}f}'


checks = []


def check(name, ok, text):
    mac('chk' + name, r'\textcolor{green!50!black}{PASS}' if ok else r'\textcolor{red}{FAIL}')
    checks.append((name, ok, text))


# ---------------------------------------------------------------- dataset / jobs
r10 = rd('runs10.txt')
runs = dict(x.split(':') for x in re.search(r'RUNS10 (\S+)', r10).group(1).strip(',').split(','))
mac('nRunsTen', len(runs))
mac('runFirst', min(runs)), mac('runLast', max(runs))
mac('runListTen', ', '.join(sorted(runs)))
mac('nListsTen', big(re.search(r'NLISTS10 (\d+)', r10).group(1)))
mac('nSegsTen', big(re.search(r'SEGS10 (\d+)', r10).group(1)))
for tag, f in (('Two', 'jobs_v2_10pct.txt'), ('One', 'jobs_v1_10pct.txt')):
    s = rd(f)
    m = re.search(r'total (\d+) \| done OK (\d+) \| FAILED (\d+)', s)
    mac('jobsTotal' + tag, big(m.group(1)))
    mac('jobsOk' + tag, big(m.group(2)))
    mac('jobsFail' + tag, m.group(3))
    ev = int(re.search(r'events processed: (\d+)', s).group(1))
    mac('events' + tag, mil(ev, 2))
    mac('size' + tag, re.search(r'size ([\d.]+) GB', s).group(1))
    if tag == 'Two':
        evTwo = ev
check('Jobs', M['jobsOk' + 'Two'] == M['jobsTotal' + 'Two'] and M['jobsFail' + 'Two'] == '0', 'all v2 10% jobs OK')
check('Events', M['eventsTwo'] == M['eventsOne'], 'same events processed in v1 and v2 (same job lists)')

sw = rd('software.txt')
mac('release', re.search(r'release: (\S+)', sw).group(1))
mac('rootVersion', re.search(r'^root (\S+)', sw, re.M).group(1))
commits = re.findall(r'commit (\w+) (.*)', sw)
with open(T / 'commits.tex', 'w') as f:
    for h, msg in reversed(commits):
        fix = re.search(r'\(d0-pvfix (fix [\w]+)\)', msg)
        desc = re.sub(r'\s*\(d0-pvfix fix \w+\)', '', msg).replace('KFParticle_nTuple: ', '').replace('KFParticle_eventReconstruction: ', '')
        desc = desc.replace('_', r'\_').replace('*', r'$^*$').replace('#', r'\#').replace('&', r'\&')
        f.write(rf'\texttt{{{h}}} & {fix.group(1) if fix else ""} & {desc} \\' + '\n')
mac('nCommits', len(commits))
mac('libMdfive', re.search(r'(\w{32})\s+\S*libkfparticle', sw).group(1)[:12])
mac('mapMdfive', re.search(r'(\w{32})\s+\S*kappa_lookup', sw).group(1)[:12])

# ---------------------------------------------------------------- smoke tests
sm = rd('smoke.txt')
blocks = re.findall(r'SMOKE run(\d+)\nA: (\d+) cand \| passOriginal (\d+) \(in stock (\d+)\) \| passRefit (\d+) \| passKF (\d+) \| only refit (\d+), only KF (\d+) \| stock cand missing (\d+)', sm)
for (run, n, po, ins, pr, pk, orf, ok_, miss), w in zip(blocks, ('A', 'B')):
    mac('smRun' + w, run)
    mac('smN' + w, n), mac('smPO' + w, po), mac('smOnlyR' + w, orf), mac('smOnlyK' + w, ok_), mac('smMiss' + w, miss)
cpu = {m[0]: (float(m[1]), int(m[2]), int(m[3]), int(m[4])) for m in re.findall(r'CPU (\S+) user ([\d.]+) size (\d+) entries (\d+) (\d+)', sm)}
mac('smStockA', cpu['stock'][2]), mac('smStockB', cpu['stock_81609'][2])
mbB = re.findall(r'^B: .*stock cand missing (\d+)', sm, re.M)
mac('smNoFiveMissA', mbB[0]), mac('smNoFiveMissB', mbB[1])
check('SmokeRepro', all(m[8] == '0' and m[2] == m[3] for m in blocks) and int(M['smPOA']) == cpu['stock'][2] and int(M['smPOB']) == cpu['stock_81609'][2],
      'smoke: every stock candidate present with passOriginal=1, passOriginal count == stock count')
mac('smCpuA', f"{100 * (cpu['v2fin'][0] / cpu['stock'][0] - 1):.0f}")
mac('smCpuB', f"{100 * (cpu['v2fin_81609'][0] / cpu['stock_81609'][0] - 1):.0f}")
mac('smBrStock', cpu['stock'][3]), mac('smBrFull', cpu['v2fin'][3]), mac('smBrSlim', cpu['v2slim'][3])
mac('smBytesCand', f"{(cpu['v2slim'][1] + cpu['v2slim_81609'][1]) / (cpu['v2slim'][2] + cpu['v2slim_81609'][2]):.0f}")
mac('smBytesStock', f"{(cpu['stock'][1] + cpu['stock_81609'][1]) / (cpu['stock'][2] + cpu['stock_81609'][2]):.0f}")
v = re.search(r'candidates: stock (\d+) closure (\d+) pvfix (\d+), in all three (\d+)', sm)
mac('valCommon', v.group(4))
dif = re.findall(r'differs: (\S+)\s+(\d+) of (\d+)', sm)
mac('valDiffBranches', len(dif))
mac('valDiffName', dif[0][0].replace('_', r'\_') if dif else 'none')
mac('valCrossEq', re.search(r'crossing == stock primary_vertex_ID for (\d+) of (\d+)', sm).group(1))
check('Branches', len(dif) == 1 and dif[0][0] == 'primary_vertex_ID', 'smoke: only primary_vertex_ID differs among the stock branches')

# ---------------------------------------------------------------- reproduction at scale, skim
rp = rd('repro.txt')
m = re.search(r'files compared (\d+) \(v1 file missing for (\d+)\) \| mismatching files (\d+)', rp)
mac('repFiles', big(m.group(1))), mac('repMissing', m.group(2)), mac('repMismatch', m.group(3))
m2 = re.search(r'v1 candidates (\d+) \| v2 candidates (\d+) \(x([\d.]+)\): passOriginal (\d+), passRefit (\d+), passKF (\d+), only refit (\d+), only KF (\d+)', rp)
n1, n2, ratio, po, pr, pk, onr, onk = m2.groups()
mac('repNone', mil(int(n1), 2)), mac('repNtwo', mil(int(n2), 2)), mac('repRatio', ratio)
mac('repPR', mil(int(pr), 2)), mac('repPK', mil(int(pk), 2)), mac('repOnlyR', mil(int(onr), 2)), mac('repOnlyK', mil(int(onk), 2))
mac('repNoneExact', big(n1))
check('Repro', m.group(3) == '0' and m.group(2) == '0' and n1 == po and m.group(1) == '5224', 'at scale: N(v2, passOriginal) == N(v1) in every file')
sk = re.search(r'skim_pvfix2: (\d+) -> (\d+) candidates \((\d+) branches\)', rd('skim.txt'))
mac('skimIn', mil(int(sk.group(1)), 2)), mac('skimOut', mil(int(sk.group(2)), 2)), mac('skimBr', sk.group(3))
check('SkimIn', sk.group(1) == n2, 'skim input == v2 candidates of the reproduction check')

# ---------------------------------------------------------------- track states, split half
ts = rd('trkstate.txt')
m = re.search(r'tracks (\d+), with a vertex id (\d+); vertex tracks (\d+), whose vertex_id == their vertex (\d+)', ts)
mac('tsTracks', big(m.group(1))), mac('tsVtxTracks', big(m.group(3)))
m = re.search(r'component along p: median \|\.\| ([\d.]+) um', ts)
mac('tsAlong', m.group(1))
m = re.search(r'\|state - vertex\| median ([\d.]+) um', ts)
mac('tsDist', m.group(1))


def split(fn):
    rows = {}
    s = rd(fn)
    mac_n = re.search(r'vertices (\d+), with >= (\d+) selected tracks (\d+)', s)
    for line in s.splitlines():
        mm = re.match(r'(pair|KF chi2cut \S+)\s+(\d+) \|\s+(\S+)\s+(\S+) \|\s+(\S+)\s+(\S+) \|\s+(\S+)\s+(\S+) \|\s+(\S+)\s+(\S+)', line)
        if mm:
            rows[mm.group(1)] = mm.groups()[1:]
    return mac_n, rows


nv, s1 = split('split_s1.txt')
_, s7 = split('split_s7.txt')
mac('splitVtx', big(nv.group(1))), mac('splitMinTrk', nv.group(2)), mac('splitSel', big(nv.group(3)))
names = {'pair': 'Pair', 'KF chi2cut 1e9': 'Kall', 'KF chi2cut 3.5': 'Kcut'}
for key, w in names.items():
    for tag, rows in (('One', s1), ('Seven', s7)):
        if key not in rows:
            continue
        r = rows[key]
        mac(f'sp{w}{tag}N', big(r[0])), mac(f'sp{w}{tag}X', r[1]), mac(f'sp{w}{tag}Z', r[2])
        if w != 'Pair':
            mac(f'sp{w}{tag}PullX', r[3]), mac(f'sp{w}{tag}PullZ', r[4]), mac(f'sp{w}{tag}Chi', r[7]), mac(f'sp{w}{tag}Use', r[8])
        else:
            if tag == 'One':
                mac('spPairCloX', r[5]), mac('spPairCloZ', r[6])
with open(T / 'splithalf.tex', 'w') as f:
    for key, lab in (('pair', 'pair average (PHSimpleVertexFinder-like)'), ('KF chi2cut 1e9', r'KF fit, no outlier rejection'),
                     ('KF chi2cut 9', r'KF fit, $\chi^2$ cut 9'), ('KF chi2cut 3.5', r'KF fit, $\chi^2$ cut 3.5')):
        for tag, rows, sc in (('1', s1, '1'), ('7', s7, '7')):
            if key not in rows or (key == 'pair' and tag == '7'):
                continue
            r = rows[key]
            pull = '--' if key == 'pair' else f'{r[3]} / {r[4]}'
            chi = '--' if key == 'pair' else r[7]
            use = '--' if key == 'pair' else r[8]
            scl = '--' if key == 'pair' else sc
            f.write(f'{lab} & {scl} & {r[0]} & {r[1]} / {r[2]} & {pull} & {chi} & {use} \\\\\n')

# ---------------------------------------------------------------- PV versions on the v2 ntuples
ps = rd('pv_shifts.txt')
mac('psN', mil(int(re.search(r'candidates (\d+)', ps).group(1)), 2))
mac('psFiles', sum(1 for _ in open(L / 'shift_files.list')))
with open(T / 'closure.tex', 'w') as f:
    for meth, lab in (('pair', 'pair refit'), ('kf', 'KF fit')):
        for mm in re.finditer(rf'CLOS {meth} (\S+) N (\d+) f1um ([\d.]+) q50 ([\d.]+) q90 ([\d.]+) dzw ([\d.]+)', ps):
            b, n, f1, q50, q90, dzw = mm.groups()
            f.write(f'{lab} & {b.replace(">=", "$\\geq$")} & {big(n)} & {float(f1) * 100:.1f} & {q50} & {q90} & {dzw} \\\\\n')
            if meth == 'pair' and b == '2':
                mac('cloPairTwoMed', q50), mac('cloPairTwoNinety', q90)
            if meth == 'pair' and b == '>=7':
                mac('cloPairManyNinety', q90), mac('cloPairManyOne', f'{float(f1) * 100:.0f}')
            if meth == 'kf' and b == '>=7':
                mac('cloKfManyMed', q50)
        f.write('\\midrule\n' if meth == 'pair' else '')
with open(T / 'shift.tex', 'w') as f:
    for meth, lab in (('pair', 'pair: refitAll $\\to$ refit'), ('kf', 'KF: kfAll $\\to$ kf')):
        for mm in re.finditer(rf'SHIFT {meth} (\d) N (\d+) q50 ([\d.]+) q90 ([\d.]+) dzw ([\d.]+)', ps):
            a, n, q50, q90, dzw = mm.groups()
            f.write(f'{lab} & {a} & {big(n)} & {q50} & {q90} & {dzw} \\\\\n')
            mac(f'sh{"Pair" if meth == "pair" else "Kf"}{["Zero", "One", "Two"][int(a)]}Med', q50)
            mac(f'sh{"Pair" if meth == "pair" else "Kf"}{["Zero", "One", "Two"][int(a)]}Dz', dzw)
check('ShiftZero', M['shPairZeroMed'] == '0.0' and M['shKfZeroMed'] == '0.0', 'no daughter attached -> refit with and without daughters identical')
st = re.search(r'STATUS N \d+ refitAll ([\d.]+) refit ([\d.]+) kfAll ([\d.]+) kf ([\d.]+) kfUse ([\d.]+)', ps)
mac('stRefitAll', f'{100 * float(st.group(1)):.2f}'), mac('stRefit', f'{100 * float(st.group(2)):.1f}')
mac('stKfAll', f'{100 * float(st.group(3)):.1f}'), mac('stKf', f'{100 * float(st.group(4)):.1f}'), mac('stKfUse', f'{100 * float(st.group(5)):.0f}')

# v1 effects (pair refit, 10% v1 ntuples)
ev = dict(l.split() for l in rd('effects10pct_values.txt').splitlines() if l.strip())
mac('efN', mil(int(ev['n']), 1))
for k, w in (('dcaMedO', 'DcaO'), ('dcaMedA', 'DcaA'), ('dcaMedR', 'DcaR'), ('lsMedO', 'LsO'), ('lsMedA', 'LsA'), ('lsMedR', 'LsR'),
             ('selTightAO', 'TightAO'), ('selTightRA', 'TightRA'), ('selPointRA', 'PointRA'), ('inBoth', 'InBoth'), ('inNone', 'InNone')):
    mac('ef' + w, ev[k])

# ---------------------------------------------------------------- pair crossing / pair cuts
for run, w in (('00079513', 'A'), ('00081609', 'B')):
    s = rd(f'pair_crossing_{run}.txt')
    mac('pcN' + w, big(re.search(r'candidates (\d+)', s).group(1)))
    m = re.search(r'same crossing: ([\d.]+) \| daughter crossing == PV crossing: ([\d.]+)', s)
    mac('pcSame' + w, f'{100 * float(m.group(1)):.1f}'), mac('pcPv' + w, f'{100 * float(m.group(2)):.1f}')
    check('Cross' + w, m.group(1) == '1.0000' and m.group(2) == '1.0000', f'run {run}: all pairs same crossing, PV same crossing')
    m = re.search(r'quantiles 5/25/50/75/95%: (\S+) \S+ \S+ \S+ (\S+)', s)
    mac('pcXlo' + w, m.group(1).replace('-', '$-$')), mac('pcXhi' + w, m.group(2))
    m = re.search(r'INTT hits: both daughters ([\d.]+)', s)
    mac('pcIntt' + w, f'{100 * float(m.group(1)):.1f}')
    m = re.search(r'\(>1: ([\d.]+)\)', s)
    mac('pcMultiPv' + w, f'{100 * float(m.group(1)):.0f}')
    m = re.search(r'nTracksOfBC\) quantiles 5/25/50/75/95%: (\d+) \d+ (\d+) \d+ (\d+)', s)
    mac('pcTrkMed' + w, m.group(2)), mac('pcTrkLo' + w, m.group(1)), mac('pcTrkHi' + w, m.group(3))
pcs = rd('pair_cuts.txt')
with open(T / 'paircuts.tex', 'w') as f:
    sect = None
    for line in pcs.splitlines():
        if line.startswith('== pointing'):
            sect = 'pointing'
        elif line.startswith('== control'):
            sect = 'control'
        mm = re.match(r'\s+nPrimaryVerticesOfBC(\S+)\s+N\s+(\d+) \| INTT ([\d.]+) \| same-coll ([\d.]+) \| both ([\d.]+) \| daughters attached: one ([\d.]+), none ([\d.]+) \| non-attached \|dz\|: median ([\d.]+), 75% ([\d.]+), 90% ([\d.]+) cm, > 1 mm ([\d.]+)', line)
        if mm and sect:
            pv, n, it, co, bo, one, none, dmed, d75, d90, gt = mm.groups()
            pvl = {'>=1': 'all', '==1': '1 PV', '>1': '$>$1 PV'}[pv]
            f.write(f'{sect} & {pvl} & {big(n)} & {100 * float(it):.1f} & {100 * float(co):.1f} & {100 * float(bo):.1f} & {10 * float(dmed):.2f} & {10 * float(d90):.1f} \\\\\n')
            key = ('Pt' if sect == 'pointing' else 'Ct') + {'>=1': 'All', '==1': 'One', '>1': 'Multi'}[pv]
            mac('cut' + key + 'Intt', f'{100 * float(it):.1f}'), mac('cut' + key + 'Coll', f'{100 * float(co):.1f}')
            mac('cut' + key + 'Gt', f'{100 * float(gt):.0f}'), mac('cut' + key + 'Dmed', f'{10 * float(dmed):.1f}'), mac('cut' + key + 'Dninety', f'{10 * float(d90):.1f}')
        if sect == 'control':
            pass

# ---------------------------------------------------------------- D0 search results
rx = re.compile(r'^([ORK])(Hi|All)([PT])\s+.*?Nsig\s+(\d+) Nctl\s+(\d+) \| raw ratio broad dchi2\s+(-?\d+) \| weighted: pol2 chi2/ndf\s+(\d+)/(\d+), broad dchi2\s+(-?\d+) \| '
                r'excess\s+(-?\d+) \+-\s+(\d+) \((-?[\d.]+) sigma\), narrow z\s+(-?[\d.]+) \| injected (\d+) -> recovered (-?\d+)%')
vname = {'O': 'Orig', 'R': 'Pair', 'K': 'Kf'}
res = {}
for var, w in (('d0cut_none', 'None'), ('d0cut_intt', 'Intt'), ('d0cut_coll', 'Coll'), ('d0cut_both', 'Both'), ('d0kc_none', 'KcNone'), ('d0kc_both', 'KcBoth')):
    s = rd(var + '.txt')
    m = re.search(r'with pair refit and KF fit (\d+) \(([\d.]+)%\), passing the extra selection (\d+) \(([\d.]+)% of these\)', s)
    mac('base' + w, mil(int(m.group(1)), 2)), mac('baseFrac' + w, m.group(2)), mac('extra' + w, mil(int(m.group(3)), 2)), mac('extraFrac' + w, f'{float(m.group(4)):.1f}')
    rows = {}
    for line in s.splitlines():
        mm = rx.match(line)
        if mm:
            g = mm.groups()
            rows[(g[0], g[1], g[2])] = dict(ns=int(g[3]), nc=int(g[4]), raw=int(g[5]), chi=int(g[6]), ndf=int(g[7]), bro=int(g[8]), ex=int(g[9]), eex=int(g[10]), z=float(g[11]),
                                            nz=float(g[12]), inj=int(g[13]), rec=int(g[14]))
    assert len(rows) == 12, (var, len(rows))
    res[w] = rows
    zs = [r['z'] for r in rows.values()]
    mac('zRms' + w, f'{statistics.pstdev(zs, 0.0) if False else (sum(z * z for z in zs) / len(zs)) ** 0.5:.2f}')
    mac('zMax' + w, f'{max(zs, key=abs):.1f}'.replace('-', '$-$'))
    mac('zAbsMax' + w, f'{max(abs(z) for z in zs):.1f}')
    mac('zGtTwo' + w, sum(1 for z in zs if abs(z) > 2))
    mac('zMean' + w, f'{sum(zs) / len(zs):+.2f}'.replace('+', '$+$').replace('-', '$-$'))
    mac('zPos' + w, sum(1 for z in zs if z > 0))
    nzs = [r['nz'] for r in rows.values()]
    mac('nzAbsMax' + w, f'{max(abs(z) for z in nzs):.1f}')
    for v_ in 'ORK':
        zv = [r['z'] for (a, b, c), r in rows.items() if a == v_]
        mac('zMax' + vname[v_] + w, f'{max(zv):.1f}'.replace('-', '$-$'))
    mac('recMin' + w, min(r['rec'] for r in rows.values())), mac('recMax' + w, max(r['rec'] for r in rows.values()))
    check('Inj' + w, all(abs(r['rec'] - 100) <= 10 for r in rows.values()), f'{var}: injected narrow peak recovered within 10% in all 12 selections')
    chis = [r['chi'] / r['ndf'] for r in rows.values()]
    mac('chiMax' + w, f'{max(chis):.1f}'), mac('chiMin' + w, f'{min(chis):.1f}')
    for (v_, rg, sel), r in rows.items():
        key = vname[v_] + rg + sel + w
        mac('ex' + key, str(r['ex']).replace('-', '$-$')), mac('eex' + key, r['eex']), mac('z' + key, f"{r['z']:.1f}".replace('-', '$-$'))
        mac('ns' + key, big(r['ns']))
selLab = {('Hi', 'P'): r'$\pt>1.5$, pointing', ('Hi', 'T'): r'$\pt>1.5$, tight', ('All', 'P'): r'all \pt, pointing', ('All', 'T'): r'all \pt, tight'}


def cell(r, bold=False):
    s = f"{r['ex']:+d} $\\pm$ {r['eex']} ({r['z']:+.1f})".replace('+', '$+$').replace('-', '$-$').replace('$$', '')
    return r'\textbf{' + s + '}' if (bold and abs(r['z']) >= 2) else s


def d0table(fn, keys):
    with open(T / fn, 'w') as f:
        for w, lab in keys:
            f.write(r'\multicolumn{4}{l}{\emph{' + lab + r'}} \\' + '\n')
            for rg, sel in (('Hi', 'P'), ('Hi', 'T'), ('All', 'P'), ('All', 'T')):
                cells = [cell(res[w][(v_, rg, sel)], True) for v_ in 'ORK']
                f.write(f'{selLab[(rg, sel)]} & ' + ' & '.join(cells) + r' \\' + '\n')


d0table('d0_uncorr.tex', (('None', 'no extra cut'), ('Intt', 'INTT on both daughters'), ('Coll', 'same collision'), ('Both', 'INTT + same collision')))
d0table('d0_kc.tex', (('KcNone', 'no extra cut'), ('KcBoth', 'INTT + same collision')))
with open(T / 'd0_counts.tex', 'w') as f:
    for rg, sel in (('Hi', 'P'), ('Hi', 'T'), ('All', 'P'), ('All', 'T')):
        cells = []
        for v_ in 'ORK':
            r = res['None'][(v_, rg, sel)]
            cells.append(f"{big(r['ns'])} / {big(r['nc'])}")
        f.write(f'{selLab[(rg, sel)]} & ' + ' & '.join(cells) + r' \\' + '\n')
with open(T / 'd0_model.tex', 'w') as f:
    for w, lab in (('None', 'uncorrected'), ('KcNone', r'$\kappa$-corrected')):
        for rg, sel in (('Hi', 'P'), ('Hi', 'T'), ('All', 'P'), ('All', 'T')):
            cells = []
            for v_ in 'ORK':
                r = res[w][(v_, rg, sel)]
                cells.append(f"{r['chi']}/{r['ndf']} & {r['raw']} & {r['bro']}")
            f.write(f'{lab}, {selLab[(rg, sel)]} & ' + ' & '.join(cells) + r' \\' + '\n')
# corrected vs uncorrected: same candidates? counts change only by D0_pT_kc (pT>1.5) and asym; report ratio of N
mac('nsChangeHiP', f"{100 * (res['KcNone'][('O', 'Hi', 'P')]['ns'] / res['None'][('O', 'Hi', 'P')]['ns'] - 1):+.1f}".replace('+', '$+$'))
check('AllPtSame', all(res['KcNone'][(v_, 'All', s_)]['ns'] == res['None'][(v_, 'All', s_)]['ns'] for v_ in 'ORK' for s_ in 'PT'),
      'all-pT selections contain the same candidates with and without the kappa correction (only D0 pT changes)')
# version consistency: number of version pairs differing by more than 2 sigma (independent-error approx.) is not computed; quote max z only

# ---------------------------------------------------------------- kappa
kc = dict(re.findall(r'^KC (\w+) (\S+)', rd('kc_stats.txt'), re.M))
mac('kcN', mil(int(kc['N']), 2))
check('KcFinite', kc['N'] == kc['Nfinite'], 'all kc-skim candidates have finite recomputed and corrected masses')
check('KcSameN', kc['N'] == sk.group(2), 'kc skim has the same candidates as the skim')
for k in ('kappaPlusMean', 'kappaPlusRms', 'kappaMinusMean', 'kappaMinusRms', 'dmRcMed', 'dmRcHw', 'dmKcMed', 'dmKcHw'):
    mac(k, kc[k].replace('-', '$-$'))
mac('dmRcTail', f"{100 * float(kc['dmRcTail']):.2f}"), mac('dmKcTail', f"{100 * float(kc['dmKcTail']):.0f}")
mac('ptRelHw', f"{100 * float(kc['ptRelHw']):.0f}"), mac('ptRelMed', f"{100 * float(kc['ptRelMed']):.0f}"), mac('ptZero', big(kc['NptZero']))
mac('ptRelTail', f"{100 * float(kc['ptRelTail']):.0f}")
vc = rd('kappa_v0check.txt')
with open(T / 'v0check.tex', 'w') as f:
    for ch, w in (('K0S', 'Ks'), ('Lambda', 'La')):
        blk = vc.split('== ' + ch)[1].split('==')[0]
        mac('vN' + w, big(re.search(r'(\d+) candidates', blk).group(1)))
        f.write(r'\multicolumn{6}{l}{\emph{' + (r'\Ks$\to\pi^+\pi^-$' if w == 'Ks' else r'$\Lambda\to p\pi^-$ (+ c.c.)') + r'}} \\' + '\n')
        for line in blk.splitlines():
            mm = re.match(r'\s+(.+?)\s+\| mu\s+([\d.]+)\+-\s*([\d.]+) sig\s+([\d.]+)\+-\s*([\d.]+) \| mu\s+([\d.]+)\+-\s*([\d.]+) sig\s+([\d.]+)\+-\s*([\d.]+) \| mu\s+([\d.]+)\+-\s*([\d.]+) sig\s+([\d.]+)\+-\s*([\d.]+) \| ([\d.]+)', line)
            if not mm:
                continue
            g = mm.groups()
            lab = g[0]
            sig_rc, sig_kc, esr, esk = float(g[7]), float(g[11]), float(g[8]), float(g[12])
            # robust flag: fitted widths away from the limits (0.2 s0, 5 s0) and errors < 30% of the width
            s0 = 6.0 if w == 'Ks' else 2.5
            ok = all(0.25 * s0 < s < 4.9 * s0 and e < 0.3 * s for s, e in ((sig_rc, esr), (sig_kc, esk)))
            ratio = f'{sig_kc / sig_rc:.2f}' if ok else '--'
            labt = lab.replace('pT ', r'$\pt$ ').replace('training runs 79509-79516', 'training runs 79509--79516')
            f.write(f'{labt} & {g[5]} & {g[7]} $\\pm$ {g[8]} & {g[9]} & {g[11]} $\\pm$ {g[12]} & {ratio} \\\\\n')
            key = {'all': 'All', 'training runs 79509-79516': 'Train', 'other runs': 'Other'}.get(lab)
            if key:
                mac(f'v{w}{key}MuRc', g[5]), mac(f'v{w}{key}MuKc', g[9]), mac(f'v{w}{key}SigRc', g[7]), mac(f'v{w}{key}SigKc', g[11]), mac(f'v{w}{key}Ratio', f'{sig_kc / sig_rc:.2f}')
                mac(f'v{w}{key}MuSt', g[1]), mac(f'v{w}{key}SigSt', g[3])

km = rd('kappa_map.txt')
for n, w in (('pT', 'Pt'), ('eta', 'Eta'), ('phi', 'Phi')):
    mm = re.search(rf'MAP {n} n (\d+) min (\S+) max (\S+)', km)
    mac('map' + w + 'N', mm.group(1)), mac('map' + w + 'Min', mm.group(2)), mac('map' + w + 'Max', mm.group(3))
mac('mapRows', big(re.search(r'MAP rows (\d+)', km).group(1)))

# ---------------------------------------------------------------- full sample (100%)
import math


def parse_search(fn, m0):
    """rows {tag: dict(ex, eex, z, nz, chi, ndf, ns)} and scan {tag: [(centre, z)]} of a search output"""
    s = rd(fn)
    rows, scans = {}, {}
    for line in s.splitlines():
        mm = re.match(r'^(\w+)\s.*Nsig\s+(\d+) Nctl\s+(\d+) \|.*?pol2 chi2/ndf\s+(\d+)/(\d+), broad dchi2\s+(-?\d+) \| excess\s+(-?\d+) \+-\s+(\d+) \((-?[\d.]+) sigma\), narrow z\s+(-?[\d.]+) \| injected (\d+) -> recovered (-?\d+)', line)
        if mm:
            g = mm.groups()
            rows[g[0]] = dict(ns=int(g[1]), nc=int(g[2]), chi=int(g[3]), ndf=int(g[4]), bro=int(g[5]), ex=int(g[6]), eex=int(g[7]), z=float(g[8]), nz=float(g[9]), rec=int(g[11]))
        mm = re.match(r'^SCAN (\w+) largest .*?;(?: centre:z/narrow_z)?(.*)$', line)
        if mm:
            scans[mm.group(1)] = [(float(t.split(':')[0]), float(t.split(':')[1].split('/')[0])) for t in mm.group(2).split()]
    for tag, r in rows.items():
        if tag in scans:
            away = [z for c, z in scans[tag] if abs(c - m0) >= 0.0399]
            r['rms'] = math.sqrt(sum(z * z for z in away) / len(away))
            r['nbig'] = sum(1 for z in away if abs(z) >= abs(r['z']))
            r['nway'] = len(away)
            r['zc'] = r['z'] / r['rms']
    return rows, scans


s9 = rd('jobs_v2_90pct.txt')
m9 = re.search(r'total (\d+) \| done OK (\d+) \| FAILED (\d+)', s9)
mac('jobsTotalNin', big(m9.group(1))), mac('jobsOkNin', big(m9.group(2))), mac('jobsFailNin', m9.group(3))
ev9 = int(re.search(r'events processed: (\d+)', s9).group(1))
mac('eventsNin', mil(ev9, 1)), mac('eventsAll', mil(ev9 + evTwo, 1)), mac('sizeNin', re.search(r'size ([\d.]+) GB', s9).group(1))
check('JobsNinety', m9.group(1) == m9.group(2) and m9.group(3) == '0', 'all v2 90% jobs OK')
rp = rd('repro90.txt')
m = re.search(r'files compared (\d+) \(v1 file missing for (\d+)\) \| mismatching files (\d+)', rp)
m2 = re.search(r'v1 candidates (\d+) \| v2 candidates (\d+) \(x([\d.]+)\).*only refit (\d+), only KF (\d+)', rp)
mac('ninFiles', big(m.group(1))), mac('ninMismatch', m.group(3)), mac('ninNone', big(m2.group(1))), mac('ninNtwo', mil(int(m2.group(2)), 1)), mac('ninRatio', m2.group(3))
mac('ninOnlyR', mil(int(m2.group(4)), 1)), mac('ninOnlyK', mil(int(m2.group(5)), 1))
check('ReproNinety', m.group(3) == '0' and m.group(2) == '0', '90%: N(v2, passOriginal) == N(stock) in every one of the files')
ska = re.search(r'skim_pvfix2: (\d+) -> (\d+)', rd('skim_all.txt'))
mac('skimAllIn', mil(int(ska.group(1)), 1)), mac('skimAllOut', mil(int(ska.group(2)), 1))
mac('candAllStock', mil(int(m2.group(1)) + int(re.search(r'v1 candidates (\d+)', rd('repro.txt')).group(1)), 1))
full = {}
for fn, w in (('fu_unc_both', 'FuUnc'), ('fu_kc_both', 'FuKc'), ('fu_kc_none', 'FuKcNone'), ('fu_kc_both_r795', 'RunA'), ('fu_kc_both_r809', 'RunB'), ('fu_kc_both_r816', 'RunC')):
    rows, scans = parse_search(fn + '.txt', 1.86484)
    assert len(rows) == 12, (fn, len(rows))
    full[w] = rows
    mm = re.search(r'passing the extra selection (\d+)', rd(fn + '.txt'))
    mac('nSel' + w, mil(int(mm.group(1)), 1))
    for tag, r in rows.items():
        k = vname[tag[0]] + tag[1:] + w
        mac('ex' + k, str(r['ex']).replace('-', '$-$')), mac('eex' + k, r['eex']), mac('z' + k, f"{r['z']:.1f}".replace('-', '$-$'))
        if 'zc' in r:
            mac('zc' + k, f"{r['zc']:.1f}".replace('-', '$-$')), mac('rms' + k, f"{r['rms']:.1f}"), mac('nbig' + k, r['nbig'])
        mac('chi' + k, f"{r['chi']}/{r['ndf']}")
    check('Inj' + w, all(abs(r['rec'] - 100) <= 10 for r in rows.values()), f'{fn}: injected peak recovered within 10% in all 12 selections')
check('RunSum', abs(sum(float(M['nSel' + w]) for w in ('RunA', 'RunB', 'RunC')) - float(M['nSelFuKc'])) < 0.15, 'the three run groups add up to the full sample (candidates)')
mac('nwayScan', full['FuKc']['OHiT']['nway'])
selLab2 = {'HiP': r'$\pt>1.5$, pointing', 'HiT': r'$\pt>1.5$, tight', 'AllP': r'all \pt, pointing', 'AllT': r'all \pt, tight'}


def fcell(r):
    s = f"{r['ex']:+d} $\\pm$ {r['eex']} ({r['z']:+.1f}; {r['zc']:+.1f})".replace('+', '$+$').replace('-', '$-$').replace('$$', '')
    return r'\textbf{' + s + '}' if abs(r['zc']) >= 2 else s


with open(T / 'full_d0.tex', 'w') as f:
    for w, lab in (('FuUnc', 'stored momenta, INTT + same collision'), ('FuKc', r'$\kappa$-corrected, INTT + same collision'), ('FuKcNone', r'$\kappa$-corrected, no extra cut')):
        f.write(r'\multicolumn{4}{l}{\emph{' + lab + r'}} \\' + '\n')
        for sel in ('HiP', 'HiT', 'AllP', 'AllT'):
            f.write(selLab2[sel] + ' & ' + ' & '.join(fcell(full[w][v + sel]) for v in 'ORK') + r' \\' + '\n')
with open(T / 'full_runs.tex', 'w') as f:
    for sel in ('HiT', 'AllT'):
        for v in 'ORK':
            cells = []
            for w in ('RunA', 'RunB', 'RunC', 'FuKc'):
                r = full[w][v + sel]
                cells.append(f"{r['ex']:+d} $\\pm$ {r['eex']}".replace('+', '$+$').replace('-', '$-$').replace('$$', ''))
            ex = [full[w][v + sel]['ex'] for w in ('RunA', 'RunB', 'RunC')]
            er = [full[w][v + sel]['eex'] for w in ('RunA', 'RunB', 'RunC')]
            wm = sum(e / s ** 2 for e, s in zip(ex, er)) / sum(1 / s ** 2 for s in er)
            chi = sum((e - wm) ** 2 / s ** 2 for e, s in zip(ex, er))
            # consistency with a common yield per event-fraction is not modelled; chi2 of equal yields (groups of similar size)
            f.write(f"{selLab2[sel]}, {['original PV', 'pair refit', 'KF fit']['ORK'.index(v)]} & " + ' & '.join(cells) + f' & {chi:.1f}' + r' \\' + '\n')
            mac('runChi' + vname[v] + sel, f'{chi:.1f}')
with open(T / 'full_scan.tex', 'w') as f:
    for v in 'ORK':
        for sel in ('HiT', 'AllT'):
            sc = parse_search('fu_kc_both.txt', 1.86484)[1][v + sel]
            near = ' & '.join(f'{zz:+.1f}' for c, zz in sc if 1.8249 < c < 1.9151).replace('+', '$+$').replace('-', '$-$').replace('$$', '')
            f.write(f"{['original PV', 'pair refit', 'KF fit']['ORK'.index(v)]}, {selLab2[sel]} & {near} \\\\\n")
scanCentres = [c for c, z in parse_search('fu_kc_both.txt', 1.86484)[1]['OHiT'] if 1.8249 < c < 1.9151]
mac('scanCentres', ' & '.join(f'{c:.3f}' for c in scanCentres))
mac('nScanCentres', len(scanCentres))

# ---------------------------------------------------------------- 3-prong charm (stock PV, 10%)
inv = rd('charm_inventory.txt')
mac('chSegs', big(re.search(r'^selected files .*\) (\d+)$', inv, re.M).group(1)))
mac('chRuns', re.search(r'^runs (\d+)', inv, re.M).group(1))
mac('chJobs', big(re.search(r'^jobs (\d+)', inv, re.M).group(1)))
ch = {}
for c, m0 in (('Dplus', 1.86966), ('Ds', 1.96835), ('Lc', 2.28646)):
    for cut in ('none', 'intt'):
        rows, scans = parse_search(f'c_{c}_{cut}.txt', m0)
        assert len(rows) == 8, (c, cut, len(rows))
        ch[(c, cut)] = rows
        s = rd(f'c_{c}_{cut}.txt')
        mm = re.search(r'candidates (\d+), after the offline base selection (\d+)', s)
        cc = {'Dplus': 'Dp', 'Ds': 'Ds', 'Lc': 'Lc'}[c] + {'none': 'None', 'intt': 'Intt'}[cut]
        mac('chN' + cc, mil(int(mm.group(1)), 2)), mac('chBase' + cc, mil(int(mm.group(2)), 2))
        for tag, r in rows.items():
            mac('chZ' + cc + tag.replace('st', 'St').replace('kc', 'Kc'), f"{r['z']:.1f}".replace('-', '$-$'))
            mac('chZc' + cc + tag.replace('st', 'St').replace('kc', 'Kc'), f"{r['zc']:.1f}".replace('-', '$-$'))
            mac('chEx' + cc + tag.replace('st', 'St').replace('kc', 'Kc'), f"{r['ex']} $\\pm$ {r['eex']}".replace('-', '$-$'))
        check('ChInj' + cc, all(abs(r['rec'] - 100) <= 10 for r in rows.values()), f'{c} {cut}: injected peak recovered within 10%')
        if c == 'Ds':
            alt = re.findall(r'^ALT (\w+) window .*?\((-?[\d.]+) sigma\)', s, re.M)
            mac('chAltMax' + cut.capitalize(), f"{max(abs(float(z)) for t, z in alt):.1f}")
with open(T / 'charm.tex', 'w') as f:
    for c, lab in (('Dplus', r'D$^+\to K^-\pi^+\pi^+$'), ('Ds', r'D$_s^+\to\phi\pi^+$'), ('Lc', r'$\Lambda_c^+\to pK^-\pi^+$')):
        f.write(r'\multicolumn{5}{l}{\emph{' + lab + r'}} \\' + '\n')
        for sel in ('HiP', 'HiT', 'AllP', 'AllT'):
            cells = []
            for cut in ('none', 'intt'):
                for v in ('st', 'kc'):
                    r = ch[(c, cut)][v + sel]
                    s_ = f"{r['ex']:+d} $\\pm$ {r['eex']} ({r['z']:+.1f}; {r['zc']:+.1f})".replace('+', '$+$').replace('-', '$-$').replace('$$', '')
                    cells.append(r'\textbf{' + s_ + '}' if abs(r['zc']) >= 2 else s_)
            f.write(selLab2[sel] + ' & ' + ' & '.join(cells) + r' \\' + '\n')

# ---------------------------------------------------------------- write
with open(T / 'numbers.tex', 'w') as f:
    for k, v in M.items():
        f.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
with open(T / 'checks.tex', 'w') as f:
    for name, ok, text in checks:
        f.write(f'\\item \\chk{name}{{}} -- {text.replace("_", chr(92) + "_").replace("%", chr(92) + "%")}\n'.replace('\\chk' + name, '\\chk' + name))
with open(L / 'numbers_summary.txt', 'w') as f:
    for k, v in M.items():
        f.write(f'{k} {v}\n')
nfail = sum(1 for c in checks if not c[1])
print(f'make_numbers: {len(M)} macros, {len(checks)} checks, {nfail} FAIL')
for c in checks:
    if not c[1]:
        print('  FAIL', c[0], c[2])
