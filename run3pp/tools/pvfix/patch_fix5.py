#!/usr/bin/env python3
"""d0-pvfix fix 5: best-combination choice compatible with the stock selection when the in-job DIRA cut is loosened.
KFParticle keeps one (mass assignment x PV) combination per track set, the one with the smallest PV-DCA chi2
(selectBestCombination). With a looser in-job DIRA cut, a combination that failed the stock DIRA cut can win and the
stock combination is lost. With KFP_PVFIX_PREFER_DIRA=<value> in the environment the choice is restricted to
combinations with DIRA (w.r.t. their PV) >= value whenever at least one exists; otherwise all combinations compete as
before. With value = the stock cut, every stock candidate is reproduced exactly. Run in offline/packages/KFParticle_sPHENIX."""

def rd(p):
    b = open(p, 'rb').read().decode()
    return b, ('\r\n' if '\r\n' in b else '\n')

def wr(p, s):
    open(p, 'wb').write(s.encode())

def sub(s, a, b, nl):
    a, b = a.replace('\n', nl), b.replace('\n', nl)
    assert s.count(a) == 1, ('not unique / missing', a[:80], s.count(a))
    return s.replace(a, b)

c, nl = rd('KFParticle_eventReconstruction.cc')
c = sub(c, '''  KFParticle smallestMassError = possibleCandidates[0];
  int bestCombinationIndex = 0;
  for (unsigned int i = 1; i < possibleCandidates.size(); ++i)
  {
''', '''  // d0-pvfix fix 5: prefer combinations passing KFP_PVFIX_PREFER_DIRA (w.r.t. their PV) when any exists
  static const float preferDira = [] {
    const char* e = std::getenv("KFP_PVFIX_PREFER_DIRA");
    float v = e ? std::atof(e) : -2.F;
    if (e)
    {
      std::cout << "d0-pvfix: best combination chosen among those with DIRA >= " << v << " when any exists" << std::endl;
    }
    return v;
  }();
  std::vector<bool> allowed(possibleCandidates.size(), true);
  if (preferDira > -1.5F && PVconstraint && !isAnInterMother)
  {
    bool any = false;
    for (unsigned int i = 0; i < possibleCandidates.size(); ++i)
    {
      allowed[i] = eventDIRA(possibleCandidates[i], possibleVertex[i]) >= preferDira;
      any = any || allowed[i];
    }
    if (!any)
    {
      allowed.assign(possibleCandidates.size(), true);
    }
  }
  unsigned int first = 0;
  while (!allowed[first])
  {
    ++first;
  }
  KFParticle smallestMassError = possibleCandidates[first];
  int bestCombinationIndex = (int) first;
  for (unsigned int i = first + 1; i < possibleCandidates.size(); ++i)
  {
    if (!allowed[i])
    {
      continue;
    }
''', nl)
if '#include <cstdlib>' not in c:
    c = sub(c, '#include "KFParticle_eventReconstruction.h"\n', '#include "KFParticle_eventReconstruction.h"\n\n#include <cstdlib>  // d0-pvfix fix 5\n#include <iostream>\n', nl)
wr('KFParticle_eventReconstruction.cc', c)
print('fix5 applied')
