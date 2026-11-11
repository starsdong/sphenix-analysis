#!/usr/bin/env python3
"""d0-pvfix fix 3b: chi2-deviation cut of the KFParticle PV fit (KFVertex::ConstructPrimaryVertex) from the environment
KFP_PVKF_CHICUT (default 3.5, the KFParticle default); the value is printed once. Run in offline/packages/KFParticle_sPHENIX."""

def rd(p):
    b = open(p, 'rb').read().decode()
    return b, ('\r\n' if '\r\n' in b else '\n')

def wr(p, s):
    open(p, 'wb').write(s.encode())

def sub(s, a, b, nl):
    a, b = a.replace('\n', nl), b.replace('\n', nl)
    assert s.count(a) == 1, ('not unique / missing', a[:80], s.count(a))
    return s.replace(a, b)

c, nl = rd('KFParticle_nTuple.cc')
c = sub(c, '''  int pvfixKFVertex(const std::vector<PVfixTrack>& t, KFParticle& out, int& nUsed)
  {''', '''  int pvfixKFVertex(const std::vector<PVfixTrack>& t, KFParticle& out, int& nUsed)
  {
    static const float chiCut = [] {
      const char* e = std::getenv("KFP_PVKF_CHICUT");
      float v = e ? std::atof(e) : 3.5F;
      std::cout << "d0-pvfix: KFParticle PV fit chi2-deviation cut " << v << std::endl;
      return v;
    }();''', nl)
c = sub(c, '''    v.ConstructPrimaryVertex(ptr.data(), (int) t.size(), flag.get(), 3.5);''', '''    v.ConstructPrimaryVertex(ptr.data(), (int) t.size(), flag.get(), chiCut);''', nl)
wr('KFParticle_nTuple.cc', c)
print('fix3b applied')
