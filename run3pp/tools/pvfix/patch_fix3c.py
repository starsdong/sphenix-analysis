#!/usr/bin/env python3
"""d0-pvfix fix 3c: track-covariance scale factor for the KFParticle PV fit from KFP_PVKF_ERRSCALE (default 1; the
covariance is multiplied by the square). Split-half tests on ana581 DSTs (tools/pvfix/pv_splithalf.C) give vertex
chi2/ndf ~ 52 and pulls 6-9 with the stored track errors, and chi2/ndf ~ 1, pulls 0.6-0.9 with a scale of 7.
Run in offline/packages/KFParticle_sPHENIX."""

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
c = sub(c, '''      std::cout << "d0-pvfix: KFParticle PV fit chi2-deviation cut " << v << std::endl;
      return v;
    }();''', '''      std::cout << "d0-pvfix: KFParticle PV fit chi2-deviation cut " << v << std::endl;
      return v;
    }();
    static const float errScale2 = [] {
      const char* e = std::getenv("KFP_PVKF_ERRSCALE");
      float v = e ? std::atof(e) : 1.F;
      std::cout << "d0-pvfix: KFParticle PV fit track covariance scaled by " << v << "^2" << std::endl;
      return v * v;
    }();''', nl)
c = sub(c, '''          cov[k++] = s->get_error(a, b);
        }
      }
      parts[i].Create(par, cov, (Int_t) s->get_charge(), -1);''', '''          cov[k++] = s->get_error(a, b) * errScale2;
        }
      }
      parts[i].Create(par, cov, (Int_t) s->get_charge(), -1);''', nl)
wr('KFParticle_nTuple.cc', c)
print('fix3c applied')
