#!/usr/bin/env python3
"""d0-pvfix fix 1: primary_vertex_ID = SvtxVertexMap key of the chosen PV (was: its beam crossing, from getPVID());
the old value is kept in the new branch primary_vertex_crossing. Run in offline/packages/KFParticle_sPHENIX."""
import io

def rd(p):
    b = open(p, 'rb').read().decode()
    return b, ('\r\n' if '\r\n' in b else '\n')

def wr(p, s):
    open(p, 'wb').write(s.encode())

def sub(s, a, b, nl):
    a, b = a.replace('\n', nl), b.replace('\n', nl)
    assert s.count(a) == 1, ('not unique / missing', a[:80], s.count(a))
    return s.replace(a, b)

h, nl = rd('KFParticle_nTuple.h')
h = sub(h, '  int m_calculated_vertex_ID{-1};\n',
        '  int m_calculated_vertex_ID{-1};\n  int m_calculated_vertex_crossing{-1};  // d0-pvfix fix 1: beam crossing of the PV (was stored as primary_vertex_ID)\n', nl)
wr('KFParticle_nTuple.h', h)

c, nl = rd('KFParticle_nTuple.cc')
c = sub(c, '    m_tree->Branch("primary_vertex_ID", &m_calculated_vertex_ID, "primary_vertex_ID/I");\n',
        '    m_tree->Branch("primary_vertex_ID", &m_calculated_vertex_ID, "primary_vertex_ID/I");\n'
        '    m_tree->Branch("primary_vertex_crossing", &m_calculated_vertex_crossing, "primary_vertex_crossing/I");\n', nl)
c = sub(c, '''    if (m_use_fake_pv_nTuple)
    {
      m_calculated_vertex_ID = -100;  // error value returned by getPVID
    }
    else
    {
      m_calculated_vertex_ID = getPVID(topNode, vertex_fillbranch);
    }
''', '''    if (m_use_fake_pv_nTuple)
    {
      m_calculated_vertex_ID = -100;  // error value returned by getPVID
      m_calculated_vertex_crossing = -100;
    }
    else
    {
      // d0-pvfix fix 1: getPVID() returns the beam crossing of the vertex, not its identifier.
      // primary_vertex_ID is now the SvtxVertexMap key (set in KFParticle_Tools::makeAllPrimaryVertices);
      // the crossing is kept in primary_vertex_crossing.
      m_calculated_vertex_ID = (int) vertex_fillbranch.Id();
      m_calculated_vertex_crossing = getPVID(topNode, vertex_fillbranch);
    }
''', nl)
wr('KFParticle_nTuple.cc', c)
print('fix1 applied')
