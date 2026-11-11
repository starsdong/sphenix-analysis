#!/bin/bash
# Clean, pinned environment for run3pp KFParticle analysis.
# Data produced with ana581 (ana581_2026p005_v001) -> use the same release.
# Do NOT use the login env (~/mysetup.sh: "new" + /sphenix/u/dongxs/install with old libs).
export TOP=/gpfs02/sphenix/user/dongxs/run3pp
unset MYINSTALL
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.581 > /dev/null
export MYINSTALL=$TOP/install
mkdir -p $MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL > /dev/null
