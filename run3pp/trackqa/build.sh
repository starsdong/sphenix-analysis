#!/bin/bash
# Build + install TrackQATree into $TOP/install (shared with eventinfotree).
# After editing Makefile.am delete $TOP/build/TrackQATree (an existing Makefile is not regenerated).
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh
mkdir -p $TOP/build/TrackQATree && cd $TOP/build/TrackQATree
if [ ! -f Makefile ]; then bash $TOP/trackqa/src/TrackQATree/autogen.sh --prefix=$MYINSTALL || exit 1; fi
make -j4 install || exit 1
echo BUILD_OK
