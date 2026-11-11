#!/bin/bash
# build + install EventInfoTree into $TOP/install (delete build/ after editing Makefile.am)
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh
mkdir -p $TOP/build/EventInfoTree && cd $TOP/build/EventInfoTree
if [ ! -f Makefile ]; then $TOP/src/EventInfoTree/autogen.sh --prefix=$MYINSTALL || exit 1; fi
make -j4 install || exit 1
echo BUILD_OK
