AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libtrackqatree.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(OFFLINE_MAIN)/include/eigen3 \
  -isystem`root-config --incdir`

pkginclude_HEADERS = TrackQATree.h

libtrackqatree_la_SOURCES = TrackQATree.cc

libtrackqatree_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lffaobjects \
  -lffarawobjects \
  -ltrackbase_historic_io \
  -ltrackbase_historic \
  -lglobalvertex_io

BUILT_SOURCES = testexternals.cc
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libtrackqatree.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
