AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libeventinfotree.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

pkginclude_HEADERS = EventInfoTree.h

libeventinfotree_la_SOURCES = EventInfoTree.cc

libeventinfotree_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lffaobjects \
  -lffarawobjects \
  -ltrackbase_historic_io \
  -lglobalvertex_io

BUILT_SOURCES = testexternals.cc
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libeventinfotree.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
