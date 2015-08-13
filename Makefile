CC = g++

ifndef BOOST_ROOT
    $(error BOOST_ROOT is undefined.  Point BOOST_ROOT at your boost installation ie. BOOST_ROOT/lib and BOOST_ROOT/include should exist)
endif


CFLAGS = -O3  -I$(BOOST_ROOT)/include
LFLAGS=  $(BOOST_ROOT)/lib/libboost_iostreams.a $(BOOST_ROOT)/lib/libboost_program_options.a -lz

all: mergeReads nxtrim



debug: CFLAGS = -Wall -g -I$(BOOST_ROOT)/include
debug: all

GIT_HASH := $(shell git describe --abbrev=4 --always )

VERSION = v0.3.2-alpha
GIT_VERSION =
ifneq "$(wildcard .git)" ""
GIT_VERSION = -$(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)$(GIT_VERSION)"' > $@


test: test.cpp fastqlib.o utilityfunc.o matepair.o
	$(CC) $(CFLAGS) test.cpp fastqlib.o utilityfunc.o matepair.o -o test   $(LFLAGS)
nxtrim: nxtrim.cpp fastqlib.o utilityfunc.o matepair.o fastqlib.o version.h
	$(CC) $(CFLAGS) nxtrim.cpp fastqlib.o utilityfunc.o matepair.o -o nxtrim  $(LFLAGS)
mergeReads: mergeReads.cpp fastqlib.o utilityfunc.o fastqlib.o githash.h version.h
	$(CC) $(CFLAGS)  mergeReads.cpp fastqlib.o utilityfunc.o -o mergeReads   $(LFLAGS)
matepair.o: matepair.cpp matepair.h fastqlib.h
	$(CC) $(CFLAGS) -c matepair.cpp
fastqlib.o: fastqlib.cpp fastqlib.h utilityfunc.h
	$(CC) $(CFLAGS) -c fastqlib.cpp
utilityfunc.o:  utilityfunc.cpp utilityfunc.h
	$(CC) $(CFLAGS) -c utilityfunc.cpp 
clean:
	rm *.o nxtrim test mergeReads version.h

