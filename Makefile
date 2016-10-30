Cxx = g++


CFLAGS = -O3 
LFLAGS=  -lz

all: mergeReads nxtrim

debug: CFLAGS = -Wall -g 
debug: all

GIT_HASH := $(shell git describe --abbrev=4 --always )

VERSION = v0.4.1
GIT_VERSION =
ifneq "$(wildcard .git)" ""
GIT_VERSION = -$(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)$(GIT_VERSION)"' > $@


unit_test: test.cpp fastqlib.o utilityfunc.o matepair.o
	$(CXX) $(CFLAGS) test.cpp fastqlib.o utilityfunc.o matepair.o -o unit_test   $(LFLAGS)
nxtrim: nxtrim.cpp fastqlib.o utilityfunc.o matepair.o fastqlib.o version.h
	$(CXX) $(CFLAGS) nxtrim.cpp fastqlib.o utilityfunc.o matepair.o -o nxtrim  $(LFLAGS)
mergeReads: mergeReads.cpp fastqlib.o utilityfunc.o fastqlib.o githash.h version.h
	$(CXX) $(CFLAGS)  mergeReads.cpp fastqlib.o utilityfunc.o -o mergeReads   $(LFLAGS)
matepair.o: matepair.cpp matepair.h fastqlib.h
	$(CXX) $(CFLAGS) -c matepair.cpp
fastqlib.o: fastqlib.cpp fastqlib.h utilityfunc.h
	$(CXX) $(CFLAGS) -c fastqlib.cpp
utilityfunc.o:  utilityfunc.cpp utilityfunc.h
	$(CXX) $(CFLAGS) -c utilityfunc.cpp 
test: nxtrim
	bash example/run_test.sh
clean:
	rm *.o nxtrim test mergeReads version.h


