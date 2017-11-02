CXX ?= g++
CC ?= gcc


CXXFLAGS ?= -O2
LFLAGS = -lz

all: nxtrim

debug: CXXFLAGS += -Wall -g
debug: all

GIT_HASH := $(shell git describe --abbrev=4 --always )

VERSION = v0.4.2
GIT_VERSION =
ifneq "$(wildcard .git)" ""
GIT_VERSION = -$(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)$(GIT_VERSION)"' > $@

OBJS=matepair.o fastqlib.o utilityfunc.o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

nxtrim: nxtrim.cpp $(OBJS) version.h
	$(CXX) $(CXXFLAGS) nxtrim.cpp $(OBJS)  $(LFLAGS) -o $@
matepair.o: matepair.cpp matepair.h fastqlib.h
fastqlib.o: fastqlib.cpp fastqlib.h utilityfunc.h
utilityfunc.o:  utilityfunc.cpp utilityfunc.h
test: nxtrim
	bash -e example/run_test.sh
ecmg: nxtrim
	cd test/;bash -e ecmg.sh
clean:
	rm $(OBJS) nxtrim test version.h
	rm -rf test/output_dir/
	rm test/*bam test/*pe.fastq.gz test/*mp.fastq.gz  test/*unknown.fastq.gz
