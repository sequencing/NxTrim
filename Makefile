CC = g++

ifndef BOOST_ROOT
    $(error BOOST_ROOT is undefined)
endif

LFLAGS = -L$(BOOST_ROOT)/lib -lz -lboost_iostreams  -lboost_program_options
CFLAGS = -O3  -I$(BOOST_ROOT)/include
#CFLAGS = -I -g  -I$(BOOST_ROOT)/include


all: mergeReads nxtrim test
test: test.cpp fastqlib.o utilityfunc.o matepair.o
	$(CC) $(CFLAGS) test.cpp fastqlib.o utilityfunc.o matepair.o -o test   $(LFLAGS) 
nxtrim: nxtrim.cpp fastqlib.o utilityfunc.o matepair.o fastqlib.o 
	$(CC) $(CFLAGS) nxtrim.cpp fastqlib.o utilityfunc.o matepair.o -o nxtrim  $(LFLAGS)
mergeReads: mergeReads.cpp fastqlib.o utilityfunc.o fastqlib.o 
	$(CC) $(CFLAGS)  mergeReads.cpp fastqlib.o utilityfunc.o -o mergeReads   $(LFLAGS)
matepair.o: matepair.cpp matepair.h fastqlib.h
	$(CC) $(CFLAGS) -c matepair.cpp
fastqlib.o: fastqlib.cpp fastqlib.h utilityfunc.h
	$(CC) $(CFLAGS) -c fastqlib.cpp
utilityfunc.o:  utilityfunc.cpp utilityfunc.h
	$(CC) $(CFLAGS) -c utilityfunc.cpp 
clean:
	rm *.o
	rm nxtrim
	rm test
	rm mergeReads
