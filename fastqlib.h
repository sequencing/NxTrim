#pragma once
#include "utilityfunc.h"
#include <zlib.h>  
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)  
using namespace std;
class fqread
{
public:
    fqread();
    fqread(int L);
    fqread(string header,string dna,string line3,string qual);
    int set(string header,string dna,string line3,string qual);
    int l;
    bool filtered,description;//if filtered is true, it failed QC
    string h,s,l3,q;
    int notN();
    int notN(int a,int b);//tells us how many non-missing bases are in this window.
    fqread mask(int a,int b);//N masks the region a,b
    fqread mask();
    fqread window(int a,int b);
    fqread window(int a);
    fqread rc();
    void print();
    int clear();
};

class readPair
{
public:
    readPair();
    readPair(fqread read1,fqread read2);
    int rc();
    void print();
    fqread r1,r2;
    int l;
    bool filtered;
};

class fastqReader
{
public:
    fastqReader(string fname);
    ~fastqReader();
    int next(fqread & r);

private:
    bool warned;
    gzFile fp;
    kseq_t *seq;
};


class fastqWriter
{
public:
    fastqWriter();
    ~fastqWriter();
    fastqWriter(string fname);
    int open(string fname);
    int write(fqread & read);
    int write(readPair & read);
private:
    gzFile fp;
    bool _stdout;//true if writing to stdout
};

class pairReader
{
public:
    pairReader(string fname1,string fname2);
    int next(readPair & r);
    int print();
    bool getPair(readPair & p);
private:
    fastqReader *f1,*f2;
};


class pairWriter
{
public:
    pairWriter();

    //interleaved
    pairWriter(string fname);
    int open(string fname);
    //separate files
    pairWriter(string fname1,string fname2);
    int  open(string fname1,string fname2);
    int write(readPair & read);
    bool separate;
private:
    fastqWriter outfile;//interleaved
    fastqWriter outfile1,outfile2;//separate files
};
