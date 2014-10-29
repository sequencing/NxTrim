#pragma once
#include "utilityfunc.h"

class fqread {
 public:
  fqread();
  fqread(int L);
  fqread(string header,string dna,string line3,string qual);
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
  int print();
  int clear();
};

class readPair {
 public:
  readPair();
  readPair(fqread read1,fqread read2);
  int rc();
  int set(fqread read1, fqread read2);
  fqread r1,r2;
  int l;
  bool filtered;
};

class fastqReader {
 public:
  fastqReader(string fname);
  fqread next();
  bool fin();

 private:
  ifile infile;
  bool warned;
};


class fastqWriter {
 public:
  fastqWriter();
  fastqWriter(string fname);
  int open(string fname);
  int write(fqread & read);
  int write(readPair & read);
 private:
  ofile outfile;
};

class pairReader {
 public:
  pairReader(string fname1,string fname2);
  readPair next();
  int print();
  bool getPair(readPair & p);

 private:
  fastqReader *f1,*f2;
};

