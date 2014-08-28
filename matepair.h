#pragma once
#include "fastqlib.h"
#include "math.h"
#include "stdint.h"

int hamming(string & s1,string & s2,int offset1,int offset2,int L,int maxd);
int overlap(string & s1,string & s2,int minoverlap,float similarity);
int partial_match(string & s1,string & a,int minoverlap,int maxdist);
int findAdapter(string & s,int minoverlap,float similarity,bool hamming);

class levenshtein {
 public:
  levenshtein(string s1);
  ~levenshtein() {delete column;}

  unsigned int *column;
  unsigned   int i, j, lastdiag, olddiag;
  unsigned int L1,L2;
  string s1;
  int distance(string & s2,int offset,int maxdist,int indel_penalty=2);
};

class matePair {
  public:
  matePair(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming);
  matePair();
  int build(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming);
  int clear();
  int trimUnknown();
  bool trimExternal(readPair & rp);
  int joinReads(fqread & r1,fqread & r2,fqread & output);
  int  set_preserve_mp(bool);
  bool joinreads,use_hamming,preserve_mp;
  fqread se;
  readPair mp,pe,unknown;
  int resolve_overhang(fqread & r1, fqread & r2,int a,int b);
  unsigned int ham_align(string & s1,string & s2);
  int minoverlap,  minlen;
  float similarity;
};

