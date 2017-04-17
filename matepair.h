#pragma once
#include "fastqlib.h"
#include "math.h"
#include "stdint.h"

extern "C" {
#include "ksw.h"
}

int hamming(string & s1,string & s2,int offset1,int offset2,int L,int maxd);
int overlap(string & s1,string & s2,int minoverlap,float similarity);



class matePair {
  public:
  matePair(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming,bool preserve_mp,bool jjmp);
  matePair();
  ~matePair();  
  int build(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming,bool preserve_mp,bool jjmp);
  int findAdapter(string & s,int minoverlap,float similarity,bool hamming);  
  int clear();
  int trimUnknown();
  bool trimExternal(readPair & rp);
  int joinReads(fqread & r1,fqread & r2,fqread & output);
  bool joinreads,use_hamming,preserve_mp,_justmp;
  fqread se;
  readPair mp,pe,unknown;
  int resolve_overhang(fqread & r1, fqread & r2,int a,int b);
  unsigned int ham_align(string & s1,string & s2);
  int minoverlap,  minlen;
  float similarity;

  //stuff for ksw
  uint8_t *adapter1_sw,*adapter2_sw,*adapterj_sw;
  int8_t sw_mat[25];  

//stuff for shredding
  vector<string> seeds;
  int nseed,seedsize;
};

//handles the output for nxtrim (which reads go to which file etc)
class nxtrimWriter { 

 public:
  nxtrimWriter();
  nxtrimWriter(string prefix,bool jmp,bool separate_read_files=false);
  int open(string prefix,bool jmp,bool separate_read_files);
  int open();
  bool setMP(bool val) {_write_mp=val;return(_write_mp);}
  bool setUN(bool val) {_write_un=val;return(_write_un);}

  int  write(matePair & m);
  int weird,n_mp,n_unk,n_se,n_pe;//counts for each virtual library
  bool _justmp,_write_mp,_write_se,_write_pe,_write_un;
  bool print_to_stdout;
  pairWriter mp_out;
  pairWriter pe_out;
  pairWriter unknown_out;
  fastqWriter se_out;
};
