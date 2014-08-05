// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Jared O'Connell
///

#pragma once
#include "utilityfunc.h"

class fqread {
 public:
  fqread();
  fqread(int L);
  fqread(string header,string dna,string line3,string qual);
  int l;
  bool filtered;//if filtered is true, it failed QC
  string h,s,l3,q;
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
};


class fastqWriter {
 public:
  fastqWriter(string fname);
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

