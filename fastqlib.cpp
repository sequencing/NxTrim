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


#include "fastqlib.h"
using namespace std;

fqread::fqread(string header,string dna,string line3,string qual){
  h=header;
  s=dna;
  l3=line3;
  q=qual;
  l=s.size();
  istringstream iss;
  iss.str(header);
  string tmp;
  iss >> tmp;
  iss >> tmp;
  if(tmp[2]=='N')
    filtered=false;
  else
    filtered=true;

  assert(dna.size()==qual.size());
}

fqread::fqread() {
  l=0;
}

int fqread::clear() {
  l=0;
  h="";
  s.clear();
  q.clear();
  l3="";  
  return(0);
}

fqread::fqread(int L) {
  l=L;
  h="";
  s.resize(L);
  q.resize(L);
  l3="";
}

readPair::readPair() {
}

readPair::readPair(fqread read1,fqread read2) {
  r1=read1;
  r2=read2;
  filtered = read1.filtered || read2.filtered;
}

int readPair::rc() {
  r1 = r1.rc();
  r2 = r2.rc();
  return(0);
}

fqread fqread::rc() {
  fqread ret(l);
  ret.h=h;
  ret.l3=l3;
  for(int i=0;i<l;i++) {
    ret.q[l-i-1]=q[i];
    char base='N';
    switch(s[i]) {
    case 'C': 
      base='G';
      break;
    case 'G': 
      base='C';
      break;
    case 'T': 
      base='A';
      break;
    case 'A': 
      base='T';         
      break;
    }
    ret.s[l-i-1]=base;
  }  
  return(ret);
}

fqread fqread::window(int a,int b) {
  return(fqread(h,s.substr(a,b-a),l3,q.substr(a,b-a)));
}

fqread fqread::window(int a) {
  return(fqread(h,s.substr(a),l3,q.substr(a)));
}

int fqread::print() {
  cout << h << endl;
  cout << s << endl;
  cout <<l3 << endl;
  cout << q << endl;
  return(0);
}


fastqReader::fastqReader(string fname){
  infile.open(fname);
  if(!infile) {
    cerr << "Problem reading "<<fname<<endl;
    exit(1);
  }
}

fastqWriter::fastqWriter(string fname){
  outfile.open(fname);
  if(!outfile) {
    cerr << "Problem writing to "<<fname<<endl;
    exit(1);
  }
}

int fastqWriter::write(fqread & read) {
  if(read.l>0) {
    outfile << read.h <<endl;
    outfile << read.s <<endl;
    outfile << read.l3 <<endl;
    outfile << read.q <<endl;
    return(1);
  }
  else  return(0);
}

int fastqWriter::write(readPair & p) {
  if(p.r1.l>0 && p.r2.l>0) {
    write(p.r1);
    write(p.r2);
    return(1);
  }
  else
    return(0);
}

bool fastqReader::fin() {
  if(infile) return(true);
  else return(false);
}

fqread fastqReader::next() {
  string s,h,l3,q;
  getline(infile,h);
  getline(infile,s);
  getline(infile,l3);
  getline(infile,q);

  return(fqread(h,s,l3,q));
}

pairReader::pairReader(string fname1,string fname2) {
  f1 = new fastqReader(fname1);
  f2 = new fastqReader(fname2);
}

readPair pairReader::next() {
  return( readPair(f1->next(),f2->next()));
}

int readPair::set(fqread read1, fqread read2) {
  r1 = read1;
  r2 = read2;
  filtered = read1.filtered || read2.filtered;  
  return(0);
}

bool pairReader::getPair(readPair & p) {
  if(f1->fin() && f2->fin()) {
    p.set(f1->next(),f2->next());
    return(f1->fin() && f2->fin());
  } else {
    return false;
  }
}
