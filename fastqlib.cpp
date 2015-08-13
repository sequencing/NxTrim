#include "fastqlib.h"
using namespace std;

fqread::fqread(string header,string dna,string line3,string qual){
  description=true;
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
  if(!iss) {
    filtered=false;
    description=false;
  }
  else{
    if(tmp.substr(1,3)==":Y:")
      filtered=true;
    else
      filtered=false;
  }
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

int fqread::notN() {
  return(notN(0,l));
}

int fqread::notN(int a,int b) {
  assert(b>a);
  int count = 0;
  for(int i=a;i<b;i++)
    if(s[i]=='N' || ((unsigned int)q[i])<((unsigned int)'%'))
      count++;
  return(b - a - count);
}

fqread fqread::mask() { //N masks the entire read
  return(mask(0,l));
}

fqread fqread::mask(int a,int b) { //N masks the region a,b
  assert(a>=0&&b<=l);
  string new_s = s;
  string new_q = q;
  for(int i=a;i<b;i++) {
    new_s[i]='N';
    new_q[i]='#';
  }
  return(fqread(h,new_s,l3,new_q));
}

fqread fqread::window(int a,int b) {
  return(fqread(h,s.substr(a,b-a),l3,q.substr(a,b-a)));
}

fqread fqread::window(int a) {
  return(fqread(h,s.substr(a),l3,q.substr(a)));
}

void fqread::print() {
  if(l>0) {
    cout << h << endl;
    cout << s << endl;
    cout <<l3 << endl;
    cout << q << endl;
  }
}


fastqReader::fastqReader(string fname){
  warned=false;
  if(! infile.open(fname) ) {
    cerr << "Problem reading "<<fname<<endl;
    exit(1);
  }
  if(!infile.good()) {
    cerr << "Problem reading "<<fname<<endl;
    exit(1);
  }
}

fastqWriter::fastqWriter(){
}

fastqWriter::fastqWriter(string fname){
  open(fname);
}

int fastqWriter::open(string fname){
  return(  outfile.open(fname) );
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
  fqread r(h,s,l3,q);
  if(!warned&&!r.description&&infile)  {
    cerr << "WARNING: no description found in read header.  Assuming read passed passed chastity/purity filters." << endl;
    warned=true;
  }    
  return(r);
}

pairReader::pairReader(string fname1,string fname2) {
  f1 = new fastqReader(fname1);
  f2 = new fastqReader(fname2);
}

readPair pairReader::next() {
  return( readPair(f1->next(),f2->next()));
}

void readPair::print() {
  r1.print();
  r2.print();
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

pairWriter::pairWriter(){
}

pairWriter::pairWriter(string fname) {
  open(fname);
}

int pairWriter::open(string fname) {
  outfile.open(fname);
  separate=false;
  return(0);
}

pairWriter::pairWriter(string fname1,string fname2) {
  open(fname1,fname2);
}

int pairWriter::open(string fname1,string fname2) {
  outfile1.open(fname1);
  outfile2.open(fname2);
  separate=true;
  return(0);
}

int pairWriter::write(readPair & p) {
  if(p.r1.l>0 && p.r2.l>0) {
    if(separate) {
      outfile1.write(p.r1);
      outfile2.write(p.r2);
    }
    else {
      outfile.write(p);
    }
    return(1);
  }
  else
    return(0);
}
