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
#include "matepair.h"
             
//string r1_external_adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
string r1_external_adapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC";
//string r2_external_adapter = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
string r2_external_adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATC";
string adapter1 = "CTGTCTCTTATACACATCT";
string adapter2 = "AGATGTGTATAAGAGACAG";
string adapterj = adapter1+adapter2;

#define DEBUG 0

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))


int hamming(string & s1,string & s2,int offset1, int offset2,int L,int maxd) {

  int d = 0;
  //  cout << s1 << endl << s2 << " " << offset1 << " " << offset2 << " " << L << endl;
  for(int i=0;i<L;i++) {
    int j1=offset1+i;
    int j2=offset2+i;

    if(j1>=0 && j1<(int)s1.size() && j2>=0 && j2<(int)s2.size() ) {
      if(s1[j1]!=s2[j2])
	d++;
      if(d>maxd) {
	d=s1.length();
	break;
      }
    }
    //    cout << j1 << " " <<j2<<" "<<s1[j1] << " " << s2[j2] << " " << d <<endl;
  }
  //  cout << "d = "<<d << endl;
  return(d);
}


//only handles substitution errors in adapter (faster)
int partial_match(string & s1,string & s2,int minoverlap,float similarity) {
  assert(s2.size()<s1.size());
  int mini=s1.size(),mind=s2.size();
  assert((int)s1.size()>=minoverlap);
  
  int start = -(s2.size()-minoverlap);
  int stop = s1.size() - minoverlap;
  //  cout << "Range: " << start << " " << stop << endl;
  for(int i=start;i<stop;i++) {
    float compare_length = s2.size();
    if( i<0) compare_length += i;
    if( i>=(int)s1.size()) compare_length -= (i-(float)s1.size());
    int maxdist = ceil ( (1.-similarity) * compare_length);
    int d = hamming(s1,s2,i,0,s2.size(),maxdist);
    //    cout << i<<" "<<maxdist<<" "<<d<<endl;

    if(d<mind&&d<=maxdist) {
      mini=i;
      mind=d;
    }
  }
  return(mini);
}


int findAdapter(string & s,int minoverlap,float similarity,bool use_hamming){
  unsigned  int L1 = s.size();
  unsigned  int L2 = adapter1.size();
  //  cout << L1 << " " << L2 << endl;
  assert(use_hamming);
  unsigned  int perfect;
  //first half of adapter
  perfect = s.find(adapter1);
  if(perfect<L1)
    return(perfect);

  //second half of adapter
  perfect = s.find(adapter2);
  if(perfect<L1)
    return(perfect-L2);
  
  int a;
  //  cout << "no exact matches, lets try partial matching"<<endl;
  if(use_hamming) a = partial_match(s,adapter1,minoverlap,similarity);
  //  else  a = partial_match(s,lev1,minoverlap,maxdist);

  if(a<(int)L1)
    return(a);

  //second half
  if(use_hamming) a = partial_match(s,adapter2,minoverlap,similarity);
  //  else a = partial_match(s,lev2,minoverlap,maxdist);

  if(a<(int)L1)
    return(a-L2);  

  //ok nothing found. return end of string.
  return(L1);
}

int overlap(string & s1,string & s2,int minoverlap,float similarity) {
  int mini=0,mind,minL;
  if(s1.size()<s2.size())  minL=s1.size();
  else minL=s2.size();
  mind=minL;

  if(minL<minoverlap)
    return(0);

  for(int i=minoverlap;i<minL;i++) {
    int maxdist = ceil( (1. - similarity) * i );
    int d = hamming(s1,s2,0,0,i,maxdist);
    if(d<mind && d<=maxdist) {
      mind=d;
      mini=i;
    }
  }
  return(mini);  
}

//checks for an overlap between r1(suffix) and r2(prefix). if there is overlap create new single read in output return 1
//if no overlap return 0
int matePair::joinReads(fqread & r1,fqread & r2,fqread & output) {
  if(r1.l<minoverlap || r2.l<minoverlap || !joinreads)
    return(0);

  int w = overlap(r1.s,r2.s,minoverlap,similarity);
  if(DEBUG>0)  cout << "w = "<<w<<endl;
  if(w==0)
    return(0);
  else {
    output.h=r1.h;
    output.l3=r1.l3;
    output.s = r1.s + r2.s.substr(w);
    output.q = r1.q + r2.q.substr(w);      
    output.l = output.s.size();
    return(1);
  }  
}

int matePair::resolve_overhang(fqread & r1, fqread & r2,int a,int b) {

  fqread tmp1 = r1.window(b,r1.l);
  fqread tmp2 = r1.window(b,r1.l).rc();

  if(a<minlen) {//preceeding dna is too small.
    pe.r1 = tmp1;
    pe.r2 = r2;
  }
  else if(joinReads(tmp2,r2,se)) {
    mp.r1=r1.window(0,a);
    mp.r2=r2;
  }
  else if( (r1.l-b)>a && !preserve_mp) {
    pe.r1=tmp1;
    pe.r2=r2;
    if(a>=minlen)
      se = r1.window(0,a);
  }
  else {
    if(tmp1.l>=minlen)
      se = tmp1;

    mp.r1=r1.window(0,a);
    mp.r2=r2;
  }

  return(0);
}


matePair::matePair() {
}

matePair::matePair(readPair& readpair,int minovl,float sim,int ml,bool jr,bool uh)
{
  build(readpair,minovl,sim,ml,jr,uh);
}

int matePair::set_preserve_mp(bool p) {
  preserve_mp=p;
  return(0);
}

int matePair::clear() {
  mp.r1.clear();
  mp.r2.clear();
  pe.r1.clear();
  pe.r2.clear();
  unknown.r1.clear();
  unknown.r2.clear();
  se.clear();
  return(0);
}

//this is just a simple routine to trim edges off matepairs where the Nx adapter was not detected
//removes perfect adapter matches on edges that are < minoverlap
int matePair::trimUnknown() {
  //start
  int a1=0,a2=0,b1=unknown.r1.l,b2=unknown.r2.l;
  for(int i=3;i<minoverlap;i++) {
    int offset = adapter2.size()-i;
    int maxdist = floor ( (1.-similarity) * i);
    /*
      cout << maxdist << " " <<  i << " " << hamming(unknown.r1.s,adapter2,0,offset,i,maxdist) << " " << unknown.r1.s.substr(0,i) << " " << adapter2.substr(offset,i) << endl;
      cout << maxdist << " " << i << " " << hamming(unknown.r2.s,adapter2,0,offset,i,maxdist) <<endl;
      cout << maxdist << " " << i << " " << hamming(unknown.r1.s,adapter1,unknown.r1.l-i,0,i,maxdist)  << " " << unknown.r1.s.substr(unknown.r1.l-i,i) << " " << adapter1.substr(0,i) << endl;
      cout << maxdist << " " << i << " " << hamming(unknown.r2.s,adapter1,unknown.r2.l-i,0,i,maxdist) <<endl;
    */
    if(hamming(unknown.r1.s,adapter2,0,offset,i,maxdist)<=maxdist)
      a1=i;
    if(hamming(unknown.r2.s,adapter2,0,offset,i,maxdist)<=maxdist)
      a2=i;
    if(hamming(unknown.r1.s,adapter1,0,offset,i,maxdist)<=maxdist)
      a1=i;
    if(hamming(unknown.r2.s,adapter1,0,offset,i,maxdist)<=maxdist)
      a2=i;
    if(hamming(unknown.r1.s,adapter1,unknown.r1.l-i,0,i,maxdist)<=maxdist)
      b1=unknown.r1.l-i;
    if(hamming(unknown.r2.s,adapter1,unknown.r2.l-i,0,i,maxdist)<=maxdist)
      b2=unknown.r2.l-i;
    if(hamming(unknown.r1.s,adapter2,unknown.r1.l-i,0,i,maxdist)<=maxdist)
      b1=unknown.r1.l-i;
    if(hamming(unknown.r2.s,adapter2,unknown.r2.l-i,0,i,maxdist)<=maxdist)
      b2=unknown.r2.l-i;

  }
  //  cout << "trimUnknown: " << a1 << " " << b1 << " " << a2 << " " << b2 << endl;
  if(a1>0||b1<unknown.r1.l)
    unknown.r1 = unknown.r1.window(a1,b1);
  if(a2>0||b2<unknown.r2.l)
    unknown.r2 = unknown.r2.window(a2,b2);
  return(0);
}

//gets rid of the rare case where external adapters are present (isize < 2*L)
bool matePair::trimExternal(readPair& rp) {
  bool found = false;
  int a,b;
  unsigned int tmp = rp.r1.s.find(r2_external_adapter);//PERFECT MATCH?
  if(tmp>=rp.r1.s.size()) //PARTIAL MATCH?
    a = partial_match(rp.r1.s,r2_external_adapter,minoverlap,similarity);
  else a = (int)tmp;
    
  tmp = rp.r2.s.find(r2_external_adapter);//PERFECT MATCH?
  if(tmp>=rp.r1.s.size()) //PARTIAL MATCH?
    b = partial_match(rp.r2.s,r1_external_adapter,minoverlap,similarity);
  else
    b = (int)tmp;

  //  OK NO ADAPTERS FOUND, LETS TRY LOOKING FOR AN OVERLAP -> PAIRED END FRAG

  if(!(a>0 && a<rp.r1.l)&&!(b>0 && b<rp.r1.l)) {
    fqread rc2 = rp.r2.rc();
    
    int mini=rp.r1.l,mind=rp.r1.l;
    for(int i=0;i<(rp.r1.l-minlen);i++) {
      int compare_length = rp.r2.l-i;
      int maxdist = ceil ( (1.-similarity) * compare_length);
      int d = hamming(rp.r1.s,rc2.s,0,i,rp.r2.l-i,maxdist);
      if(d<mind&&d<=maxdist) {
        mini = i;
        mind = d;
      }
    }
    if(mini<rp.r1.l) {
      a=rp.r1.l-mini;
      b=rp.r2.l-mini;
      //      cout <<"OVERLAP FOUND "<< a <<" " << b<<endl;
    }
  }

  if((a>0 && a<rp.r1.l)||(b>0 && b<rp.r1.l)) {
    if(a<rp.r1.l)
      pe.r1 = rp.r1.window(0,a);
    else
      pe.r1 = rp.r1;
    if(b<rp.r2.l)
      pe.r2 = rp.r2.window(0,b);
    else
      pe.r2 = rp.r2;
    found = true;
  }

  return(found);
}

//aligns s2 to s1 with sim>=sim. returns s1.size() if no alignment found
unsigned int matePair::ham_align(string & s1,string & s2) {
  assert(s1.size()>s2.size());

  int L1 = s1.size();
  int L2 = s2.size();
  assert(L2>=minoverlap);
  int maxd = ceil ( (1.-similarity) * L2);
  int mind=maxd,mini=L1;
  int d;
  //  cout << "L1 = "<<L1<<"\tL2 = "<<L2<<endl;
  for(int i=0;i<(L1-L2);i++) {
    d = hamming(s1,s2,L1-i-L2,0,L2,maxd);

    if(d<mind && d<maxd) {
    // cout << i << " " << L1-i-L2<<endl;
    // cout << s1.substr(L1-i-L2,L2)<<endl;
    // cout << s2 << endl;
    // cout<<d << endl;
      mind=d;
      mini=L1-i;
    }
  }
  //  cout << "mini = "<<mini<<" mind ="<<mind <<endl;
  return(mini);
}

//int matePair::build(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming) {
int matePair::build(readPair& readpair,int minovl,float sim,int ml,bool jr,bool uh) {
  //  assert(readpair.r1.l==readpair.r2.l);
  clear();
  preserve_mp=false;
  minoverlap=minovl;
  similarity=sim;
  minlen=ml;
  joinreads=jr;
  use_hamming=uh;
  
  int L1 = readpair.r1.l;
  int L2 = readpair.r2.l;

  int a1 = findAdapter(readpair.r1.s, minoverlap, similarity,use_hamming);
  int a2 = findAdapter(readpair.r2.s, minoverlap, similarity,use_hamming);
  
  fqread rc1 = readpair.r1.rc();
  fqread rc2 = readpair.r2.rc();

  //check for double adapter.
  if(readpair.r1.s.find(adapter1,a1+adapterj.size())<readpair.r1.s.size()||
     readpair.r1.s.find(adapter2,a1+adapterj.size())<readpair.r1.s.size()||
     readpair.r2.s.find(adapter1,a2+adapterj.size())<readpair.r2.s.size()||
     readpair.r2.s.find(adapter2,a2+adapterj.size())<readpair.r2.s.size())
    return(1);

  int b1 = a1+adapterj.size();
  int b2 = a2+adapterj.size();
  if(DEBUG>1)  cout << a1 <<  " " << b1  <<  " " <<  a2  <<  " " <<  b2 << endl;
  if(a1==L1&&b2<(L2-minlen)) {
    string overhang = rc2.s.substr(0,rc2.l-b2);
    a1 = ham_align(readpair.r1.s,overhang);
    b1 = a1+adapterj.size();
  }

  if(a2==L2&&b1<(L1-minlen)) {
    string overhang = rc1.s.substr(0,rc1.l-b1);
    a2 = ham_align(readpair.r2.s,overhang);
    b2 = a2+adapterj.size();
  }


  if(DEBUG>1)  cout << a1 <<  " " << b1  <<  " " <<  a2  <<  " " <<  b2 << endl;
  if(a1==L1 && a2==L2) {//no adapter found
    if(!joinReads(readpair.r1,rc2,se)) {
      if(!trimExternal(readpair)) {
        unknown=readPair(readpair.r1,readpair.r2);
        trimUnknown();
      }
    }
  }
  else {//adapter found.
    bool both_have_adapter = a1<L1 && a2<L2;
    bool R1_has_adapter_at_end =  a1<L1 && b1>=(L1-minlen);
    bool R2_has_adapter_at_end =  a2<L2 && b2>=(L2-minlen);
    if(a1<L1 && a2<minlen) {//r2 redundant
      if(a1>=minlen) {
        //        pe.r1 = readpair.r1.window(0,a1);  
        //        pe.r2 = readpair.r2.window(b2,b2+a1);
        se = readpair.r1.window(0,a1);  
      }
    } else if(a2<L1 && a1<minlen) {//r1 redundant
      if(a2>=minlen) {
        //        pe.r1 = readpair.r1.window(b1,b1+a2);  
        //        pe.r2 = readpair.r2.window(0,a2);    
        se = readpair.r2.window(0,a2);    
      }
    } else if(both_have_adapter||R1_has_adapter_at_end||R2_has_adapter_at_end) {
      //standard mp
      mp.r1=readpair.r1.window(0,a1);
      mp.r2=readpair.r2.window(0,a2);
      /*
      if((L1-b1)>minlen && b1<=b2)
        se = readpair.r1.window(b1);
      if((L2-b2)>minlen && b2<b1)
        se = readpair.r2.window(b2);
      */
    } 
    else if(b1<L1 && a2==L2) {
      resolve_overhang(readpair.r1,readpair.r2,a1,b1);
    } 
    else if(b2<L2 && a1==L1) {
      resolve_overhang(readpair.r2,readpair.r1,a2,b2);
    }    
  }
  return(0);
}


//LEVENSHTEIN DISTANCE CODE - REMOVED THIS
/*

  levenshtein lev1(adapter1);
  levenshtein lev2(adapter2);

  levenshtein::levenshtein(string s1) {
  L1 = s1.size();  
  this->s1 = s1;
  column = new unsigned int[L1+1];
  }

  int levenshtein::distance(string & s2,int offset,int maxdist,int indel_penalty) {
  assert(s2.size()>=s1.size());
  if(offset<0) 
  int L1 = s1.size()+offset;
  else if(offset>(int)(s2.size()-s1.size()))
  int L1 = s2.size()-offset;
  else
  int L1 = s1.size();

  int L2=L1;

  assert(L1<=s1.size() && L2<=s2.size());
  int offset1 = offset<0 ? -offset : 0;
  offset = offset<0 ? 0:offset;

  for(int j = 1; j <= L1; j++)
  column[j] = j;
  for(int i = 1; i <= L2; i++) {
  column[0] = i;
  //    cout << column[0] << " ";
  lastdiag=i-1;
  for(int j=1; j <= L1; j++) {
  //      cout << column[j] << " ";      
  olddiag = column[j];
  column[j] = MIN3(column[j] + 1, column[j-1] + 1, lastdiag + (s1[offset1+j-1] == s2[offset+i-1] ? 0 : 1));
  lastdiag = olddiag;
  }
  if(i>maxdist && column[i]>maxdist) {
  column[L1]=L1;
  break;
  }
  //    cout << endl;
  }
  return(column[L1]);
  }


  //allows for indel errors in adapter (slower)
  int partial_match(string & s1,levenshtein & lev,int minoverlap,int maxdist) {
  assert(lev.s1.size()<s1.size());
  int mini=-1,mind=lev.s1.size();
  assert((int)s1.size()>=minoverlap);
  
  int start = -(lev.s1.size()-minoverlap);
  int stop = s1.size() - minoverlap;
  for(int i=start;i<stop;i++) {
  int d;
  if(i<0 || i>(int)(s1.size()-lev.s1.size())) d = lev.distance(s1,i,maxdist/2);
  else d = lev.distance(s1,i,maxdist);
  //    cout << i << " " << d << endl;
  if(d<mind) {
  mini=i;
  mind=d;
  }
  }

  if(mind<maxdist)
  return(mini);
  else
  return(s1.size());
  }

*/
