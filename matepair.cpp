#include "matepair.h"

//nextera mp adapters
string adapter1 = "CTGTCTCTTATACACATCT";
string adapter2 = "AGATGTGTATAAGAGACAG";
string adapterj = adapter1+adapter2;
//EXTERNAL adapters
// string r1_external_adapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC";
// string r2_external_adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATC";                
string r1_external_adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
string r2_external_adapter = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

#define DEBUG 0

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))


//only handles substitution errors in adapter (faster)
int partial_match(string & s1,string & s2,int minoverlap,float similarity) {
  //  assert(s2.size()<s1.size());
  if(s2.size()>=s1.size())
    return(s1.size());
  int mini=s1.size(),mind=s2.size();
  assert((int)s1.size()>=minoverlap);
  
  int start = -(s2.size()-minoverlap);
  int stop = s1.size() - minoverlap;
  //  cout << "Range: " << start << " " << stop << endl;
  for(int i=start;i<stop;i++) {
    float compare_length = s2.size();
    if( i<0) compare_length += i;
    if( (int)s2.size() > ((int)s1.size()-i)) compare_length = (float)(s1.size()-i);
    int maxdist = ceil ( (1.-similarity) * compare_length);
    int d = hamming(s1,s2,i,0,s2.size(),maxdist);
    //    cout << i<<" "<<maxdist<<" "<<d<<endl;

    if(d<mind&&d<maxdist) {
      mini=i;
      mind=d;
      if(DEBUG>2) {
	cout << compare_length << endl;
	cout << mini << " " << mind << "<" << mind << endl;
      }
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

//returns min( hamming( s1[offset1,offset+L], s2[offset2,offset2+L] ) , maxd )
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

int overlap(string & s1,string & s2,int minoverlap,float similarity) {
  int mini=0,mind,minL;
  if(s1.size()<s2.size())  minL=s1.size();
  else minL=s2.size();
  mind=minL;

  if(minL<minoverlap)
    return(0);

  for(int i=minoverlap;i<minL;i++) {
    int maxdist = ceil( (1. - similarity) * i );
    int d = hamming(s1,s2,s1.size()-i,0,i,maxdist);
    if(d<mind && d<maxdist) {
      mind=d;
      mini=i;
    }
  }
  if(DEBUG>1) cout << "mind = "<<mind<<"\tmini = "<<mini<<endl;
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
    int offset = r1.s.size() - w;
    for(int i=0;i<w;i++) {//takes highest quality base.
      int j = offset+w;
      if( (uint8_t)r1.q[j] <  (uint8_t)r2.q[i] ) {
	output.s[j] = r2.s[i];
	output.q[j] = r2.q[i];
      }
    }
    output.l = output.s.size();
    return(1);
  }  
}

int matePair::resolve_overhang(fqread & r1, fqread & r2,int a,int b) {
  if(DEBUG>0)  cout << "Resolving overhang"<<endl;
  fqread tmp1 = r1.window(b,r1.l);
  fqread tmp2 = r1.window(b,r1.l).rc();
  if(DEBUG>1) {
    cout << r2.s <<endl;;
    cout << tmp2.s << endl;
  }
  if(a<minlen) {//preceeding dna is too small. 
    //TODO:could possibly merge for a big read here
    if(justmp) {
      mp.r1 = r1.mask();
      mp.r2 = r2;
    }
    else {
      pe.r1 = tmp1;
      pe.r2 = r2;
    }
  }
  else if(joinReads(r2,tmp2,mp.r2)) {
    mp.r1=r1.window(0,a);
  }
  else if(r1.notN(b,r1.l)>r1.notN(0,a) && !preserve_mp) {
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

matePair::matePair(readPair& readpair,int minovl,float sim,int ml,bool jr,bool uh,bool pmp,bool jmp)
{
  build(readpair,minovl,sim,ml,jr,uh,pmp,jmp);
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
  for(int i=3;i<=minoverlap;i++) {
    int offset = unknown.r1.l-i;    
    int maxd = ceil ( (1.-similarity) * i);
    if(hamming(unknown.r1.s,adapter1,offset,0,i,maxd)<maxd)
      a1=offset;
    offset = unknown.r2.l-i;    
    if(hamming(unknown.r2.s,adapter1,offset,0,i,maxd)<maxd)
      a2=offset;
  }

  if(DEBUG>0)
    if(a1>0||a2>0)
      cout << "trimUnknown: " << a1 << " " << b1 << " " << a2 << " " << b2 << endl;

  if(a1>0)
    unknown.r1 = unknown.r1.window(0,a1);
  if(a2>0)
    unknown.r2 = unknown.r2.window(0,a2);
  return(0);
}

//gets rid of the rare case where external adapters are present (isize < 2*L)
bool matePair::trimExternal(readPair& rp) {
  bool found = false;
  int a,b;

  unsigned int tmp = rp.r1.s.find(r1_external_adapter);//PERFECT MATCH?
  if(tmp>=rp.r1.s.size()) //PARTIAL MATCH?
    a = partial_match(rp.r1.s,r1_external_adapter,minoverlap,similarity);
  else a = (int)tmp;
    
  tmp = rp.r2.s.find(r2_external_adapter);//PERFECT MATCH?
  if(tmp>=rp.r1.s.size()) //PARTIAL MATCH?
    b = partial_match(rp.r2.s,r2_external_adapter,minoverlap,similarity);
  else
    b = (int)tmp;

  if(DEBUG>1) {
    if((a>0 && a<rp.r1.l)||(b>0 && b<rp.r1.l)) {
      cout << "EXTERNAL ADAPTER DETECTED " << a << " " << b << endl;
      if(a>0 && a<rp.r1.l) {
	rp.r1.window(a,rp.r1.l).print();
      }
      if(b>0 && b<rp.r1.l) {
	rp.r2.window(b,rp.r2.l).print();
      }
      rp.r1.print();
      rp.r2.print();      
    }
  }

  //  OK NO ADAPTERS FOUND, LETS TRY LOOKING FOR AN OVERLAP -> PAIRED END FRAG
  if(!(a>0 && a<rp.r1.l)&&!(b>0 && b<rp.r1.l)) {
    fqread rc2 = rp.r2.rc();
    
    int mini=rp.r1.l,mind=rp.r1.l;
    for(int i=0;i<(rp.r1.l-minlen);i++) {
      int compare_length = rp.r2.l-i;
      int maxdist = ceil ( (1.-similarity) * compare_length);
      int d = hamming(rp.r1.s,rc2.s,0,i,rp.r2.l-i,maxdist);
      if(d<mind&&d<maxdist) {
        mini = i;
        mind = d;
      }
    }
    if(mini<rp.r1.l) {
      a=rp.r1.l-mini;
      b=rp.r2.l-mini;
      if(DEBUG>1)      cout <<"OVERLAP FOUND "<< a <<" " << b<<endl;
    }
  }

  if((a>0 && a<rp.r1.l)||(b>0 && b<rp.r1.l)) {
    found = true;
    if(justmp) {
      if(a<rp.r1.l)
	mp.r1 = rp.r1.mask();
      else
	mp.r2 = rp.r2;
      if(a<rp.r2.l)
	mp.r2 = rp.r2.mask();
      else
	mp.r2 = rp.r1;
    }
    else {
      if(a<rp.r1.l)
	pe.r1 = rp.r1.window(0,a);
      else
	pe.r1 = rp.r1;
      if(b<rp.r2.l)
	pe.r2 = rp.r2.window(0,b);
      else
	pe.r2 = rp.r2;
    }
  }
  return(found);
}

//aligns s2 to s1 with sim>=sim. returns s1.size() if no alignment found
unsigned int matePair::ham_align(string & s1,string & s2) {
  if(s1.size()<s2.size())
    return(s1.size());

  int L1 = s1.size();
  int L2 = s2.size();
  assert(L2>=minoverlap);
  int maxd = ceil ( (1.-similarity) * L2);
  int mind=maxd,mini=L1;
  int d;
  for(int i=0;i<(L1-L2);i++) {
    d = hamming(s1,s2,L1-i-L2,0,L2,maxd);

    if(d<mind) {
      mind=d;
      mini=L1-i;
    }
  }
  if(d>maxd)//hit wasnt good enough
    mini=L1;
  return(mini);
}


//checks the right end of a read for partial adapter hit
int checkRight(string & s1,string & adapter,int offset,int minoverlap,float similarity) {
  assert(offset < (s1.size()-minoverlap));
  int a=s1.size();
  int mind = s1.size();
  for(int i=offset;i<(s1.size()-minoverlap);i++) {
    int compare_len = (s1.size() - i);
    int maxdist = ceil(compare_len * (1. - similarity));
    int d = hamming(s1,adapter,i,0,compare_len,maxdist);
    if(d<mind&&d<maxdist) {
      a=i;
      mind=d;
    }
  }
  return(a);    
}

//int matePair::build(readPair & readpair,int minoverlap,float similarity,int minlen,bool joinreads,bool use_hamming) {
int matePair::build(readPair& readpair,int minovl,float sim,int ml,bool jr,bool uh,bool pmp,bool jmp) {
  //  assert(readpair.r1.l==readpair.r2.l);
 clear();
  justmp=jmp;
  preserve_mp=pmp;
  minoverlap=minovl;
  similarity=sim;
  minlen=ml;
  joinreads=jr;
  use_hamming=uh;
  
  int L1 = readpair.r1.l;
  int L2 = readpair.r2.l;
  if(L1<minlen||L2<minlen) {
    cout << "readlength < minlenght"<<endl;
    return(0);
  }

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
  if(DEBUG>1) {
    cout << "L1 ="<<L1<<endl;
    cout << "L2 ="<<L2<<endl;

    cout << a1 <<  " " << b1  <<  " " <<  a2  <<  " " <<  b2 << endl;
  }


  if(a1==L1&&b2<(L2-minoverlap)) {//try to overlang the r2 overhang to r1 -> finds adapter on r1
    string overhang = rc2.s.substr(0,rc2.l-b2);
    a1 = ham_align(readpair.r1.s,overhang);
    b1 = a1+adapterj.size();
  }

  if(a2==L2&&b1<(L1-minoverlap)) {//vice-versa
    string overhang = rc1.s.substr(0,rc1.l-b1);
    a2 = ham_align(readpair.r2.s,overhang);
    b2 = a2+adapterj.size();
  }
  if(DEBUG>1)  cout << a1 <<  " " << b1  <<  " " <<  a2  <<  " " <<  b2 << endl;
  int minoverlap2 = 1; //final attempt to find unidfentifed adapters
  if(a1<L1&&a2==L2)//we know R1 has adapter. try check R2 for adapter with more liberal thresholds
    a2 = checkRight(readpair.r2.s,adapter1, L2-minoverlap, minoverlap2, similarity);
  if(a2<L2&&a1==L1)//vice-versa
    a1 = checkRight(readpair.r1.s,adapter1, L1-minoverlap, minoverlap2, similarity);

  if(DEBUG>1)  cout << a1 <<  " " << b1  <<  " " <<  a2  <<  " " <<  b2 << endl;

  if(a1==L1 && a2==L2) {//no adapter found
    // we could potentially run if(!joinReads(readpair.r1,rc2,se)) but this tends to give a lot of false joins
    // possible improvement: check for r1/r2 overlap in absence of adapter -> overlap implies PE
    if(!trimExternal(readpair)) {
      unknown=readPair(readpair.r1,readpair.r2);
      trimUnknown();
    }    
    if(DEBUG>1) cout << "CASE A"<<endl;
  }
  else {//adapter found.
    bool both_have_adapter = a1<L1 && a2<L2;
    bool R1_has_adapter_at_end =  a1<L1 && b1>=(L1-minlen);
    bool R2_has_adapter_at_end =  a2<L2 && b2>=(L2-minlen);
    if(a1<minlen&&a2<minlen) {//very short template. discard
      return(0);
    }
    else if(a1<(L1-minoverlap) && a2<minlen) {//r2 redundant
      if(justmp) 
	mp=readPair(readpair.r1.window(0,a1),readpair.r2.mask()) ;
      else
	se = readpair.r1.window(0,a1); 
      if(DEBUG>1) cout << "CASE B"<<endl;
    }
    else if(a2<(L2-minoverlap) && a1<minlen) {//r1 redundant
      if(justmp)
	mp=readPair(readpair.r1.mask(),readpair.r2.window(0,a2));
      else
	se = readpair.r2.window(0,a2);    
      if(DEBUG>1) cout << "CASE C"<<endl;
    }
    else if(a1>=(L1-minoverlap) && a2<minlen) {//obvious PE
      if(a1>=minlen && (L2-b2)>=minlen) {
	if(justmp) {
	  mp=readPair(readpair.r1.window(0,a1),readpair.r2.mask()) ;
	}
	else {
	  pe.r1 = readpair.r1.window(0,a1);  
	  pe.r2 = readpair.r2.window(b2,b2+a1);
	}        
      }
      if(DEBUG>1) cout << "CASE D"<<endl;
    } 
    else if(a2>=(L2-minoverlap) && a1<minlen) {//obvious PE
      if(a2>=minlen && (L1-b1)>=minlen) {
        if(justmp){
	  mp=readPair(readpair.r1.mask(),readpair.r2.window(0,a2));
	}
	else{
	  pe.r1 = readpair.r1.window(b1,b1+a2);  
	  pe.r2 = readpair.r2.window(0,a2);    
	}
      }
      if(DEBUG>1) cout << "CASE E"<<endl;
    } 
    else if(both_have_adapter||R1_has_adapter_at_end||R2_has_adapter_at_end) {
      //standard mp
      mp.r1=readpair.r1.window(0,a1);
      mp.r2=readpair.r2.window(0,a2);
      if(DEBUG>1) cout << "CASE F"<<endl;
      /*
      if((L1-b1)>minlen && b1<=b2)
        se = readpair.r1.window(b1);
      if((L2-b2)>minlen && b2<b1)
        se = readpair.r2.window(b2);
      */
    } 
    else if(b1<L1 && a2==L2) {
      resolve_overhang(readpair.r1,readpair.r2,a1,b1);
      if(DEBUG>1) cout << "CASE G"<<endl;
    } 
    else if(b2<L2 && a1==L1) {
      resolve_overhang(readpair.r2,readpair.r1,a2,b2);
      fqread swap1 = pe.r1;
      pe.r1 = pe.r2;
      pe.r2 = swap1;
      fqread swap2 = mp.r1;
      mp.r1 = mp.r2;
      mp.r2 = swap2;
      if(DEBUG>1) cout << "CASE H"<<endl;
    }    
  }
  return(0);
}

nxtrimWriter::nxtrimWriter(string prefix,bool jmp,bool separate) {

  n_mp=0;
  n_pe=0;
  n_se=0;
  n_unk=0;
  justmp = jmp;

  if(separate)     mp_out.open(prefix+"_R1.mp.fastq.gz", prefix+"_R2.mp.fastq.gz");
  else     mp_out.open(prefix+".mp.fastq.gz");

  if(separate)    unknown_out.open(prefix+"_R1.unknown.fastq.gz",prefix+"_R2.unknown.fastq.gz");
  else unknown_out.open(prefix+".unknown.fastq.gz");  

  if(!justmp) {
    if(separate) pe_out.open(prefix+"_R1.pe.fastq.gz",prefix+"_R2.pe.fastq.gz");
    else  pe_out.open(prefix+".pe.fastq.gz");
    se_out.open(prefix+".se.fastq.gz");
  }
}

int nxtrimWriter::write(matePair m) {
  n_mp+=mp_out.write(m.mp);
  n_unk+=unknown_out.write(m.unknown);  
  if(!justmp) {
    n_pe+=pe_out.write(m.pe);
    n_se+=se_out.write(m.se);  
  }
  return(0);
}
