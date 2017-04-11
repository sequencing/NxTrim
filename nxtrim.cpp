#include "version.h"
#include "matepair.h"
#include "fastqlib.h"
#include <getopt.h>

using namespace std;

string percent(int num,int den) {
  char buffer[100];
  sprintf(buffer,"%d / %d\t( %.2f%% )\t",num,den,100. * float(num)/float(den));
  return(buffer);
}

void usage() {
  cerr << "\nProgram:\tnxtrim" << endl;
  cerr << "Version:\t" << VERSION <<endl;
  cerr << "Contact:\tjoconnell@illumina.com\n" << endl;
  cerr << "Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.\n"<<endl;
  cerr << "Usage:\tnxtrim -1 R1.fastq.gz -2 R2.fastq.gz [options]\n" << endl;
  cerr << "Required arguments:"    <<endl;
  cerr << "  -1 [ --r1 ] arg                 read 1 in fastq format (gzip allowed)"<<endl;
  cerr << "  -2 [ --r2 ] arg                 read 2 in fastq format (gzip allowed)"<<endl;
  cerr << "Allowed options:"<<endl;
  cerr << "  -O [ --output-prefix ] arg      output prefix"<<endl;
  cerr << "  --justmp                        just creates a the mp/unknown libraries (reads with adapter at the start will be completely N masked)"<<endl;
  cerr << "  --stdout                        print trimmed reads to stdout (equivalent to justmp)"<<endl;
  cerr << "  --stdout-mp                     print only known MP reads to stdout (good for scaffolding)"<<endl;
  cerr << "  --stdout-un                     print only unknown reads to stdout"<<endl;
  cerr << "  --joinreads                     try to merge overhangs from R2 with R1 (default: no joining)"<<endl;
  cerr << "  --rf                            leave mate pair reads in RF orientation [by default are flipped into FR]"<<endl;
  cerr << "  --preserve-mp                   preserve MPs even when the corresponding PE has longer reads"<<endl;
  cerr << "  --ignorePF                      ignore chastity/purity filters in read headers"<<endl;
  cerr << "  --separate                      output paired reads in separate files (prefix_R1/prefix_r2). Default is interleaved."<<endl;
  cerr << "  -s, --similarity arg (=0.85)    The minimum similarity between strings to be considered a match.  Where hamming_distance  <=  ceiling( (1-similarity) * string_length )"<<endl;
  cerr << "  -v, --minoverlap arg (=12       The minimum overlap to be considered for matching"<<endl;
  cerr << "  -l, --minlength arg (=21)       The minimum read length to output (smaller reads will be filtered)"<<endl;
  cerr << "  -w, --smith-waterman            Use Smith-Waterman alignmnent rather than simple Hamming matching"<<endl;
  cerr << "  -x, --smith-waterman-score      minimum score for an SW match (higher: less sensivity, lower: mor sensitivity))"<<endl;    
  exit(0);
}

#define STDOUT 0
#define JUSTMP 1
#define JOINREADS 2
#define NORC 3
#define PMP 4
#define IGNOREPF 5
#define MP 6
#define UNKNOWN 7
#define SEPARATE 8
#define STDOUT_MP 9
#define STDOUT_UN 10

int main(int argc,char **argv) {
  int c;
  if(argc<2)
    usage();

  bool joinreads=false;
  bool preserve_mp=false;
  bool justmp=false;
  int minoverlap=12;
  float similarity=0.85;
  int minlen=21;
  char *r1 = NULL;
  char *r2 = NULL;
  string prefix;
  bool rc = true;
  bool ignorePF = false;
  bool write_stdout=false;
  bool write_stdout_mp=false;
  bool write_stdout_un=false;
  bool hamming = true;
  bool separate=false;
  float sw_score = 12;
  static struct option loptions[] =    {
    {"r1",1,0,'1'},	
    {"r2",1,0,'2'},	
    {"output-prefix",1,0,'O'},	
    {"stdout",0,0,STDOUT},
    {"stdout-mp",0,0,STDOUT_MP},
    {"stdout-un",0,0,STDOUT_UN},
    {"justmp",0,0,JUSTMP},
    {"joinreads",0,0,JOINREADS},
    {"rf",0,0,NORC},
    {"preserve-mp",0,0,PMP},
    {"ignorePF",0,0,IGNOREPF},
    {"mp",0,0,MP},
    {"unknown",0,0,UNKNOWN},
    {"separate",0,0,SEPARATE},
    {"smith-waterman",0,0,'w'},    
    {"similarity",1,0,'s'},
    {"minoverlap",1,0,'v'},
    {"minlength",1,0,'l'},
    {0,0,0,0}
  };
  while ((c = getopt_long(argc, argv, "1:2:O:s:v:l:x:w",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case '1': r1 = optarg; break;
      case '2': r2 = optarg; break;
      case 'O': prefix = optarg; break;
      case 's': similarity = atof(optarg); break;    
      case 'v': minoverlap = atoi(optarg); break;    
      case 'l': minlen = atoi(optarg); break;    
      case STDOUT: write_stdout=true; break;    
      case STDOUT_MP: write_stdout_mp=true; break;    
      case STDOUT_UN: write_stdout_un=true; break;    
      case JUSTMP:justmp=true; break;    
      case JOINREADS:joinreads=true; break;    
      case NORC:rc=false; break;    
      case PMP:preserve_mp=true; break;    
      case IGNOREPF:ignorePF=true; break;    
      case SEPARATE:separate=true; break;
      case 'w':hamming=false; break;
      case 'x':sw_score=atof(optarg); break;
      default: die("Unrecognised argument");
      }
  }
  if(!(r1==NULL&&r2==NULL) && !(r1!=NULL&&r2!=NULL))
    die("both --r1 and --r2 must be speicified");
  if(write_stdout && !prefix.empty())
    die("--stdout and -O are incompatible");
  if(!write_stdout && !write_stdout_mp && !write_stdout_un && prefix.empty() )
    die("one of --stdout / --stdout-mp / --stdout-un / -O must be specified");
  if(preserve_mp&&justmp) 
    die("the --preserve_mp and --justmp flags are incompatible!");
  if( (write_stdout+write_stdout_mp+write_stdout_un)>1)
    die("only one of --stdout / --stdout-mp / --stdout-un may be specified!");
    
  if(write_stdout||write_stdout_mp||write_stdout_un)
  {
      cerr << "Writing to stdout"<<endl;
      justmp=true;
  }
  else
  {
      cerr << "Output: " << prefix <<".*.fastq.gz"<<endl;
  }

  cerr << "Trimming:\nR1:\t" <<r1<<"\nR2:\t"<<r2<<endl;
  if(preserve_mp) cerr<< "--preserve-mp is on: will favour MPs over PEs" <<endl;
  if(joinreads) cerr<< "--joinreads is on: will attempt to merge R1 with R2 that proceeds an adapter" <<endl;


  if(justmp)
    preserve_mp=true;

  pairReader infile(r1,r2);

  int nodata=0;
  readPair p;
  pair<int,int> pos;
  matePair m;
  int nweird=0,npass=0,nread=0;
  bool trim_warn=true;

  nxtrimWriter out;
  if(write_stdout||write_stdout_un||write_stdout_mp)  {
    out.open();
    if(write_stdout_un) out.setMP(false);
    if(write_stdout_mp) out.setUN(false);
  }
  else  
    out.open(prefix,justmp,separate);
  if(!hamming)
  {
      similarity=sw_score;
  }
  int se_only = 0;
  while(infile.next(p))
  {

    if(p.r1.l!=p.r2.l && trim_warn)
    {
      cerr << "WARNING: reads with differing lengths. Has this data already been trimmed?"<<endl;
      trim_warn=false;
    }
    if((!p.r1.filtered && !p.r2.filtered)||ignorePF)
    {
      bool weird=m.build(p,minoverlap,similarity,minlen,joinreads,hamming,preserve_mp,justmp);
      nweird+=weird;
      if(!weird) {
	nodata+=( (m.mp.r1.l==0||m.mp.r2.l==0) && (m.pe.r1.l==0||m.pe.r2.l==0) && (m.unknown.r1.l==0||m.unknown.r2.l==0) && m.se.l==0);
	se_only += ( (m.mp.r1.l==0||m.mp.r2.l==0) && (m.pe.r1.l==0||m.pe.r2.l==0) && (m.unknown.r1.l==0||m.unknown.r2.l==0) ) && m.se.l>0;
      }
      if(rc) {
	m.mp.rc();
	m.unknown.rc();
      }
      else{
	m.pe.rc();
      }
      out.write(m);
      npass++;
    }
    nread++;
    if(nread%10000==0)
      cerr <<  "READ PAIR "<<nread<<endl;

  }
  cerr << "\nTrimming summary:"<<endl;
  cerr << percent(npass,nread) << "reads passed chastity/purity filters."<<endl;
  cerr << percent(nweird,npass) << "reads had TWO copies of adapter (filtered)."<<endl;
  npass-=nweird;
  cerr << percent(nodata,npass) << "read pairs were ignored because template length appeared less than read length"<<endl;
  npass-=nodata;
  cerr << npass << " remaining reads were trimmed"<<endl<<endl;
  cerr << percent(out.n_mp,npass) << "read pairs had MP orientation"<<endl;
  cerr << percent(out.n_pe,npass) << "read pairs had PE orientation"<<endl;
  cerr << percent(out.n_unk,npass) << "read pairs had unknown orientation"<<endl;
  cerr << percent(se_only,npass) << "were single end reads"<<endl<<endl;
  cerr << percent(out.n_se - se_only,npass) << "extra single end reads were generated from overhangs"<<endl;
  return(0);
}
