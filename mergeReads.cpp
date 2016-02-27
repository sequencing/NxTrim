#include "version.h"
#include "githash.h"
#include "fastqlib.h"
#include <getopt.h>

using namespace std;

string percent(int num,int den) {
  char buffer[100];
  sprintf(buffer,"%d / %d\t( %.2f %% )\t",num,den,100. * float(num)/float(den));
  return(buffer);
}
void die(const string & err) {
  cerr << "ERROR: "<< err << endl;
  exit(1);
}

void usage() {
  cerr << "\nProgram:\tmergeReads: Simple utility for creating an interleaved fastq file from two separate R1/R2 files."<<endl;
  cerr << "Version:\t" << VERSION <<endl;
  cerr << "Contact:\tjoconnell@illumina.com\n" << endl;
  cerr << "Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.\n"<<endl;
  cerr << "Usage:\tmergeReads -1 R1.fastq.gz -2 R2.fastq.gz [options]\n" << endl;
  cerr << "Required arguments:"    <<endl;
  cerr << "  -1 [ --r1 ] arg                 read 1 in fastq format (gzip allowed)"<<endl;
  cerr << "  -2 [ --r2 ] arg                 read 2 in fastq format (gzip allowed)"<<endl;
  cerr << "Allowed options:"<<endl;
  cerr << "  -O [ --output-prefix ] arg      output prefix"<<endl;
  //  cerr << "  --stdout                        print trimmed reads to stdout"<<endl;
  cerr << "  --rc                          reverse-complement mate-pair reads (use this if your reads are already in FR orientation)"<<endl;
  exit(0);
}

int main(int argc,char **argv) {
  if(argc<2) usage();
  int c;
  string r1,r2,prefix;
  bool rc = false;

  static struct option loptions[] =    {
    {"r1",1,0,'1'},	
    {"r2",1,0,'2'},	
    {"output-prefix",1,0,'O'},
    {"norc",0,0,'r'},
    {0,0,0,0}
  };

  while ((c = getopt_long(argc, argv, "1:2:O:s:v:l:",loptions,NULL)) >= 0) {  
    switch (c)      {
      case '1': r1 = optarg; break;
      case '2': r2 = optarg; break;
      case 'O': prefix = optarg; break;
      case 'r':rc=true;
      }
  }
  if(r1.empty()||r2.empty())
    die("both --r1 and --r2 must be speicified");
  if(prefix.empty())
    die("the -O option is required");
      
  cerr << "Merging:\nR1:\t" <<r1<<"\nR2:\t"<<r2<<endl;
  cerr << "Output: " << prefix <<".fastq.gz"<<endl;
  if(rc)
    cerr << "Reads will be reverse-complemented."<<endl;

  pairReader infile(r1,r2);
  fastqWriter outfile(prefix+".fastq.gz");

  readPair p;
  int npass=0;
  int nread=0;  

  while(infile.next(p)) {
    if(!p.filtered) {
      if(rc) {
	p.rc();
	outfile.write(p);
      }
      else outfile.write(p);
      npass++;
    }
    nread++;
    if(nread%10000==0)
      cerr <<  "READ PAIR "<<nread<<endl;
  }
  cerr << percent(npass,nread) << "reads passed chastity/purity filters."<<endl;
  return(0);
}
