#include "matepair.h"
#include "fastqlib.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

string percent(int num,int den) {
  char buffer[100];
  sprintf(buffer,"%d / %d\t( %.2f%% )\t",num,den,100. * float(num)/float(den));
  return(buffer);
}

int checkParameters(int argc,char **argv,po::variables_map & vm) {
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("r1,1", po::value<string>(), "read 1 in fastq format (gzip allowed)")
    ("r2,2", po::value<string>(), "read 2 in fastq format (gzip allowed)")
    ("output-prefix,O", po::value<string>(), "output prefix")
    ("joinreads", "try to merge overhangs from R2 with R1 (default: no joining)")
    //    ("levenshtein", "use Levenshtein distance instead of Hamming distance (slower but possibly more accurate)")
    ("rc", "reverse-complement mate-pair reads (use this if your reads are RF orientation)")
    ("preserve-mp", "preserve MPs even when the corresponding PE has longer reads")
    ("similarity", po::value<float>()->default_value(0.85), "The minimum similarity between strings to be considered a match.  Where edit_distance  <=  ceiling( (1-similarity) * string_length )  ")
    ("minoverlap", po::value<int>()->default_value(12), "The minimum overlap to be considered for matching")
    ("minlength", po::value<int>()->default_value(21), "The minimum read length to output (smaller reads will be filtered)");

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help") || argc==1) {
      cout << desc << "\n";
      exit(1);
    }

    if (!vm.count("r1") || !vm.count("r2")) {
      cout << "Missing input!"<<endl;
      exit(1);
    }

    if(!vm.count("output-prefix")) {
      cout << "Missing output file!"<<endl;
      exit(1);     
    }
    return(0);
}

int main(int argc,char **argv) {

  po::variables_map opt;
  checkParameters(argc,argv,opt);

  bool joinreads = opt.count("joinreads");
  bool preserve_mp = opt.count("preserve-mp");
  int minoverlap= opt["minoverlap"].as<int>();
  float similarity=opt["similarity"].as<float>();
  int minlen=opt["minlength"].as<int>();
  string r1 = opt["r1"].as<string>();
  string r2 = opt["r2"].as<string>();
  string prefix = opt["output-prefix"].as<string>();
  //  bool hamming = !opt.count("levenshtein");
  bool hamming = true;
  bool rc = opt.count("rc");
  cout << "Trimming:\nR1:\t" <<r1<<"\nR2:\t"<<r2<<endl;
  cout << "Output: " << prefix <<".*.fastq.gz"<<endl;
  if(hamming)  cout << "Using Hamming distance"<<endl;
  else cout << "Using Levenshtein distance\n"<<endl;
  if(preserve_mp) cout<< "--preserve-mp is on: will favour MPs over PEs" <<endl;
  if(joinreads) cout<< "--joinreads is on: will attempt to merge R1 with R2 that proceeds an adapter" <<endl;

  pairReader infile(r1,r2);
  fastqWriter mp_out(prefix+".mp.fastq.gz");
  fastqWriter pe_out(prefix+".pe.fastq.gz");
  fastqWriter se_out(prefix+".se.fastq.gz");
  fastqWriter unknown_out(prefix+".unknown.fastq.gz");

  readPair p;
  pair<int,int> pos;
  int npass=0;
  int nread=0;  
  int nweird=0;
  matePair m;
  int n_mp=0,n_pe=0,n_se=0,n_unk=0;
  bool trim_warn=true;
  while(infile.getPair(p)) {
    if( p.r1.l!=p.r2.l && trim_warn) {
      cerr << "WARNING: reads with differing lengths. Has this data already been trimmed?"<<endl;
      trim_warn=false;
      // p.r1.print();
      // p.r2.print();
      // exit(1);
    }
    if(!p.r1.filtered && !p.r2.filtered) {
      nweird+=m.build(p,minoverlap,similarity,minlen,joinreads,hamming,preserve_mp);
  
      if(rc) {
	n_pe+=pe_out.write(m.pe);
	m.mp.rc();
	n_mp+=mp_out.write(m.mp);
	m.unknown.rc();
	n_unk+=unknown_out.write(m.unknown);
      }
      else{
	m.pe.rc();
	n_pe+=pe_out.write(m.pe);
	n_mp+=mp_out.write(m.mp);
	n_unk+=unknown_out.write(m.unknown);
      }

      n_se+=se_out.write(m.se);
      npass++;
    }
    nread++;
    if(nread%10000==0)
      cout <<  "READ PAIR "<<nread<<endl;
  }
  cout << percent(npass,nread) << "reads passed chastity/purity filters."<<endl;
  cout << percent(nweird,npass) << "reads had TWO copies of adapter (filtered)."<<endl<<endl;
  cout << percent(n_mp,npass) << "read pairs had MP orientation"<<endl;
  cout << percent(n_pe,npass) << "read pairs had PE orientation"<<endl;
  cout << percent(n_se,npass) << "single end reads were generated"<<endl;
  cout << percent(n_unk,npass) << "had unknown orientation"<<endl;
  return(0);
}
