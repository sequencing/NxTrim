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
  try{
    desc.add_options()
      ("help,h", "produce help message")
      ("r1,1", po::value<string>(), "read 1 in fastq format (gzip allowed)")
      ("r2,2", po::value<string>(), "read 2 in fastq format (gzip allowed)")
      ("output-prefix,O", po::value<string>(), "output prefix")
      ("joinreads", "try to merge overhangs from R2 with R1 (default: no joining)")
      //    ("levenshtein", "use Levenshtein distance instead of Hamming distance (slower but possibly more accurate)")
      ("rc", "reverse-complement mate-pair reads (use this if your reads are RF orientation)")
      ("preserve-mp", "preserve MPs even when the corresponding PE has longer reads")
      ("justmp", "creates a mp only library (reads with adapter at the start with be completely N masked)")
      ("similarity", po::value<float>()->default_value(0.85), "The minimum similarity between strings to be considered a match.  Where edit_distance  <=  ceiling( (1-similarity) * string_length )  ")
      ("minoverlap", po::value<int>()->default_value(12), "The minimum overlap to be considered for matching")
      ("minlength", po::value<int>()->default_value(25), "The minimum read length to output (smaller reads will be filtered)");
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  }
  catch(po::unknown_option e){
    cerr << e.what() << std::endl;
    cerr << "Type nxtrim -h for help."<<endl;
    exit(1);
  }
  
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
  bool justmp = opt.count("justmp");
  int minoverlap= opt["minoverlap"].as<int>()-1;//+1 since we use < frequently
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

  if(preserve_mp&&justmp) {
    cerr << "ERROR: the --preserve_mp and --justmp flags are incompatible!" << endl;
    return(1);
  }
  if(justmp)
    preserve_mp=true;

  pairReader infile(r1,r2);


  int nodata=0;
  readPair p;
  pair<int,int> pos;
  matePair m;
  int nweird=0,npass=0,nread=0;
  bool trim_warn=true;
  nxtrimWriter out(prefix,justmp);
  while(infile.getPair(p)) {
    if( p.r1.l!=p.r2.l && trim_warn) {
      cerr << "WARNING: reads with differing lengths. Has this data already been trimmed?"<<endl;
      trim_warn=false;
    }
    if(!p.r1.filtered && !p.r2.filtered) {
      bool weird=m.build(p,minoverlap,similarity,minlen,joinreads,hamming,preserve_mp,justmp);
      nweird+=weird;
      nodata+=(!weird && (m.mp.r1.l==0&&m.mp.r2.l==0)&&  (m.pe.r1.l==0&&m.pe.r2.l==0) &&  (m.unknown.r1.l==0&&m.unknown.r2.l==0) && m.se.l==0);
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
      cout <<  "READ PAIR "<<nread<<endl;
  }
  cout << percent(npass,nread) << "reads passed chastity/purity filters."<<endl;
  cout << percent(nweird,npass) << "reads had TWO copies of adapter (filtered)."<<endl;
  npass-=nweird;
  cout << percent(nodata,npass) << "read pairs were ignored because template length appeared less than read length"<<endl<<endl;
  npass-=nodata;
  cout << percent(out.n_mp,npass) << "read pairs had MP orientation"<<endl;
  cout << percent(out.n_pe,npass) << "read pairs had PE orientation"<<endl;
  cout << percent(out.n_unk,npass) << "had unknown orientation"<<endl;
  cout << percent(out.n_se,npass) << "single end reads were generated"<<endl;
  return(0);
}
