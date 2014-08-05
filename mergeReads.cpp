#include "fastqlib.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

string percent(int num,int den) {
  char buffer[100];
  sprintf(buffer,"%d / %d\t( %.2f %% )\t",num,den,100. * float(num)/float(den));
  return(buffer);
}

int checkParameters(int argc,char **argv,po::variables_map & vm) {
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("r1,1", po::value<string>(), "read 1 in fastq format (gzip allowed)")
    ("r2,2", po::value<string>(), "read 2 in fastq format (gzip allowed)")
    ("output-prefix,O", po::value<string>(), "output prefix")
    ("rc", "reverse-complement  reads");

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help") || argc==1) {
      cout << desc << "\n";
      exit(1);
    }

    if (!vm.count("r1") || !vm.count("r2") ||  !vm.count("output-prefix")) {
      cout << "Missing input!"<<endl;
      exit(1);
    }

    return(0);
}

int main(int argc,char **argv) {

  po::variables_map opt;
  checkParameters(argc,argv,opt);

  string r1 = opt["r1"].as<string>();
  string r2 = opt["r2"].as<string>();
  string prefix = opt["output-prefix"].as<string>();
  bool rc = opt.count("rc");
  cout << "Merging:\nR1:\t" <<r1<<"\nR2:\t"<<r2<<endl;
  cout << "Output: " << prefix <<".fastq.gz"<<endl;
  if(rc)
    cout << "Reads will be reverse-complemented."<<endl;

  pairReader infile(r1,r2);
  fastqWriter outfile(prefix+".fastq.gz");

  readPair p;
  int npass=0;
  int nread=0;  

  while(infile.getPair(p)) {
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
      cout <<  "READ PAIR "<<nread<<endl;
  }
  cout << percent(npass,nread) << "reads passed chastity/purity filters."<<endl;
  return(0);
}
