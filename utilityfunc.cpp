#include "utilityfunc.h"

using namespace std;

int argmax(vector<double> & x) {
  int maxind=0;
  for(int i=1;i<(int)x.size();i++) 
    if(x[i]>x[maxind]) maxind=i;
  return(maxind);
}


int which_max(int *x,int n) {
  int maxidx = 0;
  int maxval = x[maxidx];
  for(int i=0;i<n;i++) {
    if(x[i]>maxval) {
      maxidx = i;
      maxval = x[i];
    }
  }
  return maxidx;
}



bool fileexists(string fname){
  ifstream ifile(fname.c_str());
  return ifile.good();
}

void die(const string & err) {
  cerr << "ERROR: "<< err << endl;
  exit(1);
}
