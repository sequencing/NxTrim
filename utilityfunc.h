#pragma once
#include <utility> 
#include <vector>
#include <bitset>
//#include <unordered_set>
#include <sys/stat.h>
#include <sstream> 
#include <fstream>
#include <iostream>
#include <assert.h> 
#include <stdlib.h>

using namespace std;

int which_max(int *x,int n);

//unordered_set<int> sampleIndex(int k,int n);

template <class T>
bool comparator ( const pair<T,int>& l, const pair<T,int>& r)  { return l.first < r.first; }

template <class T>
vector<int> argsort(vector<T> *list) {

  vector<pair<T,int> > tosort(list->size());
  for(int i=0;i<list->size();i++) {
    tosort[i] = make_pair( (*list)[i] ,i );
  }
  sort(tosort.begin(),tosort.end());//,comparator);
  vector<int> ret(tosort.size());
  for(int i=0;i<tosort.size();i++) {
    ret[i] = tosort[i].second;
    //    cout << tosort[i].first << "\t" << tosort[i].second << endl;
  }
  return(ret);
}

template <class T>
T **newMatrix(int nrow,int ncol) {
  T **ret = new T*[nrow];
  for(int i=0;i<nrow;i++)
    ret[i] = new T[ncol];  
  return(ret);
}


template <class T>
void delMatrix(T **mat,int nrow,int ncol) {
  for(int i=0;i<nrow;i++)
    delete[] mat[i];
  delete[] mat;
}

template <class T>
void printMatrix(T **H,int nrow,int ncol) {
  for(int i=0;i<nrow;i++) {
    for(int j=0;j<ncol;j++)
      cout << (unsigned int) H[i][j] << "\t";
    cout << endl;
  }
}


void die(const string & err);

bool fileexists(string fname);

int argmax(vector<double> & x);

