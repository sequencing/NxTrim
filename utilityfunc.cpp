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

#include "utilityfunc.h"

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


ifile::ifile() {
}

ifile::ifile(string filename , bool binary) {
  open(filename, binary);
}

ifile::~ifile() {
	close();
}

string ifile::name() {
	return file;
}

bool ifile::open(string filename, bool binary) {
	file = filename;
	string ext = filename.substr(filename.length()-2);

	if (ext == "gz") {
		fd.open(file.c_str(), ios::in | ios::binary);
		push(io::gzip_decompressor());
	// } else if (ext == "bz2") {
	// 	fd.open(file.c_str(), ios::in | ios::binary);
	// 	push(io::bzip2_decompressor());
	} else if (binary) {
		fd.open(file.c_str(), ios::in | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

bool ifile::readString(string & str) {
	int s;
	if (!read((char*)&s, sizeof(int))) return false;
	char  * buffer = new char [s + 1];
	if (!read(buffer, s)) return false;
	buffer[s] = '\0';
	str = string(buffer);
	delete [] buffer;
	return true;
}

void ifile::close() {
	 if (!empty()) reset();
	 fd.close();
}

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
ofile::ofile() {
}

ofile::ofile(string filename , bool binary) {
	open(filename, binary);
}

ofile::~ofile() {
	close();
}

string ofile::name() {
	return file;
}

bool ofile::open(string filename, bool binary) {
  if(ifstream(filename.c_str())) {
    cerr << filename << " exists! Will not overwrite.\nExiting..." << endl;
    exit(1);
  }
  file = filename;
  string ext = filename.substr(filename.length()-2);
	if (ext == "gz") {
		fd.open(file.c_str(), ios::out | ios::binary);
		push(io::gzip_compressor());
	} 
	else if (binary) {
		fd.open(file.c_str(), ios::out | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

void ofile::writeString(string & str) {
	int s = str.size();
	write((char*)&s, sizeof(int));
	write(str.c_str(), s * sizeof(char));
}

void ofile::close() {
	 if (!empty()) reset();
	 fd.close();
}

bool fileexists(string fname){
  ifstream ifile(fname.c_str());
  return ifile;
}
