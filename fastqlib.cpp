#include "fastqlib.h"


fqread::fqread(string header,string dna,string line3,string qual)
{
    set( header, dna, line3, qual)  ;
}
  
int fqread::set(string header,string dna,string line3,string qual)
{
    description=true;
    h=header;
    s=dna;
    l3=line3;
    q=qual;
    l=s.size();
    istringstream iss;
    iss.str(header);
    string tmp;
    iss >> tmp;
    iss >> tmp;
    if(!iss)
    {
	filtered=false;
	description=false;
    }
    else
    {
	if(tmp.substr(1,3)==":Y:")
	    filtered=true;
	else
	    filtered=false;
    }
    if(dna.size()!=qual.size())
    {
	cerr << header << endl;
	cerr << dna    << endl;
	cerr << line3  << endl;
	cerr << qual << endl;
	die("base length was not equal to qual length in read");
    }
    return(0);
}

fqread::fqread()
{
    l=0;
}

int fqread::clear()
{
    l=0;
    h="";
    s.clear();
    q.clear();
    l3="";  
    return(0);
}

fqread::fqread(int L)
{
    l=L;
    h="";
    s.resize(L);
    q.resize(L);
    l3="";
}

readPair::readPair()
{
}

readPair::readPair(fqread read1,fqread read2)
{
    r1=read1;
    r2=read2;
    filtered = read1.filtered || read2.filtered;
}

int readPair::rc()
{
    r1 = r1.rc();
    r2 = r2.rc();
    return(0);
}

fqread fqread::rc()
{
    fqread ret(l);
    ret.h=h;
    ret.l3=l3;
    for(int i=0;i<l;i++)
    {
	ret.q[l-i-1]=q[i];
	char base='N';
	switch(s[i])
	{
	case 'C': 
	    base='G';
	    break;
	case 'G': 
	    base='C';
	    break;
	case 'T': 
	    base='A';
	    break;
	case 'A': 
	    base='T';         
	    break;
	}
	ret.s[l-i-1]=base;
    }  
    return(ret);
}

int fqread::notN()
{
    return(notN(0,l));
}

int fqread::notN(int a,int b)
{
    assert(b>a);
    int count = 0;
    for(int i=a;i<b;i++)
	if(s[i]=='N' || ((unsigned int)q[i])<((unsigned int)'%'))
	    count++;
    return(b - a - count);
}

fqread fqread::mask()
{ //N masks the entire read
    return(mask(0,l));
}

fqread fqread::mask(int a,int b)
{ //N masks the region a,b
    assert(a>=0&&b<=l);
    string new_s = s;
    string new_q = q;
    for(int i=a;i<b;i++)
    {
	new_s[i]='N';
	new_q[i]='#';
    }
    return(fqread(h,new_s,l3,new_q));
}

fqread fqread::window(int a,int b)
{
    return(fqread(h,s.substr(a,b-a),l3,q.substr(a,b-a)));
}

fqread fqread::window(int a)
{
    return(fqread(h,s.substr(a),l3,q.substr(a)));
}

void fqread::print()
{
    if(l>0)
    {
	cout << "@"<<h << endl;
	cout << s << endl;
	cout <<l3 << endl;
	cout << q << endl;
    }
}


fastqReader::fastqReader(string fname)
{
    warned=false;
    fp = gzopen(fname.c_str(), "r");
    seq = kseq_init(fp);
    if(!fp )
    {
	cerr << "Problem reading "<<fname<<endl;
	exit(1);
    }
}
fastqReader::~fastqReader()
{
    kseq_destroy(seq);
    gzclose(fp);
}

fastqWriter::fastqWriter()
{
    fp=NULL;
}

fastqWriter::~fastqWriter()
{
    if(fp)
	assert(Z_OK==gzclose(fp));
}

fastqWriter::fastqWriter(string fname)
{
    open(fname);
}

int fastqWriter::open(string fname)
{
    if(fname=="-")
    {//STDOUT
	_stdout=true;
    }
    else
    {
	_stdout=false;
	fp = gzopen(fname.c_str(), "wb");
	if(!fp )
	{
	    cerr << "Problem reading "<<fname<<endl;
	    exit(1);
	}
    }
    return(0);
}

int fastqWriter::write(fqread & read)
{
    if(read.l>0)
    {
	if(_stdout)
	{
	    read.print();
	}
	else
	{
	    bool write_ok=true;
	    write_ok = write_ok && gzwrite(fp,"@",1)!=0;
	    write_ok = write_ok && gzwrite(fp,(char *)read.h.c_str(),read.h.size())!=0;
	    write_ok = write_ok && gzwrite(fp,"\n",1)!=0;
	    write_ok = write_ok && gzwrite(fp,(char *)read.s.c_str(),read.s.size())!=0;
	    write_ok = write_ok && gzwrite(fp,"\n",1)!=0;
	    write_ok = write_ok && gzwrite(fp,(char *)read.l3.c_str(),read.l3.size())!=0;
	    write_ok = write_ok && gzwrite(fp,"\n",1)!=0;
	    write_ok = write_ok && gzwrite(fp,(char *)read.q.c_str(),read.q.size())!=0;
	    write_ok = write_ok && gzwrite(fp,"\n",1)!=0;
	    if(!write_ok)
		die("problem writing output");
	}
	return(1);
    }
    else  return(0);
}

int fastqWriter::write(readPair & p)
{
    if(p.r1.l>0 && p.r2.l>0)
    {
	write(p.r1);
	write(p.r2);
	return(1);
    }
    else
	return(0);
}

int fastqReader::next(fqread & r)
{
    if(kseq_read(seq)<0)
	return(0);

    if(seq->comment.s==NULL)
	r.set((string)seq->name.s,(string)seq->seq.s,"+",(string)seq->qual.s);
    else
	r.set((string)seq->name.s+" "+(string)seq->comment.s,(string)seq->seq.s,"+",(string)seq->qual.s);
  
    if(!warned&&!r.description&&fp)
    {
	cerr << "WARNING: no description found in read header.  Assuming read passed passed chastity/purity filters." << endl;
	warned=true;
    }    
    return(1);
}

pairReader::pairReader(string fname1,string fname2)
{
    f1 = new fastqReader(fname1);
    f2 = new fastqReader(fname2);  
}


void readPair::print()
{
    r1.print();
    r2.print();
}


int pairReader::next(readPair & p)
{
    bool ret = f1->next(p.r1)&&f2->next(p.r2);
    if(ret) 
	p.filtered = p.r1.filtered || p.r2.filtered;
    return(ret);
}

pairWriter::pairWriter()
{
}

pairWriter::pairWriter(string fname)
{
    open(fname);
}

int pairWriter::open(string fname)
{
    outfile.open(fname);
    separate=false;
    return(0);
}

pairWriter::pairWriter(string fname1,string fname2)
{
    open(fname1,fname2);
}

int pairWriter::open(string fname1,string fname2)
{
    outfile1.open(fname1);
    outfile2.open(fname2);
    separate=true;
    return(0);
}

int pairWriter::write(readPair & p)
{
    if(p.r1.l>0 && p.r2.l>0)
    {
	if(separate)
	{
	    outfile1.write(p.r1);
	    outfile2.write(p.r2);
	}
	else
	{
	    outfile.write(p);
	}
	return(1);
    }
    else
	return(0);
}
