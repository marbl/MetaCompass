#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::ios_base;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::istringstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <set>
using std::set;

#include <cstdlib>


typedef set<string> SS;

struct Cmdopts {
  string fafile;
  string outprefix;
};

void helpmsg();
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts);
void extractseq(const string fafile, const string outprefix);
void reverseStr(string & str);

int main(int argc, char *argv[]) {

  // read command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);
  
  // extract and output genome sequences
  extractseq(cmdopts.fafile, cmdopts.outprefix);
  
  
  return 0;
}


void reverseStr(string& str) 
{ 
    int n = str.length(); 
  
    // Swap character starting from two 
    // corners 
    for (int i = 0; i < n / 2; i++) 
        std::swap(str[i], str[n - i - 1]); 
} 
// extract and output genome sequences
void extractseq(const string fafile, const string outprefix) {
  ifstream ifs(fafile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << fafile << endl;
    exit(1);
  }
  string eachline,seqid;
  istringstream iss;
  ofstream ctgfile;
  while (getline(ifs, eachline)) {
    if (eachline[0] == '>') {
		iss.clear();
		iss.str(eachline);
      	iss >> seqid;
      	seqid.erase(0, 1);
	  	reverseStr(seqid);
	  	// Find position of '_' using find() 
      	int pos = seqid.find("_"); 
		//only for pilon reads
		int pos2 = seqid.find("_", pos+1);
		
	  	// Copy substring after pos
	  	seqid = seqid.substr(pos2 + 1);
	  	//cout << "seqid:" <<seqid <<endl;
		
	  	reverseStr(seqid);    
	  	ctgfile.open(outprefix + "/" + seqid+".fasta", std::ios_base::app);
	  	ctgfile << eachline <<endl;
	  	ctgfile.close();
  	}
    else {
		ctgfile.open(outprefix + "/" + seqid+".fasta", std::ios_base::app);
		ctgfile << eachline <<endl;
		ctgfile.close();
	}
  }
}

/*
void extractseqold2(const string fafile, const SS &ids, const string outprefix) {
  ifstream ifs(fafile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << fafile << endl;
    exit(1);
  }
  string eachline, eachword, tid, seqid,tmp="";
  bool tag = 0;
  istringstream iss;
  ofstream ctgfile;
  while (getline(ifs, eachline)) {
    if (eachline[0] == '>') {
      iss.clear();
      iss.str(eachline);
      iss >> seqid;
      seqid.erase(0, 1);
	  reverseStr(seqid);
	  // Find position of '_' using find() 
      int pos = seqid.find("_"); 
	  // Copy substring after pos
	  seqid = seqid.substr(pos + 1);
	  reverseStr(seqid);    
      if (ids.find(seqid) != ids.end()) {
		  ctgfile.open(outprefix + "/" + seqid+".fasta", std::ios_base::app);
		  ctgfile << eachline <<endl;
		  ctgfile.close();
		  tag = 1;
	  }
      else
		  tag = 0;
	  tmp=seqid;
    }
    else 
		if (tag){//cout << eachline << endl;
			ctgfile.open(outprefix + "/" + tmp+".fasta", std::ios_base::app);
			ctgfile << eachline <<endl;
			ctgfile.close();
		}
	}
}

void extractseqold(const string fafile, const SS &ids, const string outprefix) {

  ifstream ifs(fafile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << fafile << endl;
    exit(1);
  }
  string eachline, eachword, tid, seqid,tmp="";
  bool tag = 0;
  istringstream iss;
  ofstream ctgfile;
  while (getline(ifs, eachline)) {
    if (eachline[0] == '>') {
      iss.clear();
      iss.str(eachline);
      iss >> seqid;
      seqid.erase(0, 1);
	  reverseStr(seqid);
	  // Find position of '_' using find() 
      int pos = seqid.find("_"); 
	  // Copy substring after pos
	  seqid = seqid.substr(pos + 1);
	  reverseStr(seqid);    
      if (ids.find(seqid) != ids.end()) {
		  //cout << eachline << endl;
		  if (tmp.compare(seqid)==0)
		  		if (ctgfile.is_open()){
					ctgfile.close();
					ctgfile << "closing file" <<endl;
				}
		  //ctgfile.open(outprefix + "/" + tmp+".fasta", std::ios_base::app);
		  ctgfile.open(outprefix + "/" + seqid+".fasta", std::ios_base::app);
		  if (ctgfile.is_open()){
			  ctgfile << "opening file" <<endl; 
			  ctgfile << eachline <<endl;
			}
  		  else cout << "Unable to open file "+ outprefix + "/" + tmp+ ".fasta"<< endl;
		  tag = 1;
	  }
      else
		  tag = 0;
	  tmp=seqid;
    }
    else 
		if (tag){//cout << eachline << endl;
			ctgfile.open(outprefix + "/" + tmp+".fasta", std::ios_base::app);
			ctgfile << eachline <<endl;
			ctgfile.close();
		}
  }
}

*/
// parse command line options
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 3) {
    helpmsg();
    exit(1);
  }

  cmdopts.fafile = argv[1];
  cmdopts.outprefix = argv[2];
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./extractContigs <fasta file> <output directory>" << endl;
  cerr << endl;
  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;
  cerr << endl;
}
