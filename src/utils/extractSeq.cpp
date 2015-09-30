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
  string idfile;
};

void helpmsg();
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts);
void getids(const string idfile, SS &ids);
void extractseq(const string fafile, const SS &ids);

int main(int argc, char *argv[]) {

  // read command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);

  // reference genome ids
  SS ids;
  getids(cmdopts.idfile, ids);
  if (ids.empty()) return 0;
  
  // extract and output genome sequences
  extractseq(cmdopts.fafile, ids);
  
  return 0;
}


// extract and output genome sequences
void extractseq(const string fafile, const SS &ids) {

  ifstream ifs(fafile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << fafile << endl;
    exit(1);
  }

  string eachline, eachword, tid, seqid;
  bool tag = 0;
  istringstream iss;
  while (getline(ifs, eachline)) {
    if (eachline[0] == '>') {
      iss.clear();
      iss.str(eachline);
      iss >> seqid;
      seqid.erase(0, 1);

      if (ids.find(seqid) != ids.end()) {
	cout << eachline << endl;
	tag = 1;
      }
      else
	tag = 0;
      
    }
    else if (tag)
      cout << eachline << endl;
  }
}


// read reference genome ids
void getids(const string idfile, SS &ids) {

  ifstream ifs(idfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << idfile << endl;
    exit(1);
  }

  string eachline, seqid;
  istringstream iss;
  while (getline(ifs, eachline)) {
    iss.clear();
    iss.str(eachline);
    iss >> seqid;
    ids.insert(seqid);
  }
}


// parse command line options
void getcmdopts(const int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 3) {
    helpmsg();
    exit(1);
  }

  cmdopts.fafile = argv[1];
  cmdopts.idfile = argv[2];
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./extractSeq <fasta file> <IDs file>" << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "        <IDs file> 1st column is sequence IDs." << endl;
  
  cerr << "Contact:" << endl;
  cerr << "        Have problems? Contact Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;
  cerr << endl;
}
