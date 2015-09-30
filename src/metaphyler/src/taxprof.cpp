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

#include <cstdlib>

typedef unsigned int Uint;
const Uint TLEV = 6;

struct Cmdopts {
  string clsffn,
         prefix,
         tnamesfn;
  float  confcut;
};

typedef map<string, string> S2S;
typedef map<string, Uint>   S2I;
typedef vector<S2I>         VS2I;

void helpmsg();
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts);
void gettnames(string tnamesfn, S2S &tid2name);
Uint abundance(const Cmdopts &cmdopts, const S2S &tid2name, VS2I &abund);
void printtaxprof(const VS2I &taxprof, Uint n, string prefix);

int main(int argc, char *argv[]) {

  // read in command line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);


  S2S tid2name;
  if (cmdopts.tnamesfn != "")
    gettnames(cmdopts.tnamesfn, tid2name);


  VS2I abund(TLEV, S2I());
  Uint totaln = abundance(cmdopts, tid2name, abund);

  printtaxprof(abund, totaln, cmdopts.prefix);
  
  return 0;
}

void printtaxprof(const VS2I &taxprof, Uint n, string prefix) {

  vector<string> levnames;
  levnames.push_back("species");
  levnames.push_back("genus");
  levnames.push_back("family");
  levnames.push_back("order");
  levnames.push_back("class");
  levnames.push_back("phylum");


  Uint i = 0;
  for (VS2I::const_iterator citer1 = taxprof.begin(); citer1 != taxprof.end(); ++citer1, ++i) {
    if (citer1->empty()) continue;

    string outfile = prefix + "." + levnames[i] + ".taxprof";
    ofstream ofs(outfile.c_str());
    ofs.setf(ios_base::fixed);
    ofs.precision(2);
    Uint sum = 0;
    ofs << "Name\t% Abundance\t# reads" << endl;
    for (S2I::const_iterator citer2 = citer1->begin(); citer2 != citer1->end(); ++citer2) {
      ofs << citer2->first << "\t" << citer2->second*100.0/n << "\t" << citer2->second << endl;
      sum += citer2->second;
    }
    if (sum < n) {
      ofs << "Other\t" << (n-sum)*100.0/n << "\t" << n-sum << endl;
    }
  }
}

Uint abundance(const Cmdopts &cmdopts, const S2S &tid2name, VS2I &abund) {

  Uint n = 0;
  ifstream ifs(cmdopts.clsffn.c_str());
  if (!ifs) {
    cerr << "Could not open file " << cmdopts.clsffn << endl;
    exit(1);
  }

  string eachline, eachword, tid;
  istringstream iss;
  while (getline(ifs, eachline)) {

    iss.clear();
    iss.str(eachline);
    iss >> eachword;

    Uint lev = 0;
    bool tag = 0;
    while (iss >> eachword) {
      ++lev;
      if (eachword == "NA") continue;

      size_t pos = eachword.find('(');
      float conf = atof(eachword.substr(pos+1, 5).c_str());
      if (conf < cmdopts.confcut) continue;
      string tname = eachword.substr(0, pos);
      S2S::const_iterator citer = tid2name.find(tname);
      if (citer != tid2name.end())
	tname = citer->second;

      S2I::iterator iter = abund[lev-1].find(tname);
      if (iter == abund[lev-1].end())
	abund[lev-1].insert(S2I::value_type(tname, 1));
      else
	++(iter->second);
      
      tag = true;
      
    }
    if (tag) ++n;
    
  }
  return n;
}

void gettnames(string tnamesfn, S2S &tid2name) {

  ifstream ifs(tnamesfn.c_str());
  if (!ifs) {
    cerr << "Could not open file " << tnamesfn << endl;
    exit(1);
  }

  string eachline, tid, tname;
  istringstream iss;
  while (getline(ifs, eachline)) {
    iss.clear();
    iss.str(eachline);
    iss >> tid >> tname;
    tid2name.insert(S2S::value_type(tid, tname));
  }
}

// parse command line options
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 5 && argc != 4) {
    helpmsg();
    exit(1);
  }

  cmdopts.confcut = atof(argv[1]);
  cmdopts.clsffn  = argv[2];
  cmdopts.prefix  = argv[3];

  argc == 5 ? cmdopts.tnamesfn = argv[4] : cmdopts.tnamesfn = "";
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./taxprof <conf. cutoff> <classification> <prefix> <taxonomy names>" << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "        <conf. cutoff>   Cutoff for confidence score. Higher means successfully classified." << endl;
  cerr << "                         Recommendation: 0.9." << endl;;
  cerr << "        <classification> Result file from program metaphylerClassify." << endl;
  cerr << "        <prefix>         Output files prefix." << endl;
  cerr << "        <taxonomy names> File: 1st column, taxonomy ID; 2nd, name." << endl;
  cerr << "                         If omitted, output will just use taxonomy IDs." << endl;;

  cerr << "Output files:" << endl;
  cerr << "        prefix.<genus|family|order|class|phylum>.taxprof." << endl;
  cerr << "                         Taxonomy profiles at each level." << endl << endl;;
  
  cerr << endl;
}
