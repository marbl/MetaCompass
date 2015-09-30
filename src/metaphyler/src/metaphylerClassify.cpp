#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::ios_base;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::istringstream;

#include <map>
using std::map;

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <string>
using std::string;

#include <algorithm>
#include <functional>
using std::greater;

#include <cstdlib>

#include <utility>
using std::pair;

const float PCTCUT = 0;

typedef unsigned short int   Usint;
typedef unsigned int         Uint;
typedef vector<string>       VS;
typedef map<string, VS>      S2VS;
typedef vector<Usint>        VSI;
typedef map<string, VSI>     S2VSI;
typedef map<string, Usint>   S2SI;
typedef map<Usint, S2VSI>    SI2S2VSI;
typedef vector<float>        VF;

Uint LENCUT = 0;

struct Cmdopts{
  VS     scorefiles;
  string taxfile;
  string blastfile;
};

void helpmsg();
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts);
void readTaxFile(string taxfile, S2VS &seq2tax, S2SI &seq2nlevs);
void getScores(string scorefile, SI2S2VSI &len2seq2scores);
void setScores(S2SI &seq2nlevs, S2VSI &seq2scores);
void classifyBLAST(string blastfile, S2SI &seq2nlevs, SI2S2VSI &len2seq2scores, S2VS &seq2tax);
VF   computeConf(string rid, Uint bit, Usint nlevs, S2VSI &seq2scores);
void printSeq2Scores(S2SI &seq2nlevs, S2VSI &seq2scores);
void printClassification(const VF& confs, const VS& tax, const string qid);
inline float average(const VF &ary);


int main(int argc, char *argv[]) {

  // read in command-line options
  Cmdopts cmdopts;
  getcmdopts(argc, argv, cmdopts);


  // read in taxonomic labels for each reference gene
  S2VS  seq2tax;
  S2SI  seq2nlevs;
  readTaxFile(cmdopts.taxfile, seq2tax, seq2nlevs);

  
  // read in cutoff file
  SI2S2VSI len2seq2scores;
  S2VSI seq2scores;
  for (VS::const_iterator citer = cmdopts.scorefiles.begin(); citer != cmdopts.scorefiles.end(); ++citer) {
    getScores(*citer, len2seq2scores);
  }

  // prepare for classification
  for (SI2S2VSI::iterator citer = len2seq2scores.begin(); citer != len2seq2scores.end(); ++citer) {
    setScores(seq2nlevs, citer->second);
    //printSeq2Scores(seq2nlevs, citer->second);
  }

  
  classifyBLAST(cmdopts.blastfile, seq2nlevs, len2seq2scores, seq2tax);
  
  return 0;
}


// the average of an array
inline float average(const VF &ary) {
  float ave = 0;
  for (VF::const_iterator citer = ary.begin(); citer != ary.end(); ++citer)
    ave += *citer;
  return ave/ary.size();
}


// read BLAST file, classify query reads
void classifyBLAST(string blastfile, S2SI &seq2nlevs, SI2S2VSI &len2seq2scores, S2VS &seq2tax) {

  ifstream blastfile_ifs(blastfile.c_str());
  if (!blastfile_ifs) {
    cerr << "Could not open file: " << blastfile << endl;
    exit(1);
  }


  string      eachline;
  set<string> seqids;   // keeps sequences that have been processed
  while (getline(blastfile_ifs, eachline)) {

    // get query read ID
    size_t pos1 = eachline.find("\t");
    string qid  = eachline.substr(0, pos1);
    if (seqids.find(qid) != seqids.end()) continue; // has been processed
    seqids.insert(qid);

    // get reference sequence ID
    size_t pos2 = eachline.find("\t", pos1+1);
    string rid  = eachline.substr(pos1+1, pos2-pos1-1);
    S2SI::const_iterator s2si_citer = seq2nlevs.find(rid);
    if (s2si_citer == seq2nlevs.end()) continue; // does not have classifier for it
    Usint nlevs = s2si_citer->second;
    if (seq2tax.find(rid) == seq2tax.end()) continue;     // does not have taxonomic label

    // get % identity
    size_t pos3 = eachline.find("\t", pos2+1);
    float  pct  = atof(eachline.substr(pos2+1, pos3-pos2-1).c_str());
    if (pct < PCTCUT) continue;

    // HSP length
    size_t pos4 = eachline.find("\t", pos3+1);
    Usint  hsp  = atoi(eachline.substr(pos3+1, pos4-pos3-1).c_str());
    if (hsp < LENCUT) continue;

    // get bit score
    size_t pos5     = eachline.rfind(" ") != string::npos ? eachline.rfind(" ") : eachline.rfind("\t"); // sometimes an extract space before bit score
    Uint bit = atoi(eachline.substr(pos5+1, eachline.size()-pos5-1).c_str());


    VF confs; // confidence scores at each level

    // if hsp length is smaller than the shortest length from available models,
    // then use it, but do not scale bit socre
    if (hsp < len2seq2scores.begin()->first)
      confs = computeConf(rid, bit, nlevs, len2seq2scores.begin()->second);

    // if hsp length is bigger than the longest length from available models,
    // then use it, scale the bit score according to length
    else if (hsp > len2seq2scores.rbegin()->first)
      confs = computeConf(rid, bit*(len2seq2scores.rbegin()->first)/hsp, nlevs, len2seq2scores.rbegin()->second);

    else {

      // iterate through all models for different read lengths
      // suppose hsp length is 150bp; we have models for 100bp and 200 bp
      // then we try classification using both models,
      // and use the one with higher average confidence score
      for (SI2S2VSI::iterator citer = len2seq2scores.begin(); citer != len2seq2scores.end(); ++citer) {

	if (hsp == citer->first) { // if exactly the same, then just use this model
	  confs = computeConf(rid, bit, nlevs, citer->second);
	  break;
	}
	
	else if (hsp < citer->first) {
	  confs = computeConf(rid, bit, nlevs, citer->second);
	  VF confs2 = computeConf(rid, bit*(citer->first)/hsp, nlevs, (--citer)->second);
	  if (average(confs) < average(confs2)) confs = confs2;
	  break;
	}
      }
    }
	if (pct >= 95) {
		confs = VF(nlevs-1, 1.000);
	}
    if (*max_element(confs.begin(), confs.end()) >= 0.001)
      printClassification(confs, seq2tax.find(rid)->second, qid);
	else {     

		printClassification(confs, seq2tax.find(rid)->second, qid);
	}
  }
}


// print out classification information
void printClassification(const VF& confs, const VS& tax, const string qid) {

  cout.setf(ios_base::fixed);
  cout.precision(3);
  cout << qid << "\t";
  for (Usint i = 0; i < confs.size(); ++i) {
    if (tax[i] == "NA") 
      cout << tax[i] << "\t";
    else
      cout << tax[i] << "(" << confs[i] << ")\t";
  }
  cout << endl;
}


// compute confidence scores at each taxonomic level
VF computeConf(string rid, Uint bit, Usint nlevs, S2VSI &seq2scores) {

  VF confs(nlevs-1, 0.0);

  if (seq2scores.find(rid) == seq2scores.end())
    return confs;

  VSI &scores = seq2scores.find(rid)->second;
  if (scores.empty())
    return confs;

  if (bit*nlevs > scores.size()) {
    bit = scores.size() / nlevs;
  }
  
  // if score is smaller than biggest score in model
  // computer conf, otherwise conf is 1
  for (int i = 0; i < nlevs-1; ++i) { // try to classify at each level
      
    Uint samen = 0, diffn = 0;
    for (int j = 0; j < nlevs; ++j) {
      size_t loc = (bit-1)*nlevs+j;
      j <= i ? samen += scores[j] - scores[loc] : diffn += scores[loc];
    }

    if (samen != 0 || diffn != 0)
      confs[i] = samen*1.0 / (samen+diffn);
  }

  return confs;
}


// suppose within a same taxonomic level
// 10 sequences > 100, and next 20 sequences > 90
// then we also set values between 90-100 to be 10
void setScores(S2SI &seq2nlevs, S2VSI &seq2scores) {
  
  for (S2SI::iterator citer = seq2nlevs.begin(); citer != seq2nlevs.end(); ++citer) {
    
    Usint nlevs = citer->second;
    S2VSI::iterator siter = seq2scores.find(citer->first);
    if (siter == seq2scores.end()) { continue;}
    for (int i = 0; i < nlevs; ++i) {
      Uint prenum = 0;
      for (int j = siter->second.size() - 1 - i; j >= 0 ; j -= nlevs)
	siter->second[j] != 0 ? prenum = siter->second[j] : siter->second[j] = prenum;
    }
  }
}


// store classification scores
void getScores(string scorefile, SI2S2VSI &len2seq2scores) {

  ifstream scorefile_ifs(scorefile.c_str());
  if (!scorefile_ifs) {
    cerr << "Could not open file: " << scorefile << endl;
    exit(1);
  }

  string eachline, seqid, eachword;
  istringstream iss;

  S2VSI *seq2scores = &len2seq2scores.begin()->second;

  Usint nlevs = 0, lev = 0;
  S2VSI::iterator siter;
  while(getline(scorefile_ifs, eachline)) {

    iss.clear();
    iss.str(eachline);

    if (eachline[0] == '#') {

      // 1st line
      Usint length = 0;
      iss >> eachword >> length;

      // 2nd line
      string blast("");
      getline(scorefile_ifs, eachline);
      iss.clear();
      iss.str(eachline);
      iss >> eachword >> blast;
      if (blast == "blastx") {
	LENCUT = 20;
	length /= 3;
      }
  
      // model with same length has been observed before
      if (len2seq2scores.find(length) != len2seq2scores.end()) { 
	cerr << "Length " << length << " has been observed before at file " << scorefile << endl;
	exit(1);
      }
      seq2scores = &len2seq2scores.insert(SI2S2VSI::value_type(length, S2VSI())).first->second;

      // 3rd line
      getline(scorefile_ifs, eachline); // third line: if this model uses normalization
      continue;
    }

    if (eachline.empty()) { // empty line, does not have scores from same taxonomic unit
      if (lev == 0)         // if it's at the lowest level, then initialize containers
	siter = seq2scores->insert(S2VSI::value_type(seqid, VSI())).first;
      ++lev;
    }
    
    else if (eachline[0] == '>') { // this line is a header
      iss >> seqid >> nlevs;
      seqid.erase(0, 1);  // extract sequence ID
      lev = 0;            // set current level to 0
      nlevs++;            // +1 to accommondate "other" level
    }

    else { // nonempty line, contains scores
      Uint score, num;
      iss >> score >> num;

      if (lev == 0)      // lowest taxonomic level, initialize scores list
	siter = seq2scores->insert(S2VSI::value_type(seqid, VSI(score*nlevs, 0))).first;

      // if the largest score is bigger than existing one, resize and initialize
      else{
	size_t newsize = score*nlevs;
	size_t oldsize = siter->second.size();
	if (newsize > oldsize) {
	  siter->second.resize(newsize);
	  fill(siter->second.begin()+oldsize, siter->second.end(), 0.0);
	}


      }

      // store all the scores
      siter->second[(score-1)*nlevs+lev] = num;
      while(iss >> score >> num)
	siter->second[(score-1)*nlevs+lev] = num;
    
      ++lev; // increase tax level
    }
  }
}


// parse command line options
void getcmdopts(int argc, char *argv[], Cmdopts &cmdopts) {

  if (argc != 4) {
    helpmsg();
    exit(1);
  }

  cmdopts.taxfile   = argv[2];
  cmdopts.blastfile = argv[3];

  
  // could be multiple files
  string filestr = argv[1]; 
  size_t prepos  = 0;
  for (Usint i = 0; i < filestr.size(); ++i) {
    if (filestr[i] == ',') {
      cmdopts.scorefiles.push_back(filestr.substr(prepos, i-prepos));
      prepos = i + 1;
    }
  }

  // if input files string is "A,B,", then no need to store last record.
  if (filestr[filestr.size()-1] != ',') 
    cmdopts.scorefiles.push_back(filestr.substr(prepos, filestr.size()-prepos));
  
}


// print out usage help message
void helpmsg() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./metaphylerClassify <classifiers> <taxonomy file> <BLAST file>" << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "      <classifiers>   Output from program blast2TaxScores." << endl;
  cerr << "                      If there are multiple files, separate them with comma(e.g., fileA,fileB)" << endl << endl;

  cerr << "      <taxonomy file> Taxonomy labels of reference sequences in the BLAST file." << endl << endl;

  cerr << "      <BLAST file>    BLAST alignment between query reads and reference sequences." << endl << endl;

  cerr << endl;

}


// read in taxonomic profile of reference sequences
void readTaxFile(string taxfile, S2VS &seq2tax, S2SI &seq2nlevs) {

  ifstream ifs(taxfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << taxfile << endl;
    exit(1);
  }

  istringstream iss;
  string eachline;
  while (getline(ifs, eachline)) {

    iss.clear();
    iss.str(eachline);

    string seqid;
    iss >> seqid;          // 1st column is sequence ID
    
    VS tlabs;
    string eachword;
    while(iss >> eachword) // other columns are taxonomic IDs 
      tlabs.push_back(eachword);
    seq2tax.insert(S2VS::value_type(seqid, tlabs));
    seq2nlevs.insert(S2SI::value_type(seqid, tlabs.size()+1));
  }
}


// for debuging purposes
void printSeq2Scores(S2SI &seq2nlevs, S2VSI &seq2scores) {
  
  for (S2SI::const_iterator citer = seq2nlevs.begin(); citer != seq2nlevs.end(); ++citer) {
    Usint nlevs = citer->second;

    S2VSI::const_iterator siter = seq2scores.find(citer->first);
    if (siter == seq2scores.end()) { continue;}
    cout << citer->first << "\t" << siter->second.size() << endl;

    for (Uint i = 0; i < nlevs; ++i) {
      for (Uint j = i; j < siter->second.size(); j += nlevs) {
	if (siter->second[j] > 0)
	  cout << j / nlevs + 1 << "\t" << siter->second[j] << "\t";
      }
      cout << endl;
    }

  }
}

