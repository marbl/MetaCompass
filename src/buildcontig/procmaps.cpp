#include <iostream>
using std::endl;
using std::cout;
using std::cerr;
using std::ios;

#include <fstream>
using std::ifstream;
using std::ofstream;
#include <sstream>
using std::istringstream;

#include <ctime>
#include <climits>
#include <cstdlib>
#include <algorithm>

#include "procmaps.hpp"
#include "memory.hpp"


// convert nucleotide to char number
/***************************************/
VC b2n(UCHAR_MAX, 4);
static void init(VC &b2n) {
  b2n['A'] = b2n['a'] = 0;
  b2n['C'] = b2n['c'] = 1;
  b2n['G'] = b2n['g'] = 2;
  b2n['T'] = b2n['t'] = 3;
  b2n['-'] = 4;
}
/***************************************/

//Read reference to file
void readrefseqfile(const string refseqfile, S2S &ref2seq) {


  ifstream refseqfiles(refseqfile.c_str());
  if (!refseqfiles) {
    cerr << "Could not open reference sequences file " << refseqfile << endl;
    exit(1);
  }


  string line(""), refid(""),refseq("");   
  refseq.reserve(10000000);
  Uint len = 0;
  while (getline(refseqfiles, line)) {

    // fasta header line
    if (line[0] == '>') {             
      
      // not the first line of file, store previous fasta record
      if (!refid.empty()) {
	refseq.resize(refseq.size());
	ref2seq.insert(S2S::value_type(refid, refseq) );
      }
      
      // reset variables
      refseq.clear();
      refid.clear();
      len = 0;
      
      // parse fasta header line to get fasta id
      istringstream iss(line);
      iss >> refid;       // only read the first word
      refid.erase(0, 1);  // remove first char '>'
    }
    
    // sequence line
    else { 
      refseq += line;       // append sequence
      len += line.size(); // increase length
    }
    
  }
  /****************************************************************************/
  
  // store the last fasta record
  refseq.resize( refseq.size() );
  ref2seq.insert(S2S::value_type(refid, refseq) );
}
//deleted comments to stdout, march 15 2016
//mutstr = eachline.substr(pos1+1, pos2-pos1-1);//CIGAR string representation of alignment
//this is a char, only analyze M: match, I: insertion, D: deletion
//VC &prof, 
//VB &ins, 
//VB &cov, 
//I2VC &pos2ins
//string mutstr
//string &seq, 
//Uint start
static void storesam(VC &prof, VB &ins, VB &cov, I2VC &pos2ins, string mutstr, string &seq, Uint start) {

  Uint n    = 0;      // tracks length of each region
  Uint locr = start;  // locus in ref
  Uint locq = 0;      // locus in query
  //ofstream ofs("newmap.sam");
      //ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << ctgend << endl;
      //ofs << contig << endl;
  for (size_t i = 0; i < mutstr.size(); ++i) {

    char achar = mutstr[i];

    if (isdigit(achar))        // this is a digit
      n = n*10 + achar - '0';

    // this is a char, only analyze M: match, I: insertion, D: deletion
    else {

  

      if (achar == 'M') {
	/////////////cout << "Match" <<start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	for (Uint i = 0; i < n; ++i) {
	  cov[locr] = true;
	  if (prof[locr*5+b2n[seq[locq]]] != UCHAR_MAX) ++prof[locr*5+b2n[seq[locq]]];
	  ++locr, ++locq;
	}
	--locr;
	//--locq;
      }

      // insertion in query sequence
      else if (achar == 'I') {

	////////////cout << "Inser" << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	if (ins[locr]) { // > 1 times see this insertion

	  /////////cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;

	  if (pos2ins.find(locr) == pos2ins.end())
	    pos2ins.insert(I2VC::value_type(locr, VC(12, 0)));
	  
	  VC &insprof = pos2ins[locr];
	  
	  for (Uint i = 0; i < n; ++i) { // store max 3 insertions
	    /////////cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	    if (i < 3) {
	  //////////    cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << i*4+b2n[seq[locq]] << "\t" << mutstr << endl;
	      if (insprof[i*4+b2n[seq[locq]]] != UCHAR_MAX) ++insprof[i*4+b2n[seq[locq]]];
	     //////////// cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << int(insprof[i*4+b2n[seq[locq]]]) << "\t" << mutstr << endl;
	    }
	    ++locq;
	  }
	}
	else { // first time see this insertion
	  ins[locr] = true;
	  locq += n;
	}
	locr++;
      }


      // deletion in query sequence
      else if (achar == 'D') {
	locr++;
	/////////////////cout << "Delet" << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	for (Uint i = 0; i < n; ++i) {
	  cov[locr] = true;
	  prof[locr*5+4] != UCHAR_MAX ? ++prof[locr*5+4] : 1;
	  ++locr;
	}
	//--locr;
      }
      else {
	locr += n + 1;
	//locq += n - 1;
      }

      n = 0;
    }
  }
}

static void compute_breadth_sam(ifstream &ifs, S2D &ref2bre, S2S &ref2seq) {

  S2VB ref2cov;
  init(ref2seq, ref2cov, 1);

  string eachline, refid;
  size_t pos1, pos2, len;
  Uint start, end;
  while (getline(ifs, eachline)) {
    if ( eachline[0] == '@') continue;
    pos1 = eachline.find("\t");           // 1
    pos1 = eachline.find("\t", pos1+1);   // 2
    if (eachline[pos1-1] == '4') continue;// unmapped
      
    pos2 = eachline.find("\t", pos1+1);   // 3
    refid = eachline.substr(pos1+1, pos2-pos1-1);

    pos1 = eachline.find("\t", pos2+1);   // 4
    start = atoi(eachline.substr(pos2+1, pos1-pos2-1).c_str())-1;
    
    pos1 = eachline.find("\t", pos1+1);   // 5
    pos1 = eachline.find("\t", pos1+1);   // 6
    pos1 = eachline.find("\t", pos1+1);   // 7
    pos1 = eachline.find("\t", pos1+1);   // 8
    pos1 = eachline.find("\t", pos1+1);   // 9
    pos2 = eachline.find("\t", pos1+1);   // 10
    len = eachline.substr(pos1+1, pos2-pos1-1).size();

    end = start+len > ref2cov[refid].size() ? ref2cov[refid].size() : start+len;
    fill(ref2cov[refid].begin()+start, ref2cov[refid].begin()+end, true);
  }
  
  for (S2VB::iterator iter = ref2cov.begin(); iter != ref2cov.end(); ++iter) {
    double cov = 0;
    for (size_t i = 0; i < iter->second.size(); ++i)
      if (iter->second[i])
	++cov;
    ref2bre.insert(S2D::value_type(iter->first, cov));
  }
  
}

static void compute_depth_sam(ifstream &ifs, S2D &ref2dep, S2S &ref2seq) {

  string eachline, refid;
  size_t pos1, pos2, len;
  while (getline(ifs, eachline)) {
      
    if ( eachline[0] == '@') continue;

    pos1 = eachline.find("\t");     // 1
    pos1 = eachline.find("\t", pos1+1); // 2
    if (eachline[pos1-1] == '4') continue; // unmapped
      
    pos2 = eachline.find("\t", pos1+1);   // 3
    refid = eachline.substr(pos1+1, pos2-pos1-1);

    pos1 = eachline.find("\t", pos2+1);   // 4
    pos1 = eachline.find("\t", pos1+1);   // 5
    pos1 = eachline.find("\t", pos1+1);   // 6
    pos1 = eachline.find("\t", pos1+1);   // 7
    pos1 = eachline.find("\t", pos1+1);   // 8
    pos1 = eachline.find("\t", pos1+1);   // 9
    pos2 = eachline.find("\t", pos1+1);   // 10
    len = eachline.substr(pos1+1, pos2-pos1-1).size();

    if (ref2dep.find(refid) == ref2dep.end())
      ref2dep.insert(S2D::value_type(refid, len));
    else 
      ref2dep.find(refid)->second += len; // calculate the sum of reads mapped
  }

  for (S2D::iterator iter = ref2dep.begin(); iter != ref2dep.end(); ++iter) {

    // could not find ref, set to 0
    if (ref2seq.find(iter->first) == ref2seq.end()) 
      iter->second = 0;
    else
      iter->second /= ref2seq.find(iter->first)->second.size(); // normalize by length
  }
  
}

void processsams(const Cmdopt &cmdopt, S2S &ref2seq, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins) {

  init(b2n);
  ofstream ofs(string(cmdopt.outprefix + "/selected_maps.sam").c_str());
  if (!ofs) {
    cerr << "Could not open sam file " << string(cmdopt.outprefix + "/selected_maps.sam").c_str() << endl;
    exit(1);
  }
     // ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << ctgend << endl;
     //ofs << contig << endl;
  
  ifstream ifs(cmdopt.mapfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << cmdopt.mapfile << endl;
    exit(1);
  }

  // use all map records
  if (cmdopt.pickref == "all") {

    // stores each map record
    size_t  pos1 = 0, pos2 = 0;
    string eachline, refid, mutstr, seq;
    Uint start = 0;
    while (getline(ifs, eachline)) {
        if (eachline[0] == '@') {
            ofs << eachline << endl;
            continue;
        }

      pos1 = eachline.find("\t");           // 1
      pos1 = eachline.find("\t", pos1+1);   // 2
      if (eachline[pos1-1] == '4') continue;// unmapped
      
      pos2 = eachline.find("\t", pos1+1);   // 3
      refid = eachline.substr(pos1+1, pos2-pos1-1);

      pos1 = eachline.find("\t", pos2+1);   // 4
      start = atoi(eachline.substr(pos2+1, pos1-pos2-1).c_str())-1;

      pos1 = eachline.find("\t", pos1+1);   // 5
      pos2 = eachline.find("\t", pos1+1);   // 6

      mutstr = eachline.substr(pos1+1, pos2-pos1-1);
      
      pos1 = eachline.find("\t", pos2+1);   // 7
      pos1 = eachline.find("\t", pos1+1);   // 8
      pos1 = eachline.find("\t", pos1+1);   // 9
      pos2 = eachline.find("\t", pos1+1);   // 10
      seq = eachline.substr(pos1+1, pos2-pos1-1);
      //////////////cout << "Ref "<< refid << "\t" <<endl;
      storesam(ref2prof[refid], ref2ins[refid], ref2cov[refid], ref2pos2ins[refid], mutstr, seq, start);
      ofs << eachline << endl;
//ofs << eachline;
    }
  }

  else {

    S2D ref2val;

    if (cmdopt.pickref == "depth")
      compute_depth_sam(ifs, ref2val, ref2seq);

    else if (cmdopt.pickref == "breadth") 
      compute_breadth_sam(ifs, ref2val, ref2seq);
      

    // rewind input file stream
    ifs.clear();
    ifs.seekg(0, ios::beg);

    size_t  pos1 = 0, pos2 = 0;
    std::string eachline(""), refid(""), mutstr(""), seq(""), preqid(""), qid(""), prerefid("");
    Uint start = 0;
    double maxval = 0.0;

    while (getline(ifs, eachline)) {

        if (eachline[0] == '@') {
            //cout << eachline << endl;
            ofs << eachline << endl;
            continue;
        }

      pos1 = eachline.find("\t");           // 1
      qid = eachline.substr(0, pos1);

      pos1 = eachline.find("\t", pos1+1);   // 2
      if (eachline[pos1-1] == '4') continue;// unmapped

      
      pos2 = eachline.find("\t", pos1+1);   // 3
      refid = eachline.substr(pos1+1, pos2-pos1-1);

      S2D::iterator iter = ref2val.find(refid);

      if (qid != preqid && !mutstr.empty()) // store previous record
       {///////////////cout << "Ref "<< prerefid << "\t" << iter->second << endl;
           storesam(ref2prof[prerefid], ref2ins[prerefid], ref2cov[prerefid], ref2pos2ins[prerefid], mutstr, seq, start);
           ofs << eachline << endl;
       }
      else if (iter->second <= maxval)      // not as good as previous one
        continue;

      maxval = iter->second;
      preqid = qid;
      prerefid = refid;
      
      pos1 = eachline.find("\t", pos2+1);   // 4
      start = atoi(eachline.substr(pos2+1, pos1-pos2-1).c_str())-1;

      pos1 = eachline.find("\t", pos1+1);   // 5
      pos2 = eachline.find("\t", pos1+1);   // 6

      mutstr = eachline.substr(pos1+1, pos2-pos1-1);
      
      pos1 = eachline.find("\t", pos2+1);   // 7
      pos1 = eachline.find("\t", pos1+1);   // 8
      pos1 = eachline.find("\t", pos1+1);   // 9
      pos2 = eachline.find("\t", pos1+1);   // 10
      seq = eachline.substr(pos1+1, pos2-pos1-1);
    }

//////////    cout << "Ref "<< refid << "\t" <<endl;
    storesam(ref2prof[refid], ref2ins[refid], ref2cov[refid], ref2pos2ins[refid], mutstr, seq, start);
    ofs << eachline;
  }
}
