#include <iostream>
using std::endl;
using std::cout;
using std::cerr;
using std::ios;

#include <fstream>
using std::ifstream;

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


// update prof and cov, for region [start, end) of seq
static void updatebaseprof(const Uint start, const Uint end, const string &seq, VC &prof, VB &cov) {
  for (Uint i = start; i < end && i < seq.size(); ++i) {
    cov[i] = true;
    Uint j = i*5+b2n[seq[i]];
    if (prof[j] != UCHAR_MAX) ++prof[j];
  }
}


// parse and store mummer-map record
static void storemap(const Bestmap &bestmap, S2VB &ref2ins, S2VB &ref2cov, S2S &ref2seq, S2I2VC &ref2pos2ins, S2VC &ref2prof) {

  string refid  = bestmap.refid,   // ref seq id
         misstr = bestmap.misstr;  // alignment string
  Uint   start  = bestmap.start;   // start in ref seq
  VB     &ins   = ref2ins[refid];  // if insert after each position
  VB     &cov   = ref2cov[refid];  // if covered of each pos
  VC     &prof  = ref2prof[refid]; // baseprof of each pos
  
  if (refid.empty()) return;
  
  // perfect match
  if (misstr.empty()) {
    updatebaseprof(start, bestmap.len + start, ref2seq[refid], prof, ref2cov[refid]);
    return;
  }

  
  Uint gapnum    = 0;      // # consecutive gaps
  Uint qgapnum   = 0;      // total # gaps in query seq
  Uint rgapnum   = 0;      // total # gaps in ref seq
  Uint premutloc = start;  // previous mutation locus in ref
  Uint preloc    = 0;      // for parsing each mutation in alignment string
  size_t loc       = misstr.find(';', preloc); // separator of mutation string

  // process mutation string
  while (loc != string::npos) {

    Uint mutloc    = atoi(misstr.substr(preloc, loc - preloc - 2).c_str()) - 1; // mutation locus in query
    Uint mutlocref = start+mutloc+qgapnum-rgapnum;                              // mutation locus in ref
    char qbase     = misstr[loc-2]; // base in query
    char rbase     = misstr[loc-1]; // base in ref

    if ('-' == qbase) ++mutlocref;
    cov[mutlocref] = true;

    //cout << premutloc << "\t" << mutlocref << "\t" << misstr << endl;
      
    // update previous exact match region, between [premutloc, mutlocref)
    updatebaseprof(premutloc, mutlocref, ref2seq[refid], prof, ref2cov[refid]);


    // gap in reference
    if ('-' == rbase) {

      //cout << premutloc << "\t" << mutlocref-1 << "\t" << misstr << endl;
	
      // 1st time see gap at this locs
      // assume it is a random error
      S2I2VC::iterator iter = ref2pos2ins.find(refid);
      if (ins[mutlocref-1] == false) {
	//cout << "true" << endl;
	ins[mutlocref-1] = true;
      }
	
      // >1 times see gap at this locus
      else {

	I2VC::iterator insiter = iter->second.find(mutlocref-1);

	if (gapnum == 0) // start of a new ins
	  if (insiter == iter->second.end())
	    insiter = iter->second.insert(I2VC::value_type(mutlocref-1, VC(12, 0))).first;

	if (insiter != iter->second.end()) {
	  VC &insprof = insiter->second;
	  if (gapnum < 3) {
	    if (insprof[gapnum*4+b2n[qbase]] != UCHAR_MAX) ++insprof[gapnum*4+b2n[qbase]];
	    //cout << int(insprof[gapnum*4+b2n[qbase]]) << endl;
	  }
	}
      }
	
      premutloc = mutlocref;
      ++rgapnum;
      ++gapnum;
	
    }

    else {

      if (prof[mutlocref*5+b2n[qbase]] != UCHAR_MAX) {
	++prof[mutlocref*5+b2n[qbase]];
      }
	
      if ('-' == qbase) ++qgapnum;

      premutloc = mutlocref+1;
      gapnum = 0;
    }

    preloc = loc + 1;
    loc = misstr.find(';', preloc);
  }

  updatebaseprof(premutloc, start+bestmap.len+qgapnum-rgapnum, ref2seq[refid], ref2prof[refid], ref2cov[refid]);
}


static void compute_depth(ifstream &ifs, S2D &ref2dep, S2S &ref2seq) {

  string eachline, refid;
  size_t pos1, pos2, len;
  while (getline(ifs, eachline)) {
    
    pos1 = eachline.find("\t")+1;
    pos2 = eachline.find("\t", pos1);
    refid = eachline.substr(pos1, pos2-pos1);

    pos1 = eachline.rfind("\t");
    pos2 = eachline.rfind("\t", pos1-1) + 1;
    len = atoi(eachline.substr(pos2, pos1-pos2).c_str());

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

static void compute_breadth(ifstream &ifs, S2D &ref2bre, S2S &ref2seq) {

  S2VB ref2cov;   // stores if each locus is covered
  init(ref2seq, ref2cov, 1); 

  string eachline, refid;
  size_t pos1, pos2, len, start;
  while (getline(ifs, eachline)) {
	
    pos1 = eachline.find("\t")+1;
    pos2 = eachline.find("\t", pos1);
    refid = eachline.substr(pos1, pos2-pos1);
    if (ref2cov.find(refid) == ref2cov.end()) continue;
	
    pos1 = eachline.find("\t", pos2+1) + 1;
    pos2 = eachline.find("\t", pos1);
    start = atoi(eachline.substr(pos1, pos2-pos1).c_str())-1;

    pos1 = eachline.rfind("\t");
    pos2 = eachline.rfind("\t", pos1-1) + 1;
    len = atoi(eachline.substr(pos2, pos1-pos2).c_str());

    fill(ref2cov[refid].begin()+start, ref2cov[refid].begin()+start+len, true);
  }

  for (S2VB::iterator iter = ref2cov.begin(); iter != ref2cov.end(); ++iter) {

    double cov = 0;
    for (size_t i = 0; i < iter->second.size(); ++i)
      if (iter->second[i])
	++cov;
    ref2bre.insert(S2D::value_type(iter->first, cov));
  }
}


void processmaps(const Cmdopt &cmdopt, S2S &ref2seq, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins) {


  init(b2n);
  
  ifstream ifs(cmdopt.mapfile.c_str());
  if (!ifs) {
    cerr << "Could not open file " << cmdopt.mapfile << endl;
    exit(1);
  }

  // use all map records
  if (cmdopt.pickref == "all") {

    // stores each map record
    Bestmap bestmap;
    size_t  pos1 = 0, pos2 = 0;
    string eachline;
    while (getline(ifs, eachline)) {
      
      pos1 = eachline.find("\t") + 1;
      pos2 = eachline.find("\t", pos1);
      bestmap.refid = eachline.substr(pos1, pos2-pos1);
    
      pos1 = eachline.find("\t", pos2+1) + 1;
      pos2 = eachline.find("\t", pos1);
      bestmap.start = atoi(eachline.substr(pos1, pos2-pos1).c_str())-1;

      pos1 = eachline.rfind("\t");
      bestmap.misstr = eachline.substr(pos1+1, eachline.size()-pos1-1);
    
      pos2 = eachline.rfind("\t", pos1-1) + 1;
      bestmap.len = atoi(eachline.substr(pos2, pos1-pos2).c_str());

      storemap(bestmap, ref2ins, ref2cov, ref2seq, ref2pos2ins, ref2prof);
    }
  }

  else {

    S2D ref2val;

    if (cmdopt.pickref == "depth")
      compute_depth(ifs, ref2val, ref2seq);

    else if (cmdopt.pickref == "breadth") 
      compute_breadth(ifs, ref2val, ref2seq);
    
    // rewind input file stream
    ifs.clear();
    ifs.seekg(0, ios::beg);


    string eachline, qid, preqid, refid;
    size_t pos1, pos2;
    double maxval = 0;
    Bestmap bestmap;
    while (getline(ifs, eachline)) {

      pos1 = eachline.find("\t") + 1;
      qid = eachline.substr(0, pos1);

      pos2 = eachline.find("\t", pos1);
      refid = eachline.substr(pos1, pos2-pos1);

      S2D::iterator iter = ref2val.find(refid);

      if (qid != preqid) {   // process previous record
	storemap(bestmap, ref2ins, ref2cov, ref2seq, ref2pos2ins, ref2prof);
	preqid = qid;
      }
      else {                // check if this ref is better
	if (iter == ref2val.end() || iter->second <= maxval)
	  continue;
      }
      maxval = iter->second;

      pos1 = eachline.find("\t", pos2+1) + 1;
      pos2 = eachline.find("\t", pos1);
      bestmap.start = atoi(eachline.substr(pos1, pos2-pos1).c_str())-1;

      pos1 = eachline.rfind("\t");
      bestmap.misstr = eachline.substr(pos1+1, eachline.size()-pos1-1);

      pos2 = eachline.rfind("\t", pos1-1) + 1;
      bestmap.len = atoi(eachline.substr(pos2, pos1-pos2).c_str());

      bestmap.refid = refid;
      
    }

    storemap(bestmap, ref2ins, ref2cov, ref2seq, ref2pos2ins, ref2prof);
  }
}


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

static void storesam(VC &prof, VB &ins, VB &cov, I2VC &pos2ins, string mutstr, string &seq, Uint start) {

  Uint n    = 0;      // tracks length of each region
  Uint locr = start;  // locus in ref
  Uint locq = 0;      // locus in query

  for (size_t i = 0; i < mutstr.size(); ++i) {

    char achar = mutstr[i];

    if (isdigit(achar))        // this is a digit
      n = n*10 + achar - '0';

    // this is a char, only analyze M: match, I: insertion, D: deletion
    else {

  

      if (achar == 'M') {
	//cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
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

	//cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	if (ins[locr]) { // > 1 times see this insertion

	  //cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;

	  if (pos2ins.find(locr) == pos2ins.end())
	    pos2ins.insert(I2VC::value_type(locr, VC(12, 0)));
	  
	  VC &insprof = pos2ins[locr];
	  
	  for (Uint i = 0; i < n; ++i) { // store max 3 insertions
	    //cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
	    if (i < 3) {
	      //cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << i*4+b2n[seq[locq]] << "\t" << mutstr << endl;
	      if (insprof[i*4+b2n[seq[locq]]] != UCHAR_MAX) ++insprof[i*4+b2n[seq[locq]]];
	      //cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << int(insprof[i*4+b2n[seq[locq]]]) << "\t" << mutstr << endl;
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
	//cout << start << "\t" << locq << "\t" << locr << "\t" << n << "\t" << mutstr << endl;
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

      storesam(ref2prof[refid], ref2ins[refid], ref2cov[refid], ref2pos2ins[refid], mutstr, seq, start);
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
    string eachline(""), refid(""), mutstr(""), seq(""), preqid(""), qid(""), prerefid("");
    Uint start = 0;
    double maxval = 0.0;
    while (getline(ifs, eachline)) {
      
      pos1 = eachline.find("\t");           // 1
      qid = eachline.substr(0, pos1);

      pos1 = eachline.find("\t", pos1+1);   // 2
      if (eachline[pos1-1] == '4') continue;// unmapped

      
      pos2 = eachline.find("\t", pos1+1);   // 3
      refid = eachline.substr(pos1+1, pos2-pos1-1);

      S2D::iterator iter = ref2val.find(refid);

      if (qid != preqid && !mutstr.empty()) // store previous record
	storesam(ref2prof[prerefid], ref2ins[prerefid], ref2cov[prerefid], ref2pos2ins[prerefid], mutstr, seq, start);
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
    storesam(ref2prof[refid], ref2ins[refid], ref2cov[refid], ref2pos2ins[refid], mutstr, seq, start);
  }
}
