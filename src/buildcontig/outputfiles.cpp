#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <algorithm>
#include "outputfiles.hpp"


const unsigned extendbp = 100;    // use upstream and downstream extendbp bases for .newref
char num2base[] = {'A', 'C', 'G', 'T', '-'};

void createNewref(Cmdopt &cmdopt, S2S &ref2seq, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins) {

  ofstream ofs(string(cmdopt.outprefix + "/newref.fasta").c_str());

  // go through each reference genome
  for (S2VC::iterator iter = ref2prof.begin(); iter != ref2prof.end(); ++iter) {

    // stores temporary contig sequence
    string contig("");
    string refid = iter->first;
    VC &vc = iter->second;
    VB &cov = ref2cov[refid];
    VB &ins = ref2ins[refid];
    bool inctg = true;
    string &refseq = ref2seq[refid];
    Uint ctgn = 0;
    Uint ctgstart = 0, ctgend = 0; // start and end of previous contig

    // iterate each locus of reference genome
    for (size_t i = 0; i < vc.size(); i += 5) {

      size_t pos = i/5;

      
      // compute depth of coverage and base with max depcov
      /*******************************************************************/
      Uint max = vc[i]; // max depth of coverage for a base
      Uint maxj = 0;    // base with max depcov
      Uint depcov = vc[i];  // total depcov
      if (!cov[pos]) max = 0; // no read aligned to this locus
      else {
	for (size_t j = 1; j < 5; ++j) {
	  if (vc[i+j] > max) {
	    max = vc[i+j];
	    maxj = j;
	  }
	  depcov += vc[i+j];
	}
      }
      /*******************************************************************/

      // depcov is high enough
      if (max >= cmdopt.mindepcov) {

	// gap in query reference, do nothing
	if (maxj == 4) {
	  //cout << pos << endl;
	  continue;
	}

	if (!inctg) {            // not in contig, start of a new contig
	  
	  if (pos - ctgend - 1 > extendbp * 2) {         // gap is too big

	    contig += refseq.substr(ctgend+1, extendbp); // extend previous contig
	    if (contig.size() >= cmdopt.minlen) {
	      ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << ctgend+extendbp << endl;
	      ofs << contig << endl;
	    }

	    // create new contig
	    contig = refseq.substr(pos - extendbp, extendbp);
	    ctgstart = pos - extendbp;
	  }

	  // gap is not big
	  else {
	    if (ctgend == 0)
	      contig += refseq.substr(ctgend, pos - ctgend);
	    else
	      contig += refseq.substr(ctgend+1, pos - ctgend - 1);
	  }
	}
	inctg = true;
	ctgend = pos;
	contig.push_back(num2base[maxj]);
	
	//cout << contig << "\t" << pos << endl;

	// if there are insertions at this locus
	/********************************************************/
	if (ins[pos] == true) {
	  
	  S2I2VC::iterator iter_ins = ref2pos2ins.find(refid);
	  if (iter_ins == ref2pos2ins.end()) continue;

	  I2VC::iterator iter_prof = iter_ins->second.find(pos);
	  if (iter_prof == iter_ins->second.end()) continue;

	  VC &insprof = iter_prof->second;

	  for (size_t k = 0; k < 3; ++k) {

	    Uint max = insprof[k*4];
	    Uint maxj = 0;
	    for (size_t j = 1; j < 4; ++j) {
	      if (insprof[k*4+j] > max) {
		max = insprof[k*4+j];
		maxj = j;
	      }
	    }

	    if (max >= depcov / 2) {
	      contig.push_back(num2base[maxj]);
	    }
	    else break;
	  }
	}
	/********************************************************/
	
      }
      
      // this is a gap between contig
      else
	inctg = false;
    }

    // last contig
    if (contig.size() >= cmdopt.minlen) {
      ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << ctgend << endl;
      ofs << contig << endl;
    }

  }
  
}




void createContig(Cmdopt &cmdopt, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins) {

  ofstream ofs(string(cmdopt.outprefix + "/contigs.fasta").c_str());

  for (S2VC::iterator iter = ref2prof.begin(); iter != ref2prof.end(); ++iter) {

    Uint ctgn = 0;
    Uint ctgstart = 0;
    string contig("");
    string refid = iter->first;
    VC &vc = iter->second;
    VB &ins = ref2ins[refid];
    VB &cov = ref2cov[refid];
    
    for (size_t i = 0; i < vc.size(); i += 5) {

      size_t pos = i/5;

      
      Uint max = vc[i];
      Uint maxj = 0;
      Uint depcov = vc[i];

      if (!cov[pos])
	max = 0;
      else {
	for (size_t j = 1; j < 5; ++j) {
	  
	  if (vc[i+j] > max) {
	    max = vc[i+j];
	    maxj = j;
	  }
	  depcov += vc[i+j];
	}
      }

      //cout << max << " " << cmdopt.mindepcov << endl;
      if (max >= cmdopt.mindepcov) {
      
	if (maxj == 4) continue; // gap in query reference

	contig.push_back(num2base[maxj]);

	if (ins[pos]) {
	  
	  S2I2VC::iterator iter_ins = ref2pos2ins.find(refid);
	  if (iter_ins == ref2pos2ins.end()) continue;

	  I2VC::iterator iter_prof = iter_ins->second.find(pos);
	  if (iter_prof == iter_ins->second.end()) continue;

	  VC &insprof = iter_prof->second;

	  for (size_t k = 0; k < 3; ++k) {

	    Uint max = insprof[k*4];
	    Uint maxj = 0;
	    for (size_t j = 1; j < 4; ++j) {
	      if (insprof[k*4+j] > max) {
		max = insprof[k*4+j];
		maxj = j;
	      }
	    }

	    if (max > depcov / 2)
	      contig.push_back(num2base[maxj]);
	    else break;
	  }
	}
      }

      else {

	//cout << i << "\t" << contig.size() << endl;
	if (contig.size() >= cmdopt.minlen) {
	  ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << pos << endl;
	  ofs << contig << endl;
	}
	if (!contig.empty())
	  contig.clear();
	ctgstart = pos+1; // +1: 0 to 1; +1: current is gap, next is true start
      }
    }

    if (contig.size() >= cmdopt.minlen) {
      ofs << ">" << refid << "_" << ctgn++ << " " << ctgstart << " " << vc.size()/5 << endl;
      ofs << contig << endl;
    }
  }
}




void printbaseprof(Cmdopt &cmdopt, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins) {


  ofstream ofs(string(cmdopt.outprefix + "/baseprof").c_str());
  
  for (S2VC::iterator iter = ref2prof.begin(); iter != ref2prof.end(); ++iter) {
    string refid = iter->first;
    VC &vc = iter->second;
    VB &ins = ref2ins[refid];
    VB &cov = ref2cov[refid];
    
    if (refid.empty()) continue;

    ofs << ">" << refid << endl;

    //cout << vc.size() << endl;

    for (size_t i = 0; i < vc.size(); i += 5) {

      size_t pos = i/5;

      if (cov[pos]) {

	ofs << pos << "\t";
	for (size_t j = 0; j < 5; ++j) {
	  ofs << vc[i+j]*1 << "\t";
	}
	ofs << endl;



	if (ins[pos] == true) {


	  S2I2VC::iterator iter_ins = ref2pos2ins.find(refid);
	  if (iter_ins == ref2pos2ins.end()) continue;

	  I2VC::iterator iter_prof = iter_ins->second.find(pos);
	  if (iter_prof == iter_ins->second.end()) continue;

	  VC &insprof = iter_prof->second;

	  for (size_t k = 0; k < 3; ++k) {

	    bool tag = false;

	    for (size_t j = 0; j < 4; ++j) {
	      if (insprof[k*4+j] > 0) {
		tag = true;
		break;
	      }
	    }

	    if (tag) {
	      ofs << pos << "." << k << "\t";
	      for (size_t j = 0; j < 4; ++j) {
		ofs << insprof[k*4+j]*1 << "\t";
	      }
	      ofs << endl;
	    }
	    else break;
	  
	  }

	}
      }
    }
  }
}












