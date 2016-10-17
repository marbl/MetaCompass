#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <new>
using std::bad_alloc;

#include <cstdlib>
#include "memory.hpp"


void init(const S2S &ref2seq,
	  S2VC &ref2vint,
	  Uint num) {

  for (S2S::const_iterator citer = ref2seq.begin();
       citer != ref2seq.end();
       ++citer) {

    Uint len = citer->second.size();       // length

    try {
      ref2vint.insert( S2VC::value_type(citer->first, VC((len)*num, 0)) );
    }
    catch (bad_alloc& ba) {
      cerr << endl << endl;
      cerr << "Something wrong when allocating memory during allocateMem()"<< endl;
      cerr << "bad_alloc caught: " << ba.what() << endl << endl;
      exit(1);
    }
    
  }
}


void init(const S2S &ref2seq,
	  S2VB &ref2vbool,
	  Uint num) {

  for (S2S::const_iterator citer = ref2seq.begin();
       citer != ref2seq.end();
       ++citer) {

    Uint len = citer->second.size();       // length

    try {
      ref2vbool.insert( S2VB::value_type(citer->first, VB((len)*num, false)) );
    }
    catch (bad_alloc& ba) {
      cerr << endl << endl;
      cerr << "Something wrong when allocating memory during allocateMem()"<< endl;
      cerr << "bad_alloc caught: " << ba.what() << endl << endl;
      exit(1);
    }
    
  }
}

void init(const S2S &ref2seq,
	  S2I2VC &ref2pos2ins) {

  for (S2S::const_iterator citer = ref2seq.begin();
       citer != ref2seq.end();
       ++citer)
    ref2pos2ins.insert( S2I2VC::value_type(citer->first, I2VC()) );
}
