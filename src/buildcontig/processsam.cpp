#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

#include <cstdlib>
#include <time.h>

#include "cmdoptions.hpp"   // command line options and usage
#include "basetypes.hpp"    // data types and structures
#include "memory.hpp"       // STL container initialization
#include "outputfiles.hpp"  // output function
#include "procmaps.hpp"     // process and store the read maps


int main (int argc, char *argv[]) {
  

  // parse command line options
  Cmdopt cmdopt;
  parsecmdoptions(argc, argv, &cmdopt);

  

  // read reference sequences
  S2S ref2seq;
  time_t nowtime = time(NULL);
  cerr << endl << "# Reading reference fasta file ... ";
  readrefseqfile(cmdopt.reffile, ref2seq);
  cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
  

  
  // compute mapping profiles of reference genomes
  /*******************************************************************/
  nowtime = time(NULL);
  cerr << "# Computing profiles for genomes ... ";

  S2VC ref2prof;              // stores the abundance of ACGT- at each locus
  init(ref2seq, ref2prof, 5); // 5 chars per locus: ACGT-
  
  S2VB ref2ins,               // if there is an insert at each locus
       ref2cov;               // if each locus is covered
  init(ref2seq, ref2ins, 1);  // 1 bit per locus
  init(ref2seq, ref2cov, 1); 
  
  S2I2VC ref2pos2ins;         // hashtable keeps inserts for each ref.
  init(ref2seq, ref2pos2ins); 

  processsams(cmdopt, ref2seq, ref2prof, ref2ins, ref2cov, ref2pos2ins);
  cerr << "success!" << " (" << time(NULL) - nowtime << " seconds)" << endl;
  /**********************************************************************/

  
  return 0;
}


// to do:
// use const iterator
// implement random ref. selection

