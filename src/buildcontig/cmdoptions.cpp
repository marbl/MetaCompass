#include <iostream>
using std::endl;
using std::cerr;

#include <cstring>
using std::strcmp;

#include <cstdlib>
using std::atoi;

#include "cmdoptions.hpp"


// print out help information
void buildcontighelp() {
  cerr << endl;

  cerr << "Usage:" << endl;
  cerr << "        ./buildcontig -m/s [read mapping file] -r [reference seq file] -o [output directory]" << endl;
  cerr << endl;

  cerr << "Options:" << endl;
  cerr << "        -m/--mapfile  read mapping file from mummer-map" << endl;
  cerr << "        -s/--samfile  read mapping file from bowtie2, BWA, etc." << endl;
  cerr << "        -r/--refseq   reference sequences used to guide genome assembly" << endl;
  cerr << "        -o/--prefix   output files prefix" << endl;
  cerr << "        -c/--mincov   minimum depth of coverage for contigs (default: 2)" << endl;
  cerr << "        -l/--minlen   minimum length for contigs (default: 100bp)" << endl;
  cerr << "        -k/--pickref  pick reference if there are multiple best ones (default: breadth)" << endl;
  cerr << "                      all      - use all of them; will produce redundant contigs" << endl;
  cerr << "                      breadth  - pick reference with highest breadth of coverage" << endl;
  cerr << "                      depth    - pick reference with highest depth of coverage" << endl;
//  cerr << "                      random   - randomly pick one" << endl;
  cerr << "        -b/--prof     output(T) or not(F) .baseprof file (default: F)" << endl;
  cerr << "        -n/--newref   output(T) or not(F) .newref file (default: F)" << endl;
  cerr << "        -u/--usedmap  output(T) or not(F) .usedmap file (default: F)" << endl;

  cerr << "Contact:" << endl;
  cerr << "        Have problems? Victoria Cepeda - vcepeda@umiacs.umd.edu" << endl;
  cerr << endl;
}



// parse command line options
void parsecmdoptions(int argc, char *argv[], Cmdopt* cmdopt) {


  if (argc < 7) {                              // no enough arguments
    buildcontighelp();
    exit(1);
  }


  bool tag = false;
  for (int i = 1; i < argc; i += 2) {        // read each option
    
    if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mapfile") == 0) {          // read mapping file
      cmdopt->mapfile = argv[i+1];
      cmdopt->filetype = "map";
    }
    
    else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--samfile") == 0) {          // read mapping file
      cmdopt->mapfile = argv[i+1];
      cmdopt->filetype = "sam";
    }
    
    else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--refseq") == 0)     // ref seq file
      cmdopt->reffile = argv[i+1];
    
    else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--prefix") == 0)     // output prefix
      cmdopt->outprefix = argv[i+1];

    else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--mincov") == 0) {   // min depcov for outputing contigs
      if (atoi(argv[i+1]) < 1)
	tag = true;
      cmdopt->mindepcov = atoi(argv[i+1]);
    }

    else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--minlen") == 0) {   // min length for outputing contigs
      if (atoi(argv[i+1]) < 1)
	tag = true;
      cmdopt->minlen = atoi(argv[i+1]);
    }

    else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--thread") == 0) {   // min length for outputing contigs
      if (atoi(argv[i+1]) < 1)
	tag = true;
      cmdopt->nump = atoi(argv[i+1]);
    }

    else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--prof") == 0) {   // print .baseprof file or not
      if (argv[i+1][0] == 'T') 
	cmdopt->printbaseprof = true;
      else if (argv[i+1][0] == 'F') 
	cmdopt->printbaseprof = false;
      else
	tag = true;
    }

    else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--newref") == 0) {   // print .newref file or not
      if (argv[i+1][0] == 'T') 
	cmdopt->printnewref = true;
      else if (argv[i+1][0] == 'F') 
	cmdopt->printnewref = false;
      else
	tag = true;
    }

    else if (strcmp(argv[i], "-u") == 0 || strcmp(argv[i], "--usedmap") == 0) {   // print .usedmap file or not
      if (argv[i+1][0] == 'T') 
	cmdopt->printusedmap = true;
      else if (argv[i+1][0] == 'F') 
	cmdopt->printusedmap = false;
      else
	tag = true;
    }

    else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--pickref") == 0) {   // how to pick ref seq
      if (strcmp(argv[i+1], "all") == 0) 
	cmdopt->pickref = "all";
      else if (strcmp(argv[i+1], "breadth") == 0) 
	cmdopt->pickref = "breadth";
      else if (strcmp(argv[i+1], "random") == 0) 
	cmdopt->pickref = "random";
      else if (strcmp(argv[i+1], "depth") == 0) 
	cmdopt->pickref = "depth";
      else
	tag = true;
    }
    
    else                                     // does not match, then something goes wrong
      tag = true;
  }


  // these parameters must be specified
  if (cmdopt->mapfile == "" || cmdopt->reffile == "" || cmdopt->outprefix == "")
    tag = true;

  if (tag) {
    buildcontighelp();
    exit(1);
  }
}
