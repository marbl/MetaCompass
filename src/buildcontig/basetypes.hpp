#ifndef BASETYPES_HPP
#define BASETYPES_HPP

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;


typedef unsigned char         Ucha;
typedef vector<Ucha>          VC;
typedef map<string, VC>       S2VC;
typedef unsigned int          Uint;
typedef map<string, string>   S2S;
typedef map<string, double>   S2D;
typedef vector<bool>          VB;
typedef map<string, VB>       S2VB;
typedef map<Uint, VC>         I2VC;
typedef map<string, I2VC>     S2I2VC;

class Cmdopt {
public:
  string mapfile,
         reffile,
         outprefix,
         pickref,
         filetype;
  
  Uint   mindepcov,       // minimum depth of coverage of contig
         minlen,          // mininum contig length
         nump;
  
  bool   printcontig,   // print contigs.fasta or not
         printbaseprof,   // print .baseprof or not
         printnewref,     // print .newref or not
         printusedmap;    // print .usedmap or not
  
  Cmdopt() : pickref("breadth"), mindepcov(2), minlen(100), nump(1), printcontig(true),
	     printbaseprof(false), printnewref(false), printusedmap(false) {};
};
  
#endif
