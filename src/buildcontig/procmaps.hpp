#include "basetypes.hpp"

// process sam file
void processsams(const Cmdopt &cmdopt, S2S &ref2seq, S2VC &ref2prof, S2VB &ref2ins, S2VB &ref2cov, S2I2VC &ref2pos2ins);

// Read reference sequence fasta file.
void readrefseqfile(const string refseqfile, S2S &ref2seq);

struct Bestmap {
  string refid,
         misstr;
  Uint   start;
  Uint   len;
  Bestmap() : refid(""), misstr(""), start(0), len(0) {}
};
