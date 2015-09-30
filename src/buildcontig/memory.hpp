#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "basetypes.hpp"

// initialze (allocate memory) various containers (data structures)

void init(const S2S &ref2seq,
	  S2VC &ref2vecint,
	  Uint num);

void init(const S2S &ref2seq,
	  S2VB &ref2vecbool,
	  Uint num);

void init(const S2S &ref2seq,
	  S2I2VC &ref2pos2ins);

#endif
