#ifndef FLAGGED_REGION_H
#define FLAGGED_REGION_H

#include <string>
#include "flag_window.h"
#include "mate_pair.h"

class flagged_region {
public:

	std::string tag;
	std::string cname;
	int start;
	int end;
	int best_value;
	int best_value_index;
	double best_p;
	int best_p_index;
	std::vector<mate_pair> links;

	flagged_region(flag_window &f);
	std::string to_str();
};

#endif