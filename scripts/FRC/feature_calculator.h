#ifndef FEATURE_CALCULATOR_H
#define FEATURE_CALCULATOR_H

#include <iostream>
#include <vector>
#include "contig.h"
#include "flagged_region.h"

class interval{
public:
	int id;
	int index;

	interval(int index, int id) :
		id(id),
		index(index)
	{}
	
	bool operator <(const interval & cur_interval) const {
		return index < cur_interval.index;
	}
};

class feature_calculator{
	
	contig cur_contig; //contig that's being flagged
	int w; //feature window size
	int s; //step size

	//total number of basepairs in regions containing each feature
	int lo_doc;
	int hi_doc;
	int lo_mpdoc;
	int hi_mpdoc;
	int short_mpdoc;
	int long_mpdoc;
	int single_doc;
	int miso_doc;


public:
	std::vector<int> get_features();
	feature_calculator(contig &cur_contig, int window_size /*cutoffs*/);		
	int calculate_feature(contig::depth_feature f, bool compare_feature = true, bool compare_p = true);
	void calculate();
};

#endif
