#include <vector>
#include <cmath>
#include <unordered_set>
#include <algorithm> 
#include "helpers.h"
#include "feature_calculator.h"

using namespace std;

int feature_calculator::calculate_feature(contig::depth_feature f, bool compare_feature, bool compare_p){
	if(!cur_contig.valid || f.len == 0 || f.depth.size() == 0) return 0;

	int len = f.depth.size();
	bounds_check(len, 0, f.len - 1);

	vector<flagged_region> flags;
	flag_window window(f, compare_feature, compare_p);
	//get all of the flagged regions
	for(int i = 0; i < f.len; i++){
		bool is_end = i == f.len-1;
		if (window.update(i, i, f.p_vals[i], f.depth[i], is_end)){
			flags.push_back(flagged_region(window));
        }
	}
	if(flags.size() == 0) return 0; //no features here
	//merge flagged regions if they're close enough together

	vector<flagged_region> combined_flags;
	flagged_region cur_flag = flags[0];
	int n_flags = flags.size();
	for(int i = 1; i < n_flags; i ++){
		if(cur_flag.end + w >= flags[i].start){ //combine these two
			cur_flag.end = flags[i].end;
		}
		else{ //don't combine them
			combined_flags.push_back(cur_flag);
			cur_flag = flags[i];
		}
	}
	combined_flags.push_back(cur_flag);
	
	//get the feature count
	int n_feats = 0;
	for(flagged_region e: combined_flags){
		n_feats += max((int)round(((double)e.end - e.start)/w),1);
	}
	return n_feats;
}

void feature_calculator::calculate(){
	lo_doc = calculate_feature(cur_contig.doc, false, false);
	hi_doc = calculate_feature(cur_contig.doc);
	lo_mpdoc = calculate_feature(cur_contig.mpdoc, false, false);
	hi_mpdoc = calculate_feature(cur_contig.mpdoc);
	short_mpdoc = calculate_feature(cur_contig.mpdoc_short);
	long_mpdoc = calculate_feature(cur_contig.mpdoc_long);
	single_doc = calculate_feature(cur_contig.single_doc);
	miso_doc = calculate_feature(cur_contig.miso_doc);
}

vector<int> feature_calculator::get_features(){
	vector<int> feats{ lo_doc, hi_doc, lo_mpdoc, hi_mpdoc, short_mpdoc, long_mpdoc, single_doc, miso_doc };
	return feats;
}

feature_calculator::feature_calculator(contig &cur_contig, int window_size) : 
	cur_contig(cur_contig),
	w(window_size)
{}

