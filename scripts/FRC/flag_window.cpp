
#include "flag_window.h"
using namespace std;

flag_window::flag_window (contig::depth_feature f, bool compare_feature, bool compare_p): 
	tag(f.tag), 
	cname(f.cname),
	cutoff(f.flag_cutoff), 
	left(-1), 
	right(-1), 
	upper_edge_cutoff(f.upper_edge_cutoff),
	lower_edge_cutoff(f.lower_edge_cutoff),
	minimum_value(f.minimum_value),
	compare_p(compare_p), 
	compare_feature(compare_feature), 
	is_open(false), 
	best_feature(-1), 
	best_p(-1)
{}
	
flag_window::flag_window (contig::event_feature f, bool compare_feature, bool compare_p): 
	tag(f.tag), 
	cname(f.cname),
	cutoff(f.flag_cutoff), 
	left(-1), 
	right(-1), 
	upper_edge_cutoff(f.upper_edge_cutoff),
	lower_edge_cutoff(f.lower_edge_cutoff),
	minimum_value(f.minimum_value),
	compare_p(compare_p), 
	compare_feature(compare_feature), 
	is_open(false), 
	best_feature(-1), 
	best_p(-1)
	//TODO should best_feature_index, etc also be -1?
{}

bool flag_window::update(int cur_left, int cur_right, double p, int feature, bool force_update){		
	bool significant = (compare_p ? p >= cutoff : p <= -cutoff) && feature >= minimum_value; //TODO find amore elegant solution to this, like passing a func
	significant = significant && (cur_left > lower_edge_cutoff && cur_right < upper_edge_cutoff);

	if(!is_open && !significant) return false; //the window is not open and we don't need to open it

	//the window is not open and we need to open it here 
	//TODO deal with overlaps here?
	else if(!is_open && significant){ 
		best_p = p;
		best_p_index = cur_left;
		best_feature = feature;
		best_feature_index = cur_left;
		left = cur_left;
		right = cur_right;
		is_open = true;
		return force_update;
	}

	else if(is_open && significant){ //the window is open and we need to extend it
		right = cur_right;

		if((compare_p && p > best_p) || (!compare_p && p < best_p)){
			best_p = p;
			best_p_index = cur_left;
		}
		if((compare_feature && feature > best_feature) || (!compare_feature && feature < best_feature)){
			best_feature = feature;
			best_feature_index = cur_left;
		}
		return force_update;
	} 
	//else if(is_open && !significant){ //the window is open and we need to close it
	is_open = false;
	
	return true;
}

string flag_window::to_out_str(){
	return tag + " " + to_string(left) + " " + to_string(right) + " " + to_string(best_p) + " " + to_string(best_p_index) + " " + to_string(best_feature) + " " + to_string(best_feature_index);
}
