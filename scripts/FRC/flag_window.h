#ifndef FLAG_WINDOW_H
#define FLAG_WINDOW_H

#include <string>
#include "contig.h"

class flag_window {
public: 
	std::string tag; //The name of the feature that this window is tracking. used for output later
	std::string cname; //The name of the contig
	double cutoff; //The significance cutoff for this feature
	int left; //The lower index of the window
	int right; //The higher index of the window
	int upper_edge_cutoff; //The higher edge cutoff for this feature
	int lower_edge_cutoff; //The lower edge cutoff for this feture
	int minimum_value; //the smallest feature value that we will consider significant
	bool compare_p; //True if we are looking for higher p values, False if we are looking for lower ones
	bool compare_feature; //True if we are looking for the largest feature value, false if we are looking for the smallest
	bool is_open; //Keeps track of whether or not this window is currently capturing a significant region. True if it contains significant values and false if it's just searching
	int best_feature; //The max (or min) value of the feature in this window. TODO should this be a double?
	double best_p; //The most significant p value in this window
	
	int best_feature_index; //The index of the best feature in this window
	int best_p_index; //The index of the most significant p value in this window

	flag_window (contig::depth_feature f, bool compare_feature = true, bool compare_p = true);
	flag_window (contig::event_feature f, bool compare_feature = true, bool compare_p = true);
	bool update(int cur_left, int cur_right, double p, int feature, bool force_update = false);
	std::string to_out_str();
};

#endif