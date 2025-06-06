#include "flagged_region.h"
#include "helpers.h"

using namespace std;

flagged_region::flagged_region(flag_window &f) :
	tag(f.tag),
	cname(f.cname),
	start(f.left),
	end(f.right),
	best_value(f.best_feature),
	best_value_index(f.best_feature_index),
	best_p(f.best_p),
	best_p_index(f.best_p_index)
{}

string flagged_region::to_str(){
	string ret = cname + " " + tag + " " + to_string(start) + " " + to_string(end) + " " + to_string(best_p) + " " + to_string(best_p_index) + " " + to_string(best_value) + " " + to_string(best_value_index);

	int len = links.size();
	if(len == 0) return ret;
	
	bounds_check(len, 0, len - 1);
	ret += "\nAssociated_links:\n";
	for(int i = 0; i < len; i ++){
		ret += links[i].to_bed_str() + "\n";
	}
	return ret;
}