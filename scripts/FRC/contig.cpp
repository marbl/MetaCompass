
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream> //TODO debug only, remove
#include <algorithm> //for sort
#include "helpers.h"
#include "contig.h"


using namespace std;

//Depth feature functions
contig::depth_feature::depth_feature(string tag, string cname, int len, int lower_edge_cutoff, int upper_edge_cutoff, double flag_cutoff, bool contig_valid, int minimum_value): 
	tag(tag),
	cname(cname),
	len(len),
	read_count(0),
	excluded(0),
	lower_edge_cutoff(lower_edge_cutoff),
	upper_edge_cutoff(upper_edge_cutoff),
	minimum_value(minimum_value),
	bounded_len(upper_edge_cutoff - lower_edge_cutoff),
	bounded_avg(-1),
	bounded_sdv(-1),
	flag_cutoff(flag_cutoff),
	depth(vector<int>(len+1)), 
	p_vals(vector<double>(len+1)),
	perbase_avg(vector<double>())
{
	valid = contig_valid && len > 0 && lower_edge_cutoff < upper_edge_cutoff && lower_edge_cutoff >= 0 && upper_edge_cutoff <= len;
}

contig::depth_feature::depth_feature():
	len(-1),
	read_count(0),
	excluded(0),
	lower_edge_cutoff(0),
	upper_edge_cutoff(0),
	minimum_value(0),
	bounded_len(0),
	bounded_avg(0),
	bounded_sdv(0),
	flag_cutoff(0),
	depth(vector<int>(0)), 
	p_vals(vector<double>(0)),
	perbase_avg(vector<double>(0)),
	valid(false)
{
	tag = string("");
	cname = string("");
}

void contig::depth_feature::update(mate_pair& mp, int value){
	if(!valid) return;

	bool in_bounds = mp.check_edge(lower_edge_cutoff, upper_edge_cutoff);
	if(in_bounds) read_count++;
	else excluded++;

	bounds_check(depth.size(), mp.a.start, mp.b.end);
	depth[mp.a.start] += value;
	depth[mp.b.end] -= value;
}

void contig::depth_feature::update(single_read& r, int value){
	if(!valid) return;

	bool in_bounds = r.check_edge(lower_edge_cutoff, upper_edge_cutoff);
	if(in_bounds) read_count++;
	else excluded++;

	bounds_check(depth.size(), r.start, r.end);
	depth[r.start] += value;
	depth[r.end] -= value;		
}

void contig::depth_feature::cumulative_sum(){
	if(!valid) return;

	bounds_check(depth.size(), 1, len);
	for(int i = 1; i < len; i++) depth[i] += depth[i-1];
}

void contig::depth_feature::get_bounded_avg_sdv(){
	if(!valid) return;

	if(upper_edge_cutoff <= lower_edge_cutoff || bounded_len <= 0) return;
	bounded_avg = 0.0;

	bounds_check(depth.size(), lower_edge_cutoff, upper_edge_cutoff-1);
	for(int i = lower_edge_cutoff; i < upper_edge_cutoff; i++){
		bounded_avg += depth[i];
	}

	bounded_avg /= bounded_len;
	bounded_sdv = 0.0;

	for(int i = lower_edge_cutoff; i < upper_edge_cutoff; i++){
		bounded_sdv += (depth[i]-bounded_avg)*(depth[i]-bounded_avg);
	}
	bounded_sdv /= bounded_len;
	bounded_sdv = sqrt(bounded_sdv);
}

//get the probability of the depth of coverage at a base using z scores, given the average and stdev of the doc on the contig
//this will flag too low or too high doc because, though the test is one sided, it uses the absolute value of the zscore
double contig::depth_feature::get_doc_p(double avg_doc, double sdv_doc, int doc){
	if(!valid) return 0;

	boost::math::normal sd(0.0, 1.0);
	if(avg_doc == 0 || sdv_doc == 0){ return 0; } //TODO do something about this case. What does this mean?? 
	double zscore = (doc - avg_doc)/sdv_doc;
	double p = -log10(1.0 - boost::math::cdf(sd, abs(zscore)));
	if(doc < avg_doc) p *= -1.0; //the negative is so that we can tell what was a positive z score and what was a negative z score (less than expected or more than expected)
	return p;
}

void contig::depth_feature::get_pvals_zscore(){
	if(!valid) return;

	if(lower_edge_cutoff != 0){
		bounds_check(p_vals.size(), 0, lower_edge_cutoff - 1);
		for(int i = 0; i < lower_edge_cutoff; i ++){
			p_vals[i] = 0;
		}
	}

	bounds_check(p_vals.size(), lower_edge_cutoff, upper_edge_cutoff - 1);
	for(int i = lower_edge_cutoff; i < upper_edge_cutoff; i++){
		p_vals[i] = get_doc_p(bounded_avg, bounded_sdv, depth[i]);
	}

	if(upper_edge_cutoff != len){
		bounds_check(p_vals.size(), upper_edge_cutoff, len - 1);
		for(int i = upper_edge_cutoff; i < len; i ++){
			p_vals[i] = 0;
		}
	}
	//old code for reference:
	//double p_doc = get_doc_p(c.doc.depth[i], c.doc.bounded_avg, c.doc.bounded_sdv);
	//double p_mpdoc = get_doc_p(c.mpdoc.depth[i], c.mpdoc.bounded_avg, c.mpdoc.bounded_sdv); 
	//double p_miso = c.miso_doc.depth[i] >= miso_cutoff ? 1.0 : 0.0; //TODO do this better
	//double p_single = c.single_doc.depth[i] >= single_cutoff ? 1.0 : 0.0; //TODO find a better method!
}

//TODO pass pval function as param to feature so we can use different ones for different features!
void contig::depth_feature::finish(){
	if(!valid) return;

	cumulative_sum();
	get_bounded_avg_sdv();
	get_pvals_zscore();
}

string contig::depth_feature::to_str(){
	string ret = tag + " " + cname + "\nlen: " + to_string(len) + "\nread_count: " + to_string(read_count) + "\nexcluded: " + to_string(excluded) + "\nlower_edge_cutoff: " + to_string(lower_edge_cutoff) + "\nupper_edge_cutoff: " + to_string(upper_edge_cutoff) + "\nbounded_len: " + to_string(bounded_len) + "\nbounded_avg: " + to_string(bounded_avg) + "\nflag_cutoff: " + to_string(flag_cutoff) + "\n";
	if(!valid) return ret;
	bounds_check(depth.size(), 0, len - 1);
	bounds_check(p_vals.size(), 0, len - 1);
	//for(int i = 0; i < len; i ++) ret += to_string(i) + " " + to_string(depth[i]) + " " + to_string(p_vals[i]) + "\n";
	return ret;
}


//Event veature functions

contig::event_feature::event_feature(string tag, string cname, int len, int lower_edge_cutoff, int upper_edge_cutoff, int smoothing_window_len, double flag_cutoff, bool contig_valid, int minimum_value):
	tag(tag),
	cname(cname),
	len(len),
	feature_events(vector<int>(len+1)),
	smoothed_feature_events(vector<int>(len+1)),
	total_events(vector<int>(len+1)), //will be filled in during set_totals. TODO make this a reference
	smoothed_total_events(vector<int>(len+1)),
	feature_event_count(0),
	total_event_count(0),
	excluded(0),
	lower_edge_cutoff(lower_edge_cutoff),
	upper_edge_cutoff(upper_edge_cutoff),
	minimum_value(minimum_value),
	smoothing_window_len(smoothing_window_len),
	lower_edge_smoothed_cutoff(lower_edge_cutoff + smoothing_window_len),
	flag_cutoff(flag_cutoff),
	p_vals(vector<double>(len+1)),
	prop(0.0)
{
	valid = contig_valid && len > 0 && lower_edge_cutoff < upper_edge_cutoff && lower_edge_cutoff >= 0 && upper_edge_cutoff <= len && lower_edge_smoothed_cutoff < upper_edge_cutoff && lower_edge_smoothed_cutoff < len; //TODO is the end inclusive or exclusive?
}

contig::event_feature::event_feature() : 
	len(-1),
	feature_events(vector<int>()),
	smoothed_feature_events(vector<int>(len+1)),
	total_events(vector<int>()), //will be filled in during set_totals. TODO make this a reference
	smoothed_total_events(vector<int>()),
	feature_event_count(0),
	total_event_count(0),
	excluded(0),
	lower_edge_cutoff(0),
	upper_edge_cutoff(0),
	minimum_value(0),
	smoothing_window_len(0),
	lower_edge_smoothed_cutoff(0),
	flag_cutoff(0),
	p_vals(vector<double>()),
	prop(0.0),
	valid(false)
{
	tag = string("");
	cname = string("");
}

void contig::event_feature::update(mate_pair& mp, int value){
	if(!valid) return;
	bool in_bounds = mp.check_edge(lower_edge_cutoff, upper_edge_cutoff);
	if(in_bounds) feature_event_count++;
	else excluded++;

	bounds_check(feature_events.size(), mp.a.start, mp.b.end);
	feature_events[mp.a.start]++;
	feature_events[mp.b.end]++;
}

void contig::event_feature::update(single_read& r, int lower_edge_cutoff, int upper_edge_cutoff){
	if(!valid) return;
	bool in_bounds = r.check_edge(lower_edge_cutoff, upper_edge_cutoff);
	if(in_bounds) feature_event_count++;
	else excluded++;

	bounds_check(feature_events.size(), r.start, r.end);
	feature_events[r.start]++;
	feature_events[r.end]++;		
}

//TODO check how bounds and edge cutoffs are being handled here, this is almost certainly wrong
//TODO make this two calls to smooth to get rid of code duplication here
void contig::event_feature::smooth(){
	if(!valid) return;
	int feature_window_count = 0;
	int total_window_count = 0;

	if(lower_edge_smoothed_cutoff != lower_edge_cutoff){
		//calculate the first window
		bounds_check(feature_events.size(), lower_edge_cutoff, lower_edge_smoothed_cutoff - 1);
		bounds_check(total_events.size(), lower_edge_cutoff, lower_edge_smoothed_cutoff - 1);
		for(int i = lower_edge_cutoff; i < lower_edge_smoothed_cutoff; i ++){
			feature_window_count += feature_events[i];
			total_window_count += total_events[i];
		}
	}

	//calculate the rest of the windows
	bounds_check(feature_events.size(), lower_edge_smoothed_cutoff, upper_edge_cutoff - 1);
	bounds_check(total_events.size(), lower_edge_smoothed_cutoff, upper_edge_cutoff - 1);
	bounds_check(feature_events.size(), lower_edge_smoothed_cutoff - smoothing_window_len, upper_edge_cutoff - smoothing_window_len - 1);
	bounds_check(total_events.size(), lower_edge_smoothed_cutoff - smoothing_window_len, upper_edge_cutoff - smoothing_window_len - 1);
	for(int i = lower_edge_smoothed_cutoff; i < upper_edge_cutoff; i++){
		feature_window_count += feature_events[i];
		total_window_count += total_events[i];
		int window_left = i-smoothing_window_len;
		feature_window_count -= feature_events[window_left];
		total_window_count -= total_events[window_left];

		smoothed_feature_events[i] = feature_window_count;
		smoothed_total_events[i] = total_window_count;
	}

}

//get the probability of observing this many specific mate pairs out of the total number, given the overall proporion of that type to the total
//This won't flag having too few unhappy mps because the test is one sided, we are getting the probability of seeing x or more events
double contig::event_feature::get_mp_p(int observed, int total){
	if(!valid) return 0;

	if(prop == 0 || total == 0) return 1; //TODO cut out of this case earlier or something. 
	double expected = prop*total;
	boost::math::poisson_distribution<> poisson(expected);

	//TODO these two statements should be equivalent so check!! If not my understanding of pdf is wrong
	//I don't know what to do with this, it doesn't work without the pdf and I don't know why
	//return -log10(1.0 - boost::math::cdf(poisson, observed) + boost::math::pdf(poisson, observed));
	//return -log10(1.0 - boost::math::cdf(poisson, observed)); //why does this return 0 instead of some epsilon
	
	return -log10(boost::math::pdf(poisson, observed));
}

//TODO think about how edges/cutoffs are handled. Shuold they be handled here to reduce issues with poisson, or should they be handled when flagging?
void contig::event_feature::get_pvals_poisson(){
	if(!valid) return;

	bounds_check(p_vals.size(), 0, lower_edge_smoothed_cutoff - 1);
	for(int i = 0; i < lower_edge_smoothed_cutoff; i ++){
		p_vals[i] = 0;
	}
	bounds_check(p_vals.size(), lower_edge_smoothed_cutoff, upper_edge_cutoff - 1);
	for(int i = lower_edge_smoothed_cutoff; i < upper_edge_cutoff; i ++){
		p_vals[i] = get_mp_p(smoothed_feature_events[i], smoothed_total_events[i]);
	}
	bounds_check(p_vals.size(), upper_edge_cutoff, len - 1);
	for(int i = upper_edge_cutoff; i < len; i ++){
		p_vals[i] = 0;
	}
}

void contig::event_feature::set_totals(event_feature &total){
	if(!valid) return;

	total_events = total.feature_events;
	total_event_count = total.feature_event_count;

	if(total_event_count == 0){
		prop=0;
	}
	else {
		prop = (double)feature_event_count/total_event_count;
	}
}

void contig::event_feature::finish(event_feature &total){
	if(!valid) return;

	set_totals(total);
	smooth();
	get_pvals_poisson();
}

string contig::event_feature::to_str(){
	string ret = tag + " " + cname + "\nlen: " + to_string(len) + "\nevent_count: " + to_string(feature_event_count) + "\nexcluded: " + to_string(excluded) + "\nlower_edge_cutoff: " + to_string(lower_edge_cutoff) + "\nupper_edge_cutoff: " + to_string(upper_edge_cutoff) + "\nflag_cutoff: " + to_string(flag_cutoff) + "\n";
	if(!valid) return ret;
	
	bounds_check(feature_events.size(), 0, len - 1);
	bounds_check(smoothed_feature_events.size(), 0, len - 1);
	bounds_check(total_events.size(), 0, len - 1);
	bounds_check(p_vals.size(), 0, len - 1);

	//for(int i = 0; i < len; i ++) ret += to_string(i) + " " + to_string(feature_events[i]) + " " + to_string(smoothed_feature_events[i]) + " " + to_string(total_events[i]) + " " + to_string(p_vals[i]) + "\n";
	return ret;
}


//Contig functions

contig::contig(double mpevent_logp_flag_cutoff, double doc_logp_flag_cutoff, int mph_smoothing_window_len, int len, string name, int mp_lower_edge_cutoff, int sr_lower_edge_cutoff, int mpdoc_minimum_value, int singleton_minimum_value, int miso_minimum_value, int split_read_minimum_value):
	name(name), 
	len(len), 
	finished(false),
	mp_upper_edge_cutoff(len-mp_lower_edge_cutoff),
	sr_upper_edge_cutoff(len-sr_lower_edge_cutoff),
	mpdoc_minimum_value(mpdoc_minimum_value), 
	singleton_minimum_value(singleton_minimum_value),
	miso_minimum_value(miso_minimum_value),
	split_read_minimum_value(split_read_minimum_value)
{
	valid = mp_lower_edge_cutoff < mp_upper_edge_cutoff && mp_upper_edge_cutoff > 0 && mp_lower_edge_cutoff >= 0 && mp_upper_edge_cutoff <= len && sr_upper_edge_cutoff > sr_lower_edge_cutoff && sr_upper_edge_cutoff > 0 && sr_lower_edge_cutoff >= 0 && sr_upper_edge_cutoff <= len;
	links = vector<mate_pair>();
	split_reads = vector<mate_pair>();

        doc = depth_feature("read_depth_of_coverage", name, len, sr_lower_edge_cutoff, sr_upper_edge_cutoff, doc_logp_flag_cutoff, valid);
        mpdoc = depth_feature("mate_pair_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, doc_logp_flag_cutoff, valid);
        mpdoc_happy = depth_feature("happy_mate_pair_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, doc_logp_flag_cutoff, valid);
        mpdoc_unhappy = depth_feature("unhappy_mate_pair_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff,doc_logp_flag_cutoff, valid, mpdoc_minimum_value);
        mpdoc_short = depth_feature("short_mate_pair_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, doc_logp_flag_cutoff, valid, mpdoc_minimum_value);
        mpdoc_long = depth_feature("long_mate_pair_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, doc_logp_flag_cutoff, valid, mpdoc_minimum_value);
        single_doc = depth_feature("singleton_depth_of_coverage", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, doc_logp_flag_cutoff, valid, singleton_minimum_value);
        miso_doc = depth_feature("misoriented_mate_pair_depth_of_coverage", name, len, 0, len, doc_logp_flag_cutoff, valid, miso_minimum_value);
        split_doc = depth_feature("split_read_depth_of_coverage", name, len, 0, len, doc_logp_flag_cutoff, valid, split_read_minimum_value);
        avg_qual = depth_feature("average_len_mapping_quality", name, len, sr_lower_edge_cutoff, sr_upper_edge_cutoff, -1, valid); //TODO should this even have a cutoff? What does this mean

        mpevents_total = event_feature("all_mate_pair_events", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, mph_smoothing_window_len, mpevent_logp_flag_cutoff, valid);
        mpevents_happy = event_feature("happy_mate_pair_events", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, mph_smoothing_window_len, mpevent_logp_flag_cutoff, valid);
        mpevents_unhappy = event_feature("unhappy_mate_pair_events", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, mph_smoothing_window_len, mpevent_logp_flag_cutoff, valid);
        mpevents_short = event_feature("short_mate_pair_events", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, mph_smoothing_window_len, mpevent_logp_flag_cutoff, valid);
        mpevents_long = event_feature("long_mate_pair_events", name, len, mp_lower_edge_cutoff, mp_upper_edge_cutoff, mph_smoothing_window_len, mpevent_logp_flag_cutoff, valid);
}

contig::contig() :
	len(-1),
	finished(false),
	valid(false),
	mp_upper_edge_cutoff(0),
	sr_upper_edge_cutoff(0),
	mpdoc_minimum_value(0), 
	singleton_minimum_value(0),
	miso_minimum_value(0),

	doc(depth_feature()),
	mpdoc(depth_feature()),
	mpdoc_happy(depth_feature()),
	mpdoc_unhappy(depth_feature()),
	mpdoc_short(depth_feature()),
	mpdoc_long(depth_feature()),
	single_doc(depth_feature()),
	miso_doc(depth_feature()),
	split_doc(depth_feature()),
	avg_qual(depth_feature()), //TODO should this even have a cutoff? What does this mean

	mpevents_total(event_feature()), 
	mpevents_happy(event_feature()),
	mpevents_unhappy(event_feature()),
	mpevents_short(event_feature()),
	mpevents_long(event_feature())
	
{
	name = string("");
	links = vector<mate_pair>();
	split_reads = vector<mate_pair>();
}

void contig::add_split_read(single_read& r){
	if(!valid || !r.valid || finished) return;
	doc.update(r);
	split_doc.update(r);
}
void contig::add_single_read(single_read& r){
	if(!valid || !r.valid || finished) return;
	single_doc.update(r);
	doc.update(r);
	avg_qual.update(r, r.qual);
}

void contig::add_mate_pair(mate_pair& mp){

	if(!valid || !mp.valid || !mp.same_contig || finished) return;
	
	doc.update(mp.a);	
	doc.update(mp.b);	
	avg_qual.update(mp.a, mp.a.qual);
	avg_qual.update(mp.b, mp.b.qual);

	//if this mate pair is misoriented
	if(!mp.orientation){ 
		miso_doc.update(mp.a); 
		miso_doc.update(mp.b);
		return;
	}

	mpdoc.update(mp);
	mpevents_total.update(mp);

	//if this mate pair is happy
	if(mp.happy){
		mpdoc_happy.update(mp);
		mpevents_happy.update(mp);
	}

	//if this mate pair is unhappy
	else{
		mpevents_unhappy.update(mp);
		mpdoc_unhappy.update(mp);

		//if this mate pair is shortened
		if(mp.zscore < 0){
			mpdoc_short.update(mp);
			mpevents_short.update(mp);
		}
		//if this mate pair is stretched
		else{
			mpdoc_long.update(mp);
			mpevents_long.update(mp);
		}
	}
}

void contig::add_links(vector<mate_pair> &contig_links, int scaffold_cutoff){
	if(!valid || contig_links.size() <= 0) return;
	//sort them by their starting position on this contig so that we can easily report the relevant links later.
	sort(contig_links.begin(), contig_links.end(), [](const mate_pair &mp1, const mate_pair &mp2){ return mp1.a.start < mp2.a.start; });
	links = contig_links;

	int len = contig_links.size();
	bounds_check(len, 0, len-1);
	for(int i = 0; i < len; i ++){
		add_single_read(contig_links[i].a);
	} 
}

void contig::add_split_reads(vector<mate_pair> &contig_split_reads){
	int len = contig_split_reads.size();
	if(!valid || len <= 0) return;
	//sort them by their starting position on this contig so that we can easily report the relevant links later.
	sort(contig_split_reads.begin(), contig_split_reads.end(), [](const mate_pair &mp1, const mate_pair &mp2){ return mp1.a.start < mp2.a.start; });
	contig_split_reads = split_reads;

	bounds_check(len, 0, len-1);
	for(int i = 0; i < len; i ++){
		add_split_read(contig_split_reads[i].a);
	} 
}

//should be called when this contig is finished being updated and it's time to get stats
void contig::finish(){
	if(!valid) return;
	finished = true; //todo finish implementing this flag
	if(len <= 0) return;
	doc.finish();
	mpdoc.finish();
	mpdoc_happy.finish();
	mpdoc_unhappy.finish();
	mpdoc_short.finish();
	mpdoc_long.finish();
	miso_doc.finish();
	split_doc.finish();
	single_doc.finish();

	mpevents_happy.finish(mpevents_total);
	mpevents_unhappy.finish(mpevents_total);
	mpevents_long.finish(mpevents_total);
	mpevents_short.finish(mpevents_total);

	//TODO figure out what to do with this. What were we doing with avg_qual??
	avg_qual.cumulative_sum();
	avg_qual.get_bounded_avg_sdv();
	avg_qual.perbase_avg.resize(len, 0.0);  
	
	bounds_check(doc.depth.size(), 0, len - 1);
	bounds_check(avg_qual.perbase_avg.size(), 0, len - 1);
	for(int i = 0; i < len; i++){
		if(doc.depth[i] != 0){ avg_qual.perbase_avg[i] = (double)avg_qual.depth[i]/doc.depth[i]; }
	}
}

string contig::to_str(){
	string ret = name + " " + to_string(len) + "\n";
	if(!valid) return ret;

	ret += doc.to_str() + "\n";
	ret += mpdoc.to_str() + "\n";
	ret += mpdoc_happy.to_str() + "\n";
	ret += mpdoc_unhappy.to_str() + "\n";
	ret += mpdoc_short.to_str() + "\n";
	ret += mpdoc_long.to_str() + "\n";
	ret += miso_doc.to_str() + "\n";
	ret += split_doc.to_str() + "\n";
	ret += single_doc.to_str() + "\n";

	ret += mpevents_happy.to_str() + "\n";
	ret += mpevents_unhappy.to_str() + "\n";
	ret += mpevents_long.to_str() + "\n";
	ret += mpevents_short.to_str() + "\n";
	return ret;
}
