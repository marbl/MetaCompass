#ifndef CONTIG_H
#define CONTIG_H

#include <string>
#include <vector>
#include "mate_pair.h"

class contig {
public:
	
	class depth_feature {
	public:

		std::string tag; //the name of this feature. Used for debugging and reports
		std::string cname; //the name of the contig
		int len; //the length of the contig
		int read_count; //the total number of reads that counted towards this feature
		int excluded; //the number of reads that were excluded because of edge cutoffs
		int lower_edge_cutoff; //the index on the contig below which we will not count this feature
		int upper_edge_cutoff; //the index on the contig above which we will not count this feature
		int minimum_value; //the minimum value to be considered significant
		int bounded_len; //the length between the cutoffs
		double bounded_avg; //the average value of this feature between the cutoffs
		double bounded_sdv; //the standard deviation of that average
		double flag_cutoff; //the p value or other cutoff that will be used to determine when the feature is significant
		
		std::vector<int> depth; //the per-base depth of this feature
		std::vector<double> p_vals; //the per-base p value of observing this depth
		std::vector<double> perbase_avg; //optional per-base average of a feature

		bool valid;

		depth_feature(std::string tag, std::string cname, int len, int lower_edge_cutoff, int upper_edge_cutoff, double flag_cutoff, bool contig_valid, int minimum_value = 0);
		depth_feature();

		void update(mate_pair& mp, int value = 1);
		void update(single_read& r, int value = 1);
		void cumulative_sum();
		void get_bounded_avg_sdv();

		//get the probability of the depth of coverage at a base using z scores, given the average and stdev of the doc on the contig
		//this will flag too low or too high doc because, though the test is one sided, it uses the absolute value of the zscore
		double get_doc_p(double avg_doc, double sdv_doc, int doc);
		void get_pvals_zscore();
		void finish();
		std::string to_str();
	};

	class event_feature {
	public:
		std::string tag; //the name of the feature
		std::string cname; //the name of the contig
		int len; //the length of the contig
		std::vector<int> feature_events; //the events of this feature type per base on the contig
		std::vector<int> smoothed_feature_events; //the number of events in the window of length smoothing_window_len starting at each index
		std::vector<int> total_events; //the total number of events of all feature types that we are comparing this feature with, including this one, per base
		std::vector<int> smoothed_total_events; //the number of total events in the wondow of length smoothing_window_len starting at each index
		int feature_event_count; //the number of events of this feature that are within the bounds of the edge cutoffs
		int total_event_count; //the number of total events that are within the bounds of the edge cutoffs
		int excluded; //the number of events that are excluded based on the edge cutoffs 
		int lower_edge_cutoff; //the index on the contig below which we will not count this feature
		int upper_edge_cutoff; //the index on the contig above which we will not count this feature
		int minimum_value; //the minimum value needed to be considered significant
		int smoothing_window_len; //the length of the window that will be used to group events for poisson calculations
		int lower_edge_smoothed_cutoff; //the index on the contig below which we will not count this feature, accounting for smoothing
		double flag_cutoff; //the p value above which we will consider this feature significant
		std::vector<double> p_vals; //the per window p vals of the windows starting at each index
		double prop; //the proportion of events of this type to the total number of comparable events over the whole contig

		bool valid;

		event_feature(std::string tag, std::string cname, int len, int lower_edge_cutoff, int upper_edge_cutoff, int smoothing_window_len, double flag_cutoff, bool contig_valid, int minimum_value = 0);
		event_feature();
		void update(mate_pair& mp, int value = 1);
		void update(single_read& r, int lower_edge_cutoff, int upper_edge_cutoff);
		void smooth();
		//get the probability of observing this many specific mate pairs out of the total number, given the overall proporion of that type to the total
		//This won't flag having too few unhappy mps because the test is one sided, we are getting the probability of seeing x or more events
		double get_mp_p(int observed, int total);
		void get_pvals_poisson();
		void set_totals(event_feature &total);
		void finish(event_feature &total);
		std::string to_str();
	};

	std::string name; //the contig name from the FASTA file
	int len; //the contig length
	bool finished; //keeps track of whether or not this contig is finished. It can only be updated if not finished and stats can only be calculated if it is not.
	bool valid; //is this a valid contig? If the edge cutoffs overlap, it isn't TODO remove?
	int mp_upper_edge_cutoff; //TODO make the constructor not an initializer list so these can be locally created in the constructor. They don't need to be stored here.
	int sr_upper_edge_cutoff;
	int mpdoc_minimum_value; //minimum number of unhappy/shor/long mps needed to flag a region
	int singleton_minimum_value; //minimum number of singletons needed to flag a region
	int miso_minimum_value; //minimum number of misoriented pairs needed to flag a region
	int split_read_minimum_value;

	depth_feature doc; //depth of coverage, ignoring mate pair information
	depth_feature mpdoc; //mate-pair depth of coverage, including the insert
	depth_feature mpdoc_happy; //depth of coverage of concordant mate pairs
	depth_feature mpdoc_unhappy; //depth of coverage of discordant mate pairs
	depth_feature mpdoc_short; //depth of coverage of short mate pairs
	depth_feature mpdoc_long; //depth of coverage of long mate pairs
	depth_feature single_doc; //deoth of coverage of reads with no mate on the same contig
	depth_feature miso_doc; //depth of coverage of reads that are part of misoriented mate pairs
	depth_feature split_doc; //doc of split reads
    	depth_feature avg_qual; //used only to keep track of the per-base average mapping quality. This is not used for flagging!
	
	event_feature mpevents_total; //Represents the "events" at each base, where an event is the start or end of a mate pair. This is used only to keep track of the total number of mate pair events and isn't used for flagging. If you call smooth or flag on this feature it will crash!
	event_feature mpevents_happy; //Happy mate pair events at each base
	event_feature mpevents_unhappy; //Unhappy mate pair events at each base
	event_feature mpevents_short; //Short mate pair events at each base
	event_feature mpevents_long; //Long mate pair events at each base 

	std::vector<mate_pair> links; // mate pairs that link this contig to other contigs and that are not scaffold links
	std::vector<mate_pair> split_reads; //split reads where at least one part of the read is on this contig
	contig(double mpevent_logp_flag_cutoff, double doc_logp_flag_cutoff, int mph_smoothing_window_len, int len = 0, std::string name = std::string(""), int mp_lower_edge_cutoff = -1, int sr_lower_edge_cutoff = -1, int mpdoc_minimum_value = 0, int singleton_minimum_value = 0, int miso_minimum_value = 0, int split_read_minimum_value = 0);
	contig();

	void add_single_read(single_read& r);
	void add_split_read(single_read& r);
	void add_mate_pair(mate_pair& mp);
	void add_links(std::vector<mate_pair> &contig_links, int scaffold_cutoff = 0);
	void add_split_reads(std::vector<mate_pair> &split_reads);
	void finish(); 	//should be called when this contig is finished being updated and it's time to get stats
	std::string to_str();

};

#endif
