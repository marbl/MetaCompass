#ifndef BED_PARSER_H
#define BED_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <map>
#include <vector>


#include "single_read.h"
#include "mate_pair.h"
#include "contig.h"

//TODO throw if the bed file isn't correctly sorted!
//TODO handle qual of *
class BED_parser{

	std::ifstream bed_file;
	std::string filename;
	std::map<std::string, int> clens;
	bool check_cnames;

	std::string correct_direction_a;
	std::string correct_direction_b;
	double mpevent_logp_flag_cutoff;
	double doc_logp_flag_cutoff;
	double mp_lower_edge_cutoff; 
	double sr_lower_edge_cutoff;
	double mp_size_zscore_cutoff;
	double libsize_avg;
	double libsize_sdv;
	int mpdoc_minimum_value;
	int singleton_minimum_value;
	int miso_minimum_value;
	int split_read_minimum_value;
	std::unordered_set<std::string> filtered_cnames;

	bool carry_read;
	bool parsed_single_read;
	bool parsed_mate_pair;
	bool parsing_contigs;
	bool parsing_mps;
	bool parsing_reads;
	single_read buf_read;
	mate_pair buf_pair;
	bool contig_requested;
	bool has_next_flag;

	single_read prev_read;
	std::string cur_cname;
	single_read first_read;
	int cur_clen;

	//TODO move this somewhere else!
	//info needed for contig objects
	int mph_smoothing_window_len;

	void contig_init();
	void read_init();

	std::vector<std::string> get_tokenized_line();

public:

	void reopen();
	bool has_next();

	void set_filter(std::unordered_set<std::string> filter_cnames);
	
	
	//Get the next read as a read object
	//DO NOT use this with other get_next functions, it will mess up the reading frame of the file. Only use one kind of get_next function per BED_parser!
	single_read get_next_read();

	std::vector<single_read> get_n_reads(int n);
	//returns a vector of all the reads as single read objects
	//Be careful! This will read the whole file at once, which could be very large. If memory is a problem, use get_next_read for streaming input! Get one read, process that one, then get the next one.
	//DO NOT use this with other get_next functions, it will mess up the reading frame of the file. Only use one kind of get_next function per BED_parser!
	std::vector<single_read> get_all_reads();

	//get the next mate pair. It will skip over singletons. This assumes that mates will be next to each other in the file (interleaved), so the results depend on how the file is sorted.
	//For example, if you sort only by read name, it will return the next pair regardless of there the reads map. But if you sort first by contig and then by read name,
	//it will return the next mate pair where both mates are on the same contig.
	//DO NOT use this with other get_next functions, it will mess up the reading frame of the file. Only use one kind of get_next function per BED_parser!
	mate_pair get_next_mate_pair(bool same_contig = true);

	//get all the mate pairs. Warning: This could be very large! You can use get_next_mate_pair for streaming input
	std::vector<mate_pair> get_all_mate_pairs(bool same_contig = true);

	//returns the next contig in the bed file, optionally the next one that's in the set filtered_cnames
	//This is useful because the file could be very large and you might not want all the contigs in memory at once. This allows for streaming the input one contig at a time. 
	//If memory is an issue, use this function to get one contig at a time, process that contig, then get another one.
	//When there are no more contigs in the file, it will return an empty contig, with a name of "" and a length of -1. That's how you know it's done. TODO: make it indicate this in a better way
	//DO NOT use this with other get_next functions, it will mess up the reading frame of the file. Only use one kind of get_next function per BED_parser!
	contig get_next_contig();

	//returns a list of all the contigs in this BED file, optionally just those in the set filtered_cnames
	//Be careful! This will read the whole BED file into memory at once, which could be very large. If memory is a problem, use get_next_contig for streaming input. Get one contig, process that contig, then get the next one.
	//DO NOT use this with other get_next functions, it will mess up the reading frame of the file. Only use one kind of get_next function per BED_parser!
	std::vector<contig> get_all_contigs();

	//these are to be used when you want the parser to return contig objects. It needs all these cutoffs to compute the contig information
	//TODO Make them have defaults of no cutoff for edges and something common for the mp lens!
	//TODO this constructor is pretty specific to VALET. Make a generic one for others to use!
	BED_parser(std::string filename, std::map<std::string, int> clens, std::string correct_direction_a, std::string correct_direction_b, double mpevent_logp_flag_cutoff, double doc_logp_flag_cutoff, double mp_lower_edge_cutoff, double sr_lower_edge_cutoff, double mp_size_zscore_cutoff, double libsize_avg, double libsize_sdv, int mpdoc_minimum_value = 0, int singleton_minimum_value = 0, int miso_minimum_value = 0, int split_read_minimum_value = 0, std::unordered_set<std::string> filtered_cnames = std::unordered_set<std::string>());

	//these are only to be used when getting reads or mate pairs. If you want information on the mate-pair happiness (like if it's too long or short) then you can provide the mp distribution info. 
	//You can't make contig objects without all the other cutoff information yet!
	BED_parser(std::string filename, std::string correct_direction_a = "", std::string correct_direction_b = "", double mp_size_zscore_cutoff = 0.0, double libsize_avg = 0, double libsize_sdv = 0, std::unordered_set<std::string> filtered_cnames = std::unordered_set<std::string>());
	BED_parser(std::string filename, double score_cutoff); //constructor for reading in split reads from bed file
	void close();
};

#endif