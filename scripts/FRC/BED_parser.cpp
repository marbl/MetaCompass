#include "BED_parser.h"
#include <iostream> //TODO debug only
#include <iterator>
#include <sstream>
#include "helpers.h"

using namespace std;

void BED_parser::reopen(){
	close();
	bed_file.clear();
	if(parsing_contigs) contig_init();
	if(parsing_reads) read_init();
}

void BED_parser::contig_init(){
	bed_file.open(filename, ifstream::in);
	//TODO have this actually throw something useful and output it
	if(!bed_file){
		cout << "Could not open BED file: " << filename << "\n";
		throw 1; 
	}
	parsing_contigs = true;
	check_cnames = filtered_cnames.size() != 0;
	has_next_flag = true; //TODO check for empty file?
	first_read = single_read();
	cur_clen = 0;
	cur_cname = string("");
}

void BED_parser::read_init(){
	bed_file.open(filename, ifstream::in);
	//TODO have this actually throw something useful and output it
	if(!bed_file){
		cout << "Could not open BED file: " << filename << "\n";
		throw 1; 
	}
	parsing_reads = true;
	check_cnames = filtered_cnames.size() != 0;
	has_next_flag = true;
	prev_read = single_read();
	buf_read = single_read();
	buf_pair = mate_pair();
	parsed_single_read = false;
	parsed_mate_pair = false;
	carry_read = false;
}


bool BED_parser::has_next(){ /*cout << "has_next: " << has_next_flag << "\n";*/ return has_next_flag; }

//TODO move this comment
//returns the next contig in the bed file, optionally the next one that's in the set filtered_cnames
//This is useful because the file could be very large and you might not want all the contigs in memory at once. This allows for streaming the input one contig at a time. 
//If memory is an issue, use this function to get one contig at a time, process that contig, then get another one.
//When there are no more contigs in the file, it will return an empty contig, with a name of "" and a length of -1. That's how you know it's done. TODO: make it indicate this in a better way
//TODO what does it do if there are no valid reads in the file?
//TODO this needs to check that the file is sorted correctly instead of just crashing. For example, if we return a contig we should not see reads belonging to that contig later!

vector<string> BED_parser::get_tokenized_line(){
	string line;
	if(getline(bed_file, line)){
		//cout <<"tokenizedline: <" << line << ">\n";
    		istringstream buf(line);
    		istream_iterator<string> beg(buf), end;
		vector<string> tokens(beg, end); 
		return tokens;
	}
	has_next_flag = false;
	return vector<string>();
}
	
single_read BED_parser::get_next_read(){
	if(!parsed_single_read){
		parsed_single_read = true;
		get_next_read();
	}

	string buf;

	single_read next_read = buf_read;
	while(has_next_flag){

		vector<string> data = get_tokenized_line();
		if(!has_next_flag) break;

		//TODO check for bad input

		bounds_check(data.size(), 0, 5);
		string cname = data[0];
		int start = stoi(data[1]);
		int end = stoi(data[2]);
		string name = data[3];
		int qual = stoi(data[4]);
		string direction = data[5];

		string cigar = "";
		if(data.size() >= 7){
			cigar = string(data[6]);
		}

		//this read is from a contig that's not in our list
		if(check_cnames && filtered_cnames.count(cname) == 0){ continue; }
		int clen = -1;
		if(clens.count(cname) != 0) clen = clens[cname];
		
		buf_read = single_read(cname, start, end, name, qual, direction, clen, cigar);   
		return next_read;
	}
	return buf_read;
}

vector<single_read> BED_parser::get_n_reads(int n){
	vector<single_read> reads;
	int ctr = 0;
	while(has_next_flag && ctr < n){
		reads.push_back(get_next_read());
		ctr++;
	}
	return reads;
}

vector<single_read> BED_parser::get_all_reads(){
	vector<single_read> reads;
	while(has_next_flag){
		reads.push_back(get_next_read());
	}
	return reads;
}

mate_pair BED_parser::get_next_mate_pair(bool same_contig){
	if(!parsed_mate_pair){
		parsed_mate_pair = true;
		get_next_mate_pair();
	}
	mate_pair next_pair = buf_pair;

	while(has_next_flag){
		//get the next read
		single_read cur_read = get_next_read();

		//is this the mate of the previous read?
		if(!carry_read){
			prev_read = cur_read;
			carry_read = true;
			continue;
		}
		if(prev_read.is_mate(cur_read)){
			buf_pair = mate_pair(prev_read, cur_read, correct_direction_a, correct_direction_b, libsize_avg, libsize_sdv, mp_size_zscore_cutoff, same_contig);
			carry_read = false;
			return next_pair;
		}
		prev_read = cur_read;
		//skipping over the carried read because it's not part of a pair
	}
	return next_pair;
}

vector<mate_pair> BED_parser::get_all_mate_pairs(bool same_contig){
	vector<mate_pair> mps;
	while(has_next_flag){
		mps.push_back(get_next_mate_pair(same_contig));
	}
	return mps;
}

//TODO what does it do if there are no valid contigs in the file?
contig BED_parser::get_next_contig(){

	//is the current contig one that needs to be returned?
	contig_requested = (!check_cnames || filtered_cnames.count(cur_cname) != 0) && cur_cname != "";
	
	contig cur_contig = contig(mpevent_logp_flag_cutoff, doc_logp_flag_cutoff, mph_smoothing_window_len, cur_clen, cur_cname, mp_lower_edge_cutoff, sr_lower_edge_cutoff, mpdoc_minimum_value, singleton_minimum_value, miso_minimum_value, split_read_minimum_value);
	single_read prev_read = first_read;

	while(has_next_flag){
		single_read cur_read = get_next_read();
		string cname = cur_read.cname;

		//we are starting a new contig, so return the current contig and then erase it and start a new one
		if(cname != cur_cname){
			//cout << "starting new contig: <" << cur_cname << "> closing: <" << cname << "\n";	
			//start a new contig
			cur_cname = cname;
			first_read = cur_read;
			carry_read = true;
			cur_clen = -1;
			if(clens.count(cur_cname) != 0) { cur_clen = clens[cur_cname]; }
			else {}//TODO handle this case!!
			
			//return the previous contig if it was requested
			if(contig_requested){ 
				cur_contig.add_single_read(prev_read); //TODO make sure this won't crash on an empty read or empty contig
				//cur_contig.finish();
				//return t_contig; //todo debug
				//cout << "returning <" << cur_contig.name <<">\n";
				return cur_contig; 
			}
			else{
				contig_requested = (!check_cnames || filtered_cnames.count(cur_cname) != 0) && cur_cname != "";
				cur_contig = contig(mpevent_logp_flag_cutoff, doc_logp_flag_cutoff, mph_smoothing_window_len, cur_clen, cur_cname, mp_lower_edge_cutoff, sr_lower_edge_cutoff);
				//t_contig = contig(); //todo debug
				prev_read = first_read; //todo will this cause double freeing?
				continue;
			}
		}

		//we don't want to bother making a whole object for a contig we don't care about. 
		if(!contig_requested) continue;

		//cout << "attatching read <" << cur_read.name << ", " << cur_read.cname << "> to <" << cur_cname << ">\n"; 	
		//if there is no read to compare this read to, compare the next read to this one next iteration
		if(!carry_read){	
			prev_read = cur_read;
			carry_read = true;
			continue;
		}

		//if this read is the mate of the previous read, add the pair to the contig
		if(prev_read.is_mate(cur_read)){
			carry_read = false;
			mate_pair mp(prev_read, cur_read, correct_direction_a, correct_direction_b, libsize_avg, libsize_sdv, mp_size_zscore_cutoff);
			//cout << "attatching mp <" << mp.to_str() << "> to <" << cur_cname << ">\n";
			cur_contig.add_mate_pair(mp);
			continue;
		}

		//these reads aren't a pair, so we'll compare the current read to the next one next iteration. The previous read is unmated, so add it to the contig as a single read
		cur_contig.add_single_read(prev_read);
		//cout << "attatching single read <" << prev_read.to_str() << "> to <" << cur_cname << ">\n";
		prev_read = cur_read;
	}
	
	//make sure we get the last contig in the file!
	if(contig_requested){ 
		if(carry_read){ cur_contig.add_single_read(prev_read); }
		//cur_contig.finish();
		//return t_contig; //todo debug
		//cout << "returning last contig: <" << cur_contig.name << ">\n";
		return cur_contig;
	}

	//there's no valid contig so return an empty one
	//TODO handle this case!
	return contig();
}

//returns a list of all the contigs in this BED file, optionally just those in the set filtered_cnames
//Be careful! This will read the whole BED file into memory at once, which could be very large. If memory is a problem, use get_next_contig for streaming input. Get one contig, process that contig, then get the next one.
vector<contig> BED_parser::get_all_contigs(){ 
	vector<contig> contigs;
	while(has_next_flag){ //TODO make SURE this won't ever get stuck
		contigs.push_back(get_next_contig());
	}
	return contigs;
}

void BED_parser::set_filter(unordered_set<string> filter_cnames){
	filtered_cnames = filter_cnames;
}

//TODO get the cutoffs out of the constructor or make a constructor that doesn't require them!
BED_parser::BED_parser(string filename, map<string, int> clens, string correct_direction_a, string correct_direction_b, double mpevent_logp_flag_cutoff, double doc_logp_flag_cutoff, double mp_lower_edge_cutoff, double sr_lower_edge_cutoff, double mp_size_zscore_cutoff, double libsize_avg, double libsize_sdv, int mpdoc_minimum_value, int singleton_minimum_value, int miso_minimum_value, int split_read_minimum_value, unordered_set<string> filtered_cnames) : 
	filename(filename), 
	clens(clens), 
	correct_direction_a(correct_direction_a),
	correct_direction_b(correct_direction_b),
	mpevent_logp_flag_cutoff(mpevent_logp_flag_cutoff), 
	doc_logp_flag_cutoff(doc_logp_flag_cutoff), 
	mp_lower_edge_cutoff(mp_lower_edge_cutoff), 
	sr_lower_edge_cutoff(sr_lower_edge_cutoff),
	mp_size_zscore_cutoff(mp_size_zscore_cutoff), 
	libsize_avg(libsize_avg), 
	libsize_sdv(libsize_sdv),
	mpdoc_minimum_value(mpdoc_minimum_value),
	singleton_minimum_value(singleton_minimum_value),
	miso_minimum_value(miso_minimum_value),
	split_read_minimum_value(split_read_minimum_value),
	filtered_cnames(filtered_cnames)
{ 
	parsing_reads = false;
	parsing_contigs = false;
	contig_init(); 
}


BED_parser::BED_parser(string filename, string correct_direction_a, string correct_direction_b, double mp_size_zscore_cutoff, double libsize_avg, double libsize_sdv, unordered_set<string> filtered_cnames):
	filename(filename),
	correct_direction_a(correct_direction_a),
	correct_direction_b(correct_direction_b),
	mp_size_zscore_cutoff(mp_size_zscore_cutoff),
	libsize_avg(libsize_avg),
	libsize_sdv(libsize_sdv),
	filtered_cnames(filtered_cnames)
{ 
	parsing_reads = false;
	parsing_contigs = false;
	read_init(); 
}

void BED_parser::close(){
	bed_file.close();
}
