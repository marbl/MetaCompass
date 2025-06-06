#include <cmath>
#include "single_read.h"
#include "helpers.h"
#include <iostream> //For cout TODO remove

using namespace std;

void single_read::tokenize_cigar(){
	if(cigar == "" || cigar.length() <= 0 || cigar_lens.size() != 0) return;
	int len = cigar.length();
	bounds_check(len, 0, len-1);
	string cur_num = "";
	for(int i = 0; i < len; i ++){
		//we're still gathering digits
		if(isdigit(cigar[i])) cur_num += cigar[i];
		//parse this group
		else if (isalpha(cigar[i]) || cigar[i] == '='){
			cigar_lens.push_back(stoi(cur_num));
			cur_num = "";
			cigar_labels.push_back(cigar[i]);
		}
		else{
			cout << "invalid character in CIGAR string: " << cigar[i] << " in read " << name <<"\n";
			cigar_lens = vector<int>(0);
			cigar_labels = vector<char>(0);
			return;
		}
	}
}

//note: this is not efficient for recursive splitting because the cigar strings have to be parsed again
//If recursive splitting is needed, add code to pass down the parsed cigar vectors
vector<single_read> single_read::split(int breakpoint, int cigar_breakpoint){
	vector<single_read> ret;
	if(!valid || breakpoint < 0) return ret;
	breakpoint = min(breakpoint, len);

	//get the names of the new reads
	string left_name = name + "_1";
	string right_name = name + "_2";

	//get the cigar strings if needed
	string left_cigar = "";
	string right_cigar = "";
	int cigar_len = cigar.length();
	cigar_breakpoint = min(cigar_len, cigar_breakpoint);
	if(cigar != "" && cigar_len > 0){
		if(cigar_breakpoint == -1){
			int consume_left = 0;
			int consume_right = 0;
			int total_left = 0;
			int total_right = 0;
			assert(cigar_labels.size() == cigar_lens.size()); //TODO make this more useful, fail gracefully
			int n_cigars = cigar_labels.size();
			cigar_breakpoint = n_cigars - 1;
			bounds_check(n_cigars, 0, n_cigars - 1);
			bounds_check(cigar_lens.size(), 0, n_cigars - 1);
			for(int i = 0; i < n_cigars; i ++){
				if(cigar_label_is_consume(cigar_labels[i])){
					consume_right += cigar_lens[i];
				}
				total_right += cigar_lens[i];
			}
			for(int i = 0; i < n_cigars; i ++){
				if(cigar_label_is_consume(cigar_labels[i])){
					consume_left += cigar_lens[i];
					consume_right -= cigar_lens[i];
				}
				total_left += cigar_lens[i];
				total_right -= cigar_lens[i];

				if(consume_left >= breakpoint){
					cigar_breakpoint = i;
					break;
				}
			}

			bounds_check(n_cigars, 0, cigar_breakpoint - 1);
			bounds_check(n_cigars, cigar_breakpoint, n_cigars - 1);
			for(int i = 0; i < cigar_breakpoint; i ++){
				left_cigar += to_string(cigar_lens[i]);
				left_cigar += cigar_labels[i];
			}

			if(consume_left == breakpoint){
				left_cigar += to_string(cigar_lens[cigar_breakpoint]);
				left_cigar += cigar_labels[cigar_breakpoint];
			}	

			else{
				int difference = consume_left - breakpoint;
				left_cigar += to_string(cigar_lens[cigar_breakpoint] - difference);
				left_cigar += cigar_labels[cigar_breakpoint];
				right_cigar += to_string(difference);
				right_cigar += cigar_labels[cigar_breakpoint];
			}	
			for(int i = cigar_breakpoint + 1; i < n_cigars; i ++){
				right_cigar += to_string(cigar_lens[i]);
				right_cigar += cigar_labels[i];
			}
		}

		else{
			left_cigar = cigar.substr(0, cigar_breakpoint);
			right_cigar = cigar.substr(cigar_breakpoint, cigar.length());
		}
	}

	//split the sequence information if needed
	string left_seq = "";
	string right_seq = "";
	int seq_len = seq.length();
	int seq_break = min(breakpoint, seq_len);
	if(seq != "" && seq.length() > 0){
		left_seq = seq.substr(0, seq_break);
		right_seq = seq.substr(seq_break, seq_len);
	}

	//split the quality string if needed
	string left_pbqual = "";
	string right_pbqual = "";
	int pbqual_len = pb_qual.length();
	int pbqual_break = min(pbqual_len, breakpoint);
	if(pb_qual != "" && pbqual_len > 0){
		left_pbqual = seq.substr(0, pbqual_break);
		right_pbqual = seq.substr(pbqual_break, pbqual_len);
	}

	//TODO get rid of *this, it's just asking for trouble
	single_read left(*this, left_name, left_cigar, left_seq, left_pbqual, 0, breakpoint, breakpoint);
	single_read right(*this, right_name, right_cigar, right_seq, right_pbqual, breakpoint, len, len - breakpoint);
	ret.push_back(left);
	ret.push_back(right);

	return ret;
}

bool single_read::operator < (const single_read &r) const
{
        return (name.compare(r.name) < 0);
}

bool single_read::cigar_label_is_match(char label){
	return label == 'M' || label == '=';
}

bool single_read::cigar_label_is_consume(char label){
	return label == 'M' || label == '=' || label == 'X' || label == 'I' || label == 'S'; 
}

vector<single_read> single_read::cigar_split(){
	vector<single_read> ret;
	if((cigar == "" && cigar_lens.size() == 0) || (cigar_lens.size() != cigar_labels.size()) || !valid) return ret;

	//we want to put a cutoff for the number of matches so that we can differentiate the case where the read doesn't really match at all from the 
	//case where one side matches but the other side doesn't. I chose 50%
	double minimum_match_ratio = 0.5; //TODO move this
	int min_split_size = 0.33*(double)len; //TODO move this
	tokenize_cigar();

	//find the position in the cigar that maximizes the sum of the ratios of matched bases to all bases on one side and unmatched bases to all bases on the other side.
	assert(cigar_labels.size() == cigar_lens.size()); //TODO make this more useful, fail gracefully
	int n_cigars = cigar_labels.size();
	bounds_check(n_cigars, 0, n_cigars - 1);
	bounds_check(cigar_lens.size(), 0, n_cigars - 1);
	int match_left = 0;
	int match_right = 0;
	int mismatch_left = 0;
	int mismatch_right = 0;
	int consumed = 0;
	int total_right = 0;
	double max_sum = 0;
	int max_sum_cigar_index = -1;
	int max_sum_read_index = -1;

	for(int i = 0; i < n_cigars; i ++){
		if(cigar_label_is_consume(cigar_labels[i])){
			consumed += cigar_lens[i];
		}
		if (cigar_label_is_match(cigar_labels[i])){
			match_right += cigar_lens[i];
		}
		else{
			mismatch_right += cigar_lens[i];
		}
		total_right += cigar_lens[i];
	}

	int total_match = match_right;
	int total_mismatch = mismatch_right;

	for(int i = 0; i < n_cigars; i ++){
		if (cigar_label_is_match(cigar_labels[i])){
			match_left += cigar_labels[i];
			match_right -= cigar_labels[i];
		}
		else{
			mismatch_left += cigar_labels[i];
			mismatch_right -= cigar_labels[i];
		}

		if(cigar_label_is_consume(cigar_labels[i])){
			consumed += cigar_lens[i];
		}

		if(consumed < min_split_size) continue;
		if(len - consumed < min_split_size) break;

		double match_ratio_l = (double)match_left/(total_match);
		double mismatch_ratio_l = (double)mismatch_left/(total_mismatch);
		double match_ratio_r = (double)match_right/(total_match);
		double mismatch_ratio_r = (double)mismatch_right/(total_mismatch);

		if(match_ratio_l >= minimum_match_ratio && match_ratio_l + mismatch_ratio_r > max_sum){
			max_sum = match_ratio_l + mismatch_ratio_r;
			max_sum_cigar_index = i;
			max_sum_read_index = consumed;
			max_sum_cigar_index = i;
		}

		if(match_ratio_r >= minimum_match_ratio && match_ratio_r + mismatch_ratio_l > max_sum){
			max_sum = match_ratio_r + mismatch_ratio_l;
			max_sum_cigar_index = i;
			max_sum_read_index = consumed;
			max_sum_cigar_index = i;
		}
	}

	//we found a breakpoint that looks promising!
	if(max_sum_read_index != -1 && max_sum_cigar_index != -1){
		ret = split(max_sum_read_index, max_sum_cigar_index);
		return ret;
	}

	//we didn't find a breakpoint that looks promising. Maybe this read is mapped to the wrong location entirely, or maybe it's very poor quality, etc
	//in this case we'll split the first and last thirds.
	vector<single_read> left = split(len/3);
	ret.push_back(left[0]);
	vector<single_read> right = split(2*(len/3));
	ret.push_back(right[1]);
	return ret;
}

bool single_read::check_edge(int lo, int hi){ return valid && !(start < lo || end > hi); }

bool single_read::is_mate(single_read r){
	if(!valid || !r.valid || name.length() == 0 || r.name.length() == 0) return false;
	string name1 = name.substr(0, name.length() - 1);
	string name2 = r.name.substr(0, r.name.length() - 1);
	string m1 = name.substr(name.length()-1);
	string m2 = r.name.substr(r.name.length()-1);
	return  name1 == name2 && ((m1 =="1" && m2 == "2") || (m2 == "1" && m1 == "2"));
}

bool single_read::check_direction(){
	return valid && (start < end);
	//return (direction == "+" && start < end) || (direction == "-" && end < start);
}

//get output info on a read whose mate is not on the same contig
string single_read::get_umph(double avg, double sdv, double cutoff){
	if(!valid) return "";
	//TODO:assumes that this read could be either end, verify that this is correct
	//TODO instead of recalculating the z score every time, just get what the cutoff values are
	double dist = direction == "+" ? end : clen - start;
	double zscore = (dist - avg)/sdv;
	//flag is 1 if this read is happy
	string flag = (abs(zscore) < cutoff) ? "0" : "1";
	return flag + "\t" + to_str() + "\t" + to_string(zscore);
}

string single_read::to_str(){
	string ret = "read_name: " + name + " len: " + to_string(len);
	if(qual != -1){
		ret += " qual: " + to_string(qual);
	}
	if(score != -1){
		ret += " score: " + to_string(score);
	}
	if(cname != ""){
		ret += " cname: " + cname;
	}
	if(clen != -1){
		ret += " clen: " + to_string(clen);
	}
	if(start != -1){
		ret += " start: " + to_string(start);
	}
	if(end != -1){
		ret += " end: " + to_string(end);
	}
	if(direction != ""){
		ret += " direction: " + direction;
	}
	if(cigar != ""){
		ret += " cigar: " + cigar;
	}
	if(seq != ""){
		ret += " seq: " + seq;
	}
	if(fastq_spacer != ""){
		ret += " fastq_spacer: " + fastq_spacer;
	}
	if(pb_qual != ""){
		ret += " pb_qual: " + pb_qual;
	}

	return ret;
}

string single_read::to_bed_str(){
	string ret = cname + " " + to_string(start) + " " + to_string(end) + " " + name;
	if(qual != -1){
		ret += " " + to_string(qual);
	}
	else if(score != -1){
		ret += " " + to_string(score);
	}
	ret += " " + direction;
	if(cigar != "") ret += " " + cigar;
	return ret;
}

string single_read::to_fq_str(){
	return "@" + name + "\n" + seq + "\n" + fastq_spacer + "\n" + pb_qual;
}

single_read::single_read(string cname, int start1, int end1, string name, int qual, string direction, int clen, string cigar):
	name(name), 
	cname(cname), 
	direction(direction), 
	cigar(cigar),
	clen(clen),
	start(min(start1, end1)), 
	end(max(start1, end1)), 
	len(end - start), 
	qual(qual),
	valid(true)
{
	score = -1;
	seq = string("");
	fastq_spacer = string("");
	cigar_labels = vector<char>();
	cigar_lens = vector<int>();
	pb_qual = string("");
}

single_read::single_read(string name, string seq, string fastq_spacer, string pb_qual) :
	name(name), 
	seq(seq), 
	fastq_spacer(fastq_spacer), 
	pb_qual(pb_qual),
	valid(true)
{
	len = seq.length();
	cname = string("");
	direction = string("");
	cigar = string("");
	cigar_labels = vector<char>();
	cigar_lens = vector<int>();
	clen = -1;
	start = -1;
	end = -1;
	qual = -1;
	score = -1;	
}

//constructor for splitting a read at a breakpoint
single_read::single_read(single_read r, std::string name, std::string cigar, std::string seq, std::string pb_qual, int start, int end, int len) :
	name(name),
	cigar(cigar),
	seq(seq),
	pb_qual(pb_qual),
	start(start),
	end(end),
	len(len)
{
	cname = r.cname;
	clen = r.clen;
	direction = r.direction;
	cigar_labels = vector<char>();
	cigar_lens = vector<int>();
	fastq_spacer = "+"; //the read name has changed so the spacer should change, but I'm not going to case it out
	qual = -1; //mapping quality is no longer relevant when a read is being split
	score = -1;
	valid = r.valid;
}

single_read::single_read()
{
	name = string("");
	cname = string("");
	direction = string("");
	cigar = string("");
	cigar_labels = vector<char>();
	cigar_lens = vector<int>();
	seq = string("");
	fastq_spacer = string("");
	pb_qual = string("");
	clen = -1;
	start = -1;
	end = -1;
	len = -1;
	qual = -1;
	score = -1;
	valid = false;
}
