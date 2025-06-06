#ifndef SINGLE_READ_H
#define SINGLE_READ_H

#include <string>
#include "vector"

class single_read{
	public:
	std::string name; //read name
	std::string cname; //contig name that it maps to
	std::string direction; //mapping orientation
	std::string cigar; //cigar string
	std::vector<char> cigar_labels; //for storing tokenized cigar
	std::vector<int> cigar_lens; //for storing tokenized cigar
	std::string seq; //sequence information
	std::string fastq_spacer; //the spacer for the fastq file
	std::string pb_qual; //the per-base quality string
	int clen; //length of the contig it maps to
	int start; //start position of mapping
	int end; //end position of mapping
	int len; //read len
	int qual; //mapping quality
	int score; //alignment score. Right now this is read in as qual because bed files can only have either qual or score. Assign it to this field manually. //TODO find a better way to handle this
	bool valid; //is this read valid?

	single_read(std::string cname, int start1, int end1, std::string name, int qual, std::string direction, int clen, std::string cigar = std::string("")); //constructor for BED file
	single_read(std::string name, std::string seq, std::string fastq_spacer, std::string pb_qual); //constructor for FASTQ file
	single_read(single_read r, std::string name, std::string cigar, std::string seq, std::string pb_qual, int start, int end, int len); //constructor for splitting reads
	single_read(); //default constructor (invalid read)
	bool operator < (const single_read & r) const;	
	
	std::vector<single_read> split(int breakpoint, int cigar_pos = -1);	//Splits a single read at the specified position in the read. Returns a vector containing two reads. The breakpoint is exclusive on the first read (vector[0]) and inclusive on the second one (vector[1])
	bool cigar_label_is_match(char label);
	bool cigar_label_is_consume(char label);
	std::vector<single_read> cigar_split();	 //finds the "best" place to split a read based on the cigar string. See the documentation for what this is for
	bool check_edge(int lo, int hi);
	bool is_mate(single_read r);
	bool check_direction();
	std::string get_umph(double avg, double sdv, double cutoff);
	std::string to_str();
	std::string to_bed_str();
	std::string to_fq_str();
	void tokenize_cigar();
};

#endif
