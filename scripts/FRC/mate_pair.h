#ifndef MATE_PAIR_H
#define MATE_PAIR_H

#include <string>
#include "single_read.h"

class mate_pair{
	public:
		single_read a;
		single_read b;
		std::string correct_direction_a;
		std::string correct_direction_b;

		bool same_contig;
		bool orientation;
        int dist;
		bool happy;
		double zscore;
		bool valid;

        mate_pair (single_read a, single_read b, std::string correct_direction_a, std::string correct_direction_b, double avg, double sdv, double cutoff, bool same_contig = true);
		mate_pair (single_read a, single_read b, bool same_contig = false); //to be used if they are not on the same contig
		mate_pair();
		void orient();
        bool check_orientation();
		int get_dist();
		bool check_edge(int lo, int hi);
		void reverse();
		std::string to_str();
		std::string to_bed_str();
};

#endif
