#include <cmath>
#include <iostream> //TODO debug only, remove
#include "mate_pair.h"

using namespace std;

mate_pair::mate_pair(){
	a = single_read();
	b = single_read();
	correct_direction_a = string("");
	correct_direction_b = string("");
	same_contig = false;
	orientation = false;
	dist = -1;
	happy = false;
	zscore = 0;
	valid = false;
}

mate_pair::mate_pair(single_read a, single_read b, bool same_contig) : 
	a(a), 
	b(b), 
	same_contig(same_contig),
	valid(true)
{}

mate_pair::mate_pair (single_read a, single_read b, string correct_direction_a, string correct_direction_b, double avg, double sdv, double cutoff, bool same_contig):
	a(a), 
	b(b), 
	correct_direction_a(correct_direction_a), 
	correct_direction_b(correct_direction_b), 
	same_contig(same_contig),
	dist(-1), 
	happy(false), 
	zscore(0),
	valid(true) 
	{
		if(!same_contig){ same_contig = a.cname == b.cname; } //this is really slow, make sure to mark it false only when needed
		if(same_contig){
			orient();
			orientation = check_orientation();

			//TODO: epsilon checking
			if(orientation){
					dist = get_dist();
					zscore = (dist - avg)/sdv;
					happy = abs(zscore) < cutoff;
			}
		}
	}

//TODO do this with pointers instead of copying every time
//put these reads in order, a is + and b is -
void mate_pair::reverse(){
		if(!valid) return;
		single_read swap(a);
		a = single_read(b);
		b = single_read(swap);
}
void mate_pair::orient(){
	if(!valid) return;
	if(!same_contig || correct_direction_a == correct_direction_b){ return; }
	if(a.direction == correct_direction_b && b.direction == correct_direction_a){ reverse(); }
}

bool mate_pair::check_orientation(){
	//cout << "Orientation check for " << a.name << " " << b.name << ": <" <<a.direction << "> <" << b.direction << "> <" << correct_direction_a << "> <" << correct_direction_b << "> <" << a.start << "> <" << b.start << "> <" << a.end << "> <" << b.end << "> <" << (a.direction == correct_direction_a && b.direction == correct_direction_b && a.start < b.start && a.end < b.end) << ">\n"; 
	if(!valid) return false;
	return same_contig && a.direction == correct_direction_a && b.direction == correct_direction_b && a.start <= b.start && a.end <= b.end;
}

int mate_pair::get_dist(){
		if(!same_contig || !valid){ return -1; }
		return max(a.end, b.end) - min(a.start, b.start);
}

//check that the mate pair is not too close to the end of the contig, assumes they are in order
bool mate_pair::check_edge(int lo, int hi){
		if(!valid) return false;
		return same_contig && !(a.start < lo || b.end > hi);
}

string mate_pair::to_str(){
		//flag is 1 if this read is happy
		//1 for true and 0 for false
		//apparently c++ has some weird rules for concatenating char*, hence the string("/t")
		return (happy ? "1" : "0") + string("\t") + a.to_str() + "\t" + b.to_str() + "\t" + (orientation ? "1" : "0") + "\t" + to_string(dist) + "\t" + to_string(zscore);
}

string mate_pair::to_bed_str(){
	return a.to_bed_str() + "\n" + b.to_bed_str();
}	
