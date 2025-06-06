#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <exception>
#include <assert.h>
#include <unordered_map>

#include "helpers.h"
#include "contig.h"
#include "flag_window.h"
#include "flagged_region.h"
#include "FAI_parser.h"
#include "BED_parser.h"
#include "feature_calculator.h"

using namespace std;

int main(int argc, char** argv){

	//TODO input checking
	int arg = 1;

	for(int i = 0; i < argc; i ++) cout << "(" << i << ", " << argv[i] << "), ";
	cout << "\n";

	//Files:
	cout << "bed <" << argv[arg] << ">\n";
	string bed_filename = string(argv[arg]);
	arg++;

	cout << "fai <" << argv[arg] << ">\n";
	string fai_filename = string(argv[arg]);
	arg++;

	cout << "out <" << argv[arg] << ">\n";
	FILE* feature_file = fopen(argv[arg], "w");
	if(!feature_file){
		cout << "could not open feature output file: " << argv[arg] << "\n";
		return 0;
	}
	arg++;

	//Data properties:
	cout << "lib_avg <" << argv[arg] << ">\n";	
	double libsize_avg = stod(argv[arg]); //The average library size, accounting for edge effects. Calculated in estimate_libsize.cpp
	arg++;
	
	cout << "lib_sdv <" << argv[arg] << ">\n";
	double libsize_sdv = stod(argv[arg]); //The standard deviation of the library size, accounting for edge effects. Calculated in estimate_libsize.cpp
	arg++;
	
	cout << "dir_a <" << argv[arg] << ">\n";
	string correct_direction_a = string(argv[arg]); //The orientation that we expect the first mate in a pair to be. "+" or "-"
	arg++;

	cout << "dir_b <" << argv[arg] << ">\n";
	string correct_direction_b = string(argv[arg]); //The orientation that we expect the second matr in a pair to be. "+" or "-"
	arg++;

	//Statistical Parameters
	cout << "mps_z<" << argv[arg] << "!!!!!!!!!!!!!!!!!!!!!!!!!!>\n";
	double mp_size_zscore_cutoff = stod(argv[arg]); //The zscore of mate-pair length above which we will call a mate pair too long and below which we will call it too short
	arg++;
	
	cout << "<" << argv[arg] << ">\n";
	double mp_edge_z_cutoff = stod(argv[arg]); //We will exclude mate pairs that are located within avg + z * stdev basepairs of contig edges 
	arg++;
	
	cout << "<" << argv[arg] << ">\n";
	double doc_logp_flag_cutoff = stod(argv[arg]); //significance cutoff for single read features
	arg++;
	
	cout << "args parsed!\n";

	int mp_lower_edge_cutoff = (int)(libsize_avg + mp_edge_z_cutoff*libsize_sdv);
	
	int sr_lower_edge_cutoff = 0;

	cout << "parsing FAI!\n";
	FAI_parser faiparse(fai_filename);
	map<string, int> clens = faiparse.parse();
	faiparse.close();
	cout << "FAI Parsed!\n";
	
	//TODO make these parameters
	int mpdoc_minimum_value = 2;
	int singleton_minimum_value = 2;
	int miso_minimum_value = 2;
	int split_read_minimum_value = 2;
	int w = 1000;

	BED_parser bedparse(bed_filename, clens, correct_direction_a, correct_direction_b, 0, doc_logp_flag_cutoff, mp_lower_edge_cutoff, sr_lower_edge_cutoff, mp_size_zscore_cutoff, libsize_avg, libsize_sdv, mpdoc_minimum_value, singleton_minimum_value, miso_minimum_value, split_read_minimum_value);
	unordered_map<string, vector<int>> feats_by_contig;
	
	
	int ctr = 0;
	while(bedparse.has_next()){
		ctr++;
		contig cur_contig = bedparse.get_next_contig();
		cout << "Analyzing Contig "<< ctr << ": "<< cur_contig.name<<"\n";
		cur_contig.finish();
		feature_calculator featfinder(cur_contig, w);
		featfinder.calculate();
		vector<int> feats = featfinder.get_features();
		feats_by_contig[cur_contig.name].push_back(cur_contig.len);
		for(int e: feats) feats_by_contig[string(cur_contig.name)].push_back(e);
	}
	bedparse.close();

	unordered_map<string, vector<int>>::iterator it = feats_by_contig.begin();

        fprintf(feature_file, "Contig_Name\tContig_Len\tLo_DOC\tHi_DOC\tLo_MP_DOC\tHi_MP_DOC\tShort_MP_DOC\tLong_MP_DOC\tSingleton_DOC\tMisoriented_DOC\n");	
	while(it != feats_by_contig.end())
	{
		fprintf(feature_file, "%s", it->first.c_str());
		for(int e: it->second) fprintf(feature_file, "\t%d", e);
		fprintf(feature_file, "\n");
    	it++;
	}

	fclose(feature_file);
}
