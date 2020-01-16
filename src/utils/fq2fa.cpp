#include <iostream>
#include "stdc++.h"
using namespace std; 
using std::endl;
using std::cout;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istringstream;
#include <fstream>
#include <sstream>
#include <array>
#include <ctime>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <iterator>
#include <string>
#include <string>

static void show_usage(string name)
{
    std::cerr << "Usage: " << name << " <option(s)> "
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
	          << "\t-i,--input \t\tSpecify fastq input\n"	  
              << "\t-o,--output \t\tSpecify fasta output\n"
              << std::endl;
}

void read_fastq_file(string fasta,string fastq)
{
	//input files
	ifstream infile(fasta);
	if (!infile.is_open()) {
		cerr << "Could not open input fastq file " << fasta << endl;
		exit(1);
	}
	
	//output file
	ofstream outfile(fastq);
	if (!outfile.is_open()) {
		cerr << "Could not open output fasta file " << fastq << endl;
		exit(1);
	}

	//string line;
	size_t  count = 0;
	string line;	
	while (getline(infile, line))
	{
		count+=1;
		// Process str
		if (count % 4==1){
		line.erase(0, 1);
		outfile << ">"<< line <<endl;
		}
		if (count % 4==2){
		outfile << line <<endl;
		}		
	}
	outfile.close();

}

int main(int argc, char *argv[])//, Cmdopt*) 
{ 
	// Check the number of parameters
    if (argc < 5) {//in out
        show_usage(argv[0]);
        return 1;
    }
    std::vector <std::string> sources;
    std::string fastqfile,fastafile;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-i") || (arg == "--input")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                fastqfile =argv[i+1];//argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
				i++;
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--input option requires one argument." << std::endl;
                return 1;
            }  
        } else if ((arg == "-o") || (arg == "--output")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                fastafile =argv[i+1];//argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
				i++;
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--output option requires one argument." << std::endl;
                return 1;
            }  
        } else {
            sources.push_back(argv[i]);
        }
    }

	read_fastq_file(fastqfile,fastafile);

	return 0;  
}
