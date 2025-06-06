#include <iostream>
#include <cstring>
#include <cstdlib>
#include <regex>
#include <sys/stat.h>
#include "cmdopt.hpp"

using namespace std;

bool isPostiveNumber(const std::string &x) {
    if (strcmp(x.c_str(), "0") == 0) return false;
    std::regex e ("^\\d+");
    if (std::regex_match (x,e)) return true;
    else return false;
}

bool isPathExist(const std::string &s) {
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}

// set default values
Cmdopt::Cmdopt() {
    min_shared_kmer = 1;
    reads = "";
    output_folder = "";
}   

void Cmdopt::refCullingHelp(int status){
    cerr << "Usage:" << endl;
    cerr << "        ./refcull"
         << " -r [read kmer file]" 
         << " -c [reference kmer folder_1]"
         << " -c [reference kmer folder_2"
         << " ... -c [reference kmer folder_n]"
         << " -m [min shared kmers]"
         << " -o [output folder]"
         << endl;
    cerr << "Options:" << endl;
    cerr << "        --help  "
         << endl; 
    cerr << "        -r/--reads  "
         << "The kmer file (output of Jellyfish) of all reads. "
         << "Provide only 1 read kmer file for all the reads."
         << endl; 
    cerr << "        -c/--reference-folder  "
         << "Folder that contains the kmer files (output of Jellyfish) of reference candidates. "
         << "Can provide 1 or many reference kmer folders. Each folder can contain many kmer files."
         << endl;
    cerr << "        -m/--min-shared-kmer  "
         << "Minimum number of shared kmers. "
         << "Reference will be chosen if sharing at least this number of kmers with the reads. "
         << "(default: 1)"
         << endl;
    cerr << "        -o/--output-folder  "
         << "Output folder that will contain the output file of ref culling."
         << endl;
    cerr << endl;
    exit(status);
}

void Cmdopt::parseCmdOptions(int argc, char* argv[], Cmdopt* cmdopt){
    if (argc == 2 && strcmp(argv[1], "--help") == 0)
        refCullingHelp(0);
    
    for (int i = 1; i < argc - 1; i++){
        if (strcmp(argv[i], "--reads") == 0 ||
            strcmp(argv[i], "-r") == 0) {
            if (!isPathExist(argv[++i])) {
                cerr << "Error: Read kmer file was not found" << endl;
                refCullingHelp(1);
            }   
            cmdopt->reads = argv[i];
        }
        if (strcmp(argv[i], "--reference-folder") == 0 ||
            strcmp(argv[i], "-c") == 0) {
            if (!isPathExist(argv[++i])) {
                cerr << "Error: Reference candidate kmer folder was not found" << endl;
                refCullingHelp(1);
            }   
            cmdopt->reference_folders = argv[i];
        }
        if (strcmp(argv[i], "--min-shared-kmer") == 0 ||
            strcmp(argv[i], "-m") == 0) {
            if (!isPostiveNumber(argv[++i])) {
                cerr << "Error: Min shared kmers should be a positive number" << endl;
                refCullingHelp(1);
            }   
            cmdopt->min_shared_kmer = atoi(argv[i]);
        }
        if (strcmp(argv[i], "--output-folder") == 0 ||
            strcmp(argv[i], "-o") == 0) {
            if (!isPathExist(argv[++i])) {
                cerr << "Error: Output folder was not found" << endl;
                refCullingHelp(1);
            } 
            cmdopt->output_folder = argv[i];
        } 
         if (strcmp(argv[i], "--help") == 0) 
             refCullingHelp(1);
    }
    
    // Read kmer file and reference kmer folder must be specified
    if (cmdopt->reads == ""){
        cerr << "Error: Read kmer file were not provided" << endl;
        refCullingHelp(1);
    }
    if (cmdopt->reference_folders == ""){
        cerr << "Error: Reference candidate kmer folder was not provided" << endl;
        refCullingHelp(1);
    }
    if (cmdopt->output_folder == ""){
        cerr << "Error: Output folder was not provided" << endl;
        refCullingHelp(1);
    }
}

