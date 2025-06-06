#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <dirent.h> // for processing directory
#include <unordered_map>
#include <vector>
#include <set>
#include <bits/stdc++.h>
#include "cmdopt.hpp"
#include "process_map.hpp"
using namespace std;

// map to store the kmers of reads
unordered_map<string, int> reads_kmers_map = {};

vector<string> chosen_references; /* store the chosen min set of reference candidates */

vector<string> ref_files; /* store all reference filenames */

KmerSet shared_kmers_set; /* For each pair in the set:
                            + 1st component: num of shared kmers b/w read and 
                            the reference candidate
                            + 2nd component: reference candidate name */

int main (int argc, char *argv[]) {
    time_t nowtime;
    time_t time_start = time(NULL);
    cout << "START REFERENCE CULLING AT TIME: "
        << ctime(&time_start)
        << endl;
    
    // parse command options
    Cmdopt cmdopt = Cmdopt();
    cmdopt.parseCmdOptions(argc, argv, &cmdopt);

    // Check if output folder exists
    if (!opendir(cmdopt.output_folder.c_str())){
        cerr << "Error: Cannot find output folder " 
            << cmdopt.output_folder
            << endl;
        exit(1);
    }

    // Check if reference kmer folder exists
    string reference_folder = cmdopt.reference_folders.c_str();
    auto ref_folder = opendir(cmdopt.reference_folders.c_str());
    if (!ref_folder){
        cerr << "Error: Cannot find reference candidate kmer folder" 
            << reference_folder
            << endl;
        exit(1);
    }

    // 1. Need to load all read kmers first
    // 2. For each reference:
    // 2a. Load reference kmers
    // 2b. Check for matches in read kmers and record
    // 2c. Keep track of current best kmer, its # matches and kmers
    // 2d. De-allocate memory for current reference kmers if not current best match

    /**********************************************************************/
    /****************** STEP 1: STORE READ KMERS IN MAP *******************/
    /**********************************************************************/
    
    // Store the kmers of all reads in reads_kmers_map
    // Key: kmer of read, value: corresponding occurence
    store_kmers_in_map(cmdopt.reads, &reads_kmers_map);
    cout << "FINISHED LOADING ALL READ KMERS!" << endl;

    /**********************************************************************/
    /**************************** FINISH STEP 1 ***************************/
    /**********************************************************************/

    /**********************************************************************/
    /******************** STEP 2:  REFERENCE CULLING **********************/
    /**********************************************************************/

    while (auto ref_file = readdir(ref_folder)) {
        
        auto ref_name = ref_file->d_name;
        
        // skip hidden files && skip non-regular file
        if (!ref_name || ref_name[0] == '.'){
            continue;
        }
        if (ref_file->d_type != DT_REG){
            continue;
        }

        ref_files.push_back(ref_name);
    }

    string ref_name;

    string curr_ref;
    int curr_ref_idx;
    int curr_num_matches;
    unordered_map<string, int> curr_ref_kmers = {};

    string temp_ref = "";
    int temp_ref_idx;
    int temp_num_matches = 0;
    unordered_map<string, int> temp_ref_kmers = {};

    int inner = 0;
    int outer = 0;

    while (true){

        // cout << "Start Outer: " << outer << endl;
        inner = 0;

        curr_ref = "";
        curr_ref_idx = -1;
        curr_num_matches = 0;
        curr_ref_kmers.clear();

        /* Read in references one by one */
        for (int i = 0; i < ref_files.size(); i++) {

            ref_name = ref_files[i];


            // cout << "Start Inner: " << inner << endl;
            
            // cout << "Filename: " << ref_name << endl;

            // TODO: grabs all files in reference_folder but may read non-kmer (non-fasta) files.
            // store_kmers_in_map only looks for beginning ">" on line[0] and then seq on line[1].
            // if a non-kmer file happens to have this structure, it will be read in.
            //
            // To mitigate, make sure that folder contains only ref kmer files.
            string ref_path = reference_folder + "/" + ref_name;

            // Store kmers of the reference in temp_ref_kmers
            // Key: kmer, value: kmer occurrence
            temp_ref = ref_name;
            temp_ref_idx = i;
            temp_ref_kmers.clear();
            store_kmers_in_map(ref_path, &temp_ref_kmers);
            temp_num_matches = count_shared_kmers(temp_ref_kmers, reads_kmers_map);

            if (temp_num_matches > curr_num_matches && temp_num_matches >= cmdopt.min_shared_kmer){
                curr_ref = temp_ref;
                curr_ref_idx = temp_ref_idx;
                curr_num_matches = temp_num_matches;
                curr_ref_kmers = temp_ref_kmers;
            }
            // cout << "End Inner: " << inner << endl;
            inner += 1;
        }

        // closedir(ref_folder);

        if (curr_ref == "" || curr_num_matches == 0 || curr_ref_kmers.size() == 0){
            break;
        }else{
            // 1. add ref to set of chosen refs
            // 2. remove chosen ref kmers from set of read kmers
            // 3. remove chosen ref from list of references
            // cout << "For file: " << curr_ref << "; Number of shared kmers was: " << curr_num_matches << endl;
            chosen_references.push_back(curr_ref);
            delete_shared_kmers_in_read(curr_ref_kmers, reads_kmers_map);
            ref_files.erase(ref_files.begin() + curr_ref_idx);

            cout << " + Ref chosen: " << curr_ref 
                << ". Shared kmers = " << curr_num_matches << endl;

        }

        // cout << "End Outer: " << outer << endl;
        outer += 1;
    }

    /**********************************************************************/
    /**************************** FINISH STEP 2 ***************************/
    /**********************************************************************/

    /**********************************************************************/
    /********** STEP 3: WRITE THE MIN SET OF REFERENCES TO A FILE ********/
    /**********************************************************************/

    string output_file = cmdopt.output_folder + "/min_reference_candidates.txt";

    // If the file doesn't exists, create and write to it
    // If the file exists, truncate and write to it
    ofstream out_file;
    out_file.open(output_file, ios::out | ios::trunc);
    for (string reference : chosen_references){
        // remove ".fasta", ".fq", etc from string

        // Vector of string to save tokens
        vector <string> tokens;
        // stringstream class check1
        stringstream check1(reference);
        string intermediate;
        
        // Tokenizing w.r.t. '.'
        while(getline(check1, intermediate, '.'))
        {
            tokens.push_back(intermediate);
        }

        string tmp = "";
        string result = "";
        if (tokens.size() == 2){

            result = tokens[0];

        }else{

            for (int i = 0; i < tokens.size() - 1; i++){
                tmp = tmp + tokens[i] + ".";
            }
            result  = tmp.substr(0, tmp.size() - 1); // remove last '.'
        }
        
        out_file << result << endl;
    }
    out_file.close();
    
    cout << "FINISH REFERENCE CULLING. REFERENCE CULLING TAKES A TOTAL OF "
         << time(NULL) - time_start
         << " SECONDS"
         << endl;
    /**********************************************************************/
    /**************************** FINISH STEP 3 ***************************/
    /**********************************************************************/
    return 0;
}

