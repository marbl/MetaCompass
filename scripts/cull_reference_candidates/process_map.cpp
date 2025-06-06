#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <dirent.h> // for processing directory
#include "process_map.hpp"

// Read the kmers from kmer_file and store in kmer_map
// Key of map: kmer sequence, value of map: corresponding occurrence
void store_kmers_in_map(string kmer_file,
                        unordered_map<string, int> * kmer_map){

    ifstream in_file;
    in_file.open(kmer_file);

    if (in_file.fail()){
        exit(1);
    }

    string line;
    while (getline(in_file, line)) {
        if (line[0] == '>') { // line of kmer occurrence
            int occurrence = atoi(line.substr(1).c_str());
            getline(in_file, line); // line of kmer sequence
            (*kmer_map)[line] = occurrence;
        }
    }

    in_file.close();
}

int count_shared_kmers(unordered_map<string, int>& ref_kmer_map,
                       unordered_map<string, int>& read_kmer_map){

    if (&ref_kmer_map == NULL || &read_kmer_map == NULL) return 0;

    int shared_kmers = 0;

    // Find the appropriate maps for "iteration" and "look up"
    unordered_map<string, int>* map_for_iter;
    unordered_map<string, int>* map_for_lookup;

    // TODO: Reference kmers should always be iterated over.
    if (ref_kmer_map.size() <= read_kmer_map.size()){
        map_for_iter = &ref_kmer_map;
        map_for_lookup = &read_kmer_map;

    } else {
        map_for_iter = &read_kmer_map;
        map_for_lookup = &ref_kmer_map;
    }
    
    /*Find the shared kmers between 2 maps: map_for_iter & map_for_lookup 
    (i.e., read kmer map & reference kmer map) */
    for (auto iter : (*map_for_iter)){
        const string& kmer = iter.first;

        if ((*map_for_lookup).find(kmer) != (*map_for_lookup).end()){

            int read_kmer_occur = read_kmer_map[kmer];
            int ref_kmer_occur = ref_kmer_map[kmer];

            if (ref_kmer_occur <= read_kmer_occur){
                shared_kmers += ref_kmer_occur;
            } else {
                shared_kmers += read_kmer_occur;
            }
        }
    }

    return shared_kmers;
}


// Find number of shared kmers between read kmer map and reference kmer map
// Update the number of shared kmers in shared_kmer_map
void count_num_shared_kmer (unordered_map<string, int>& ref_kmer_map,
                            unordered_map<string, int>& read_kmer_map,
                            KmerSet& shared_kmers_set,
                            const string& ref_name) {
    if (&ref_kmer_map == NULL || &read_kmer_map == NULL || ref_kmer_map.size() == 0 || read_kmer_map.size() == 0) return;
    
    int shared_kmers = 0;
    // Find the appropriate maps for "iteration" and "look up"
    unordered_map<string, int>* map_for_iter;
    unordered_map<string, int>* map_for_lookup;
    // TODO: Reference kmers should always be iterated over.
    if (ref_kmer_map.size() <= read_kmer_map.size()){
        map_for_iter = &ref_kmer_map;
        map_for_lookup = &read_kmer_map;
        // cout << "REFERENCE SET SMALLER!" << endl;
    } else {
        map_for_iter = &read_kmer_map;
        map_for_lookup = &ref_kmer_map;
    }
    
    /*Find the shared kmers between 2 maps: map_for_iter & map_for_lookup 
    (i.e., read kmer map & reference kmer map) */
    for (auto iter : (*map_for_iter)){
        const string& kmer = iter.first;
        if ((*map_for_lookup).find(kmer) != (*map_for_lookup).end()){
            int read_kmer_occur = read_kmer_map[kmer];
            int ref_kmer_occur = ref_kmer_map[kmer];
            // cout << "FOUND A MATCH!" << endl;
            if (ref_kmer_occur <= read_kmer_occur){
                shared_kmers += ref_kmer_occur;
            } else {
                shared_kmers += read_kmer_occur;
            }
        }
    }

    
    shared_kmers_set.insert(make_pair(shared_kmers, ref_name));
}

// Delete from the read_kmer_map the shared kmers between the ref_kmer_map and the read_kmer_map
void delete_shared_kmers_in_read (unordered_map<string, int>& ref_kmer_map,
                                  unordered_map<string, int>& read_kmer_map){
    
    if (&ref_kmer_map == NULL || &read_kmer_map == NULL || ref_kmer_map.size() == 0 || read_kmer_map.size() == 0) return;

    // Find the appropriate maps for "iteration" and "look up"
    unordered_map<string, int>* map_for_iter;
    unordered_map<string, int>* map_for_lookup;
    if (ref_kmer_map.size() <= read_kmer_map.size()){
        map_for_iter = &ref_kmer_map;
        map_for_lookup = &read_kmer_map;
        // cout << "REFERENCE SET SMALLER111!" << endl;
    } else {
        map_for_iter = &read_kmer_map;
        map_for_lookup = &ref_kmer_map;
    }

    /* Delete the shared kmers in read_kmer_map */
    for (auto iter : (*map_for_iter)){
        const string& kmer = iter.first;
        if ((*map_for_lookup).find(kmer) != (*map_for_lookup).end()){
            int read_kmer_occur = read_kmer_map[kmer];
            int ref_kmer_occur = ref_kmer_map[kmer];
            // cout << "FOUND A MATCH111!" << endl;

            if (ref_kmer_occur < read_kmer_occur){
                read_kmer_map[kmer] = read_kmer_occur - ref_kmer_occur;
            } else {
                if (!read_kmer_map.empty()) read_kmer_map.erase(kmer);
            }
        }
    }
}

// Add the chosen_ref to list of chosen_references
// Remove chosen_ref from references
void process_chosen_ref(unordered_map<string, unordered_map<string, int> *>& references,
                    vector<string>& chosen_references,
                    const string& chosen_ref){
    chosen_references.push_back(chosen_ref);
    delete references[chosen_ref];
    if (!references.empty()) references.erase(chosen_ref);
}
