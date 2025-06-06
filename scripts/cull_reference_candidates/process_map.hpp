#include <unordered_map>
#include <set>
#include <vector>
#include <cstring>
using namespace std;

using KmerPair = pair<int, string>;
using KmerSet = set<KmerPair, greater<KmerPair>>; /*
                            Pairs in the set are in descending order
                            For each pair in the set:
                            + 1st component: num of shared kmers b/w read and 
                            the reference candidate
                            + 2nd component: reference candidate name */

// Read the kmers from kmer_file and store in kmer_map
// Key of map: kmer sequence, value of map: corresponding occurrence
void store_kmers_in_map(string kmer_file,
                        unordered_map<string, int> * kmer_map);

// Find number of shared kmers between read kmer map and reference kmer map
// Then return number of share kmers
int count_shared_kmers(unordered_map<string, int>& ref_kmer_map,
                       unordered_map<string, int>& read_kmer_map);


// Find number of shared kmers between read kmer map and reference kmer map
// Update the number of shared kmers in shared_kmer_set
void count_num_shared_kmer (unordered_map<string, int>& ref_kmer_map,
                            unordered_map<string, int>& read_kmer_map,
                            KmerSet& shared_kmers_set,
                            const string& ref_name);

// Delete from the read_kmer_map the shared kmers between the ref_kmer_map and the read_kmer_map
void delete_shared_kmers_in_read (unordered_map<string, int>& ref_kmer_map,
                                  unordered_map<string, int>& read_kmer_map);



// Add the chosen_ref to chosen_references
// Remove chosen_ref from references
void process_chosen_ref(unordered_map<string, unordered_map<string, int> *>& references,
                    vector<string>& chosen_references,
                    const string& chosen_ref);
