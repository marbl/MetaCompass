+ Check [here](https://gitlab.umiacs.umd.edu/mpop/metacompass/-/wikis/Software/Testing#dataset-1) for an analysis of this dataset.
#### Testing Steps:
+ Step 1: Create executable if you have not yet created it
   + You are now in `scripts/cull_reference_candidates/test/dataset1`
   + `cd ../../`: Navigate back to cull_reference_candidates to compile.
   + `g++ -std=c++11 cmdopt.cpp cull_reference_candidate.cpp process_map.cpp -o refcull` to compile.
   + Our executable is `refcull`
+ Step 2: Run test
   + `cd test/dataset1`: Navigate to test folder to run test
   + run `./set1_script.sh`
   + Result of ref culling is in `./min_reference_candidates.txt`

#### Content of each file/folder
+ `set1_original_reads` : folder of the original reads
+ `set1_read_kmer.txt` : kmer file of all the reads in set1_original_reads 
+ `set1_original_refs` : folder of all the original reference candidates
+ `set1_ref_kmer_folders` : kmer folder of all the references in set1_original_refs
+ `set1_script.sh` : script to run the test of dataset 1

