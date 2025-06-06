#!/bin/env python

import glob
import time
import argparse
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

# def get_read_kmers(kmer_files):

#     read_files = glob.glob(kmer_files + '/*')
#     kmers = []
#     kmer_counts = []

#     for read_file in read_files:
#         print('Reading in read_file: ', read_file)
#         with open(read_file) as kmer_file:  # Will close handle cleanly
#             for count, sequence in SimpleFastaParser(kmer_file):

#                 try:
#                     kmers.append(sequence)
#                     kmer_counts.append(int(count))
#                 except:
#                     print('Kmer: ', sequence)
#                     print('Count: ', count)

#     s1 = pd.Series(kmers, name='kmer')
#     s2 = pd.Series(kmer_counts, name='counts')
#     df = pd.DataFrame(dict(kmer=s1, reads=s2))
#     return df

def get_read_kmers(kmer_file_path):

    kmers = []
    kmer_counts = []

    with open(kmer_file_path) as kmer_file:  # Will close handle cleanly
        for count, sequence in SimpleFastaParser(kmer_file):

            try:
                # kmers[sequence] = int(count)
                kmers.append(sequence)
                kmer_counts.append(int(count))
            except:
                print('Kmer: ', sequence)
                print('Count: ', count)

    # kmers['index'] = 'reads'
    s1 = pd.Series(kmers, name='kmer')
    s2 = pd.Series(kmer_counts, name='counts')
    df = pd.DataFrame(dict(kmer=s1, reads=s2))
    # print(df)
    # df = df.set_index(['kmer'])
    # print(df)
    return df

def get_kmers_to_keep(read_kmers_set, ref_kmer_file_path):

    kmers = []

    with open(ref_kmer_file_path) as kmer_file:  # Will close handle cleanly
        for count, sequence in SimpleFastaParser(kmer_file):

            try:
                kmers.append(sequence)
            except:
                print('Kmer: ', sequence)
                print('Count: ', count)

    return read_kmers_set.intersection(set(kmers))

def get_ref_kmers(kmer_file_path, ref_name):

    kmers = []
    kmer_counts = []

    with open(kmer_file_path) as kmer_file:  # Will close handle cleanly
        for count, sequence in SimpleFastaParser(kmer_file):


            try:
                kmers.append(sequence)
                kmer_counts.append(int(count))
            except:
                print('Kmer: ', sequence)
                print('Count: ', count)

    s1 = pd.Series(kmers, name='kmer')
    s2 = pd.Series(kmer_counts, name='counts')
    df = pd.DataFrame(dict(kmer=s1, ref_name=s2))
    df = df.rename({'ref_name' : ref_name}, axis=1)
    return df


def get_min_set(kmer_df, all_refs, min_matches):

    min_set = []

    while True:
        print(kmer_df)

        # Case 1. No more read kmers
        if sum(kmer_df['reads']) == 0:
            break

        curr_best_overlap = None
        curr_best_ref = None
        curr_best_ref_count = 0

        # get next best ref
        for ref in all_refs:
            
            leftover = kmer_df[ref] - kmer_df['reads']
            leftover[leftover < 0] = 0
            temp_overlap = kmer_df[ref] - leftover
            temp_count = sum(temp_overlap)

            if temp_count >= min_matches and temp_count > curr_best_ref_count:
                curr_best_overlap = temp_overlap
                curr_best_ref = ref
                curr_best_ref_count = temp_count

        # Case 2. No more ref kmers
        if curr_best_ref is None:
            break

        kmer_df['reads'] = kmer_df['reads'] - curr_best_overlap
        kmer_df.clip(lower=0)
        kmer_df = kmer_df.drop(curr_best_ref, axis=1)
        all_refs.remove(curr_best_ref)

        min_set.append(curr_best_ref)

        print('+ Ref Chosen: ' + curr_best_ref + '. Shared kmers = ' + str(curr_best_ref_count))

    return min_set


def write_min_set(min_set, out_dir):

    out_file_path = out_dir + '/min_reference_candidates.txt'
    out_file = open(out_file_path, 'w')

    for ref in min_set:
        out_file.write(ref + '\n')

    out_file.close()


def main():

    parser = argparse.ArgumentParser(description='Reference Culling')
    
    parser.add_argument('-m', '--min_matches', help='Minimum number of kmer matches to consider a reference for assembly.', required=True)
    parser.add_argument('-r', '--read_kmers', help='Path to file containing read kmers, in fasta format.', required=True)
    parser.add_argument('-c', '--ref_kmers', help='Path to directory containing reference kmer files, all in fasta format.', required=True)
    parser.add_argument('-o', '--out', help='Path to output directory. Minimum list of references will be written here.', required=True)


    start_time = time.time()

    ######### parse args
    args = parser.parse_args()
    min_matches = int(args.min_matches)
    read_kmers_file = args.read_kmers
    ref_kmers_files = glob.glob(args.ref_kmers + '/*.fasta')
    out_dir = args.out

    ######### read kmers into dataframe
    all_refs = []

    # TODO: convert kmer counts to int32 
    kmer_df = get_read_kmers(read_kmers_file)
    print('kmer_df: ', kmer_df)

    ### Need to keep read_kmers that are in reference_kmers only
    read_kmers_set = set(kmer_df['kmer'])
    kmers_to_keep = set()
    print('kmers to keep: ', kmers_to_keep)

    for kmer_file_path in ref_kmers_files:
        ref_name = kmer_file_path.split('/')[-1]
        kmers_to_keep.update(get_kmers_to_keep(read_kmers_set, kmer_file_path))
        print('kmers to keep: ', kmers_to_keep)


    ### Remove read_kmers not in reference_kmers
    print('Number of rows without kmers removed:', len(kmer_df.index))
    kmers_to_keep = list(kmers_to_keep)
    kmer_df = kmer_df.set_index('kmer')
    kmer_df = kmer_df.T
    print(kmer_df)
    kmer_df = kmer_df[kmers_to_keep]
    kmer_df = kmer_df.T
    kmer_df.reset_index(inplace=True)
    kmer_df = kmer_df.rename(columns = {'index':'kmer'})
    print(kmer_df)
    print('Number of rows after kmers removed:', len(kmer_df.index))
    del kmers_to_keep # free memory
    del read_kmers_set # free memory


    for kmer_file_path in ref_kmers_files:
        ref_name = kmer_file_path.split('/')[-1]
        ref_kmers = get_ref_kmers(kmer_file_path, ref_name)
        all_refs.append(ref_name)

        # TODO: Every time we do this join, NaN values take place for kmers not present in the newly add ref_kmer column
        # NaN value dtype is float64. Must first fillna(0) to convert these NaN's to 0.0 (float64).
        # Then, convert them all to int32. There is no integer implementation of NaN values in Pandas at the moment.
        kmer_df = kmer_df.join(ref_kmers.set_index('kmer'), on='kmer')
        print(kmer_df)

    print('FINISHED LOADING ALL REFERENCE KMERS!')
    del ref_kmers_files # free memory
    kmer_df = kmer_df.set_index('kmer')
    kmer_df = kmer_df.fillna(0)
    print(kmer_df)
    ######### get min set of references
    min_set = get_min_set(kmer_df, all_refs, min_matches)
    print(min_set)

    ######### write min set
    write_min_set(min_set, out_dir)

    end_time = time.time()
    print('FINISH REFERENCE CULLING. REFERENCE CULLING TAKES A TOTAL OF ' + str(end_time - start_time) + ' SECONDS')

if __name__ == "__main__":
    main()