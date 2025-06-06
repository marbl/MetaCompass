#!/bin/env python

import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isdir, join

def parse_data(data):
    data_file = open(data, 'r') # input data
    out_dir = sys.argv[2]      # output folder

    data = data_file.readlines()
    data_file.close()

    header = data[0].strip().split('\t')[1:]

    data_list = []

    for i in range(1, len(data), 1):

        line = data[i].split('\t')

        length = float(line[1])
        lo_doc = float(line[2])
        hi_doc = float(line[3])
        lo_mp_doc = float(line[4])
        hi_mp_doc = float(line[5])
        short_mp_doc = float(line[6])
        long_mp_doc = float(line[7])
        singleton_doc = float(line[8])
        misoriented_doc = float(line[9])

        data_list.append((length, lo_doc, hi_doc, lo_mp_doc, hi_mp_doc, short_mp_doc, long_mp_doc, singleton_doc, misoriented_doc))

    return data_list, header

def graph_features(feature, data, header, total_length, out_dir):
    # assume data is sorted correctly - in ascending order of the feature

    x = [] # feature values
    y = [] # coverage %

    curr_x = 0
    curr_y = 0

    for i in range(len(data)):

        contig_features = data[i]

        # if contig_features[feature] < 1:
        #     curr_y += contig_features[0] / total_length
        #     continue

        curr_y += contig_features[feature]
        curr_x += contig_features[0]

        x.append(curr_x)
        y.append(curr_y)

    plt.figure()
    plt.step(x, y)
    plt.ylabel('Feature Threshold: ' + header[feature])
    plt.xlabel('Cumulative Contig Length')
    plt.show()
    plt.savefig(out_dir + '/' + header[feature] + '_fig.png')
    plt.close()

def graph_all(name, data, out_dir, log_file):
    print('Graphing: ', name)

    log_file.write('######## Calculating features for ' + name + ' contigs ########\n')

    data_list, header = parse_data(data)

    data_array = np.array(data_list)
    total_length = sum(data_array[:,0])
    log_file.write('Total Length: ' + str(total_length) + '\n')
    N50 = total_length / 2
    log_file.write('N50: ' + str(N50) + '\n')

    data_list.sort() # sort via contig length - 1st element of tuple
    data_list = data_list[::-1]

    for i in range(1, len(header), 1):

        # data_list.sort(key = lambda x: x[i]) # sort by feature values in ascending order
        graph_features(i, data_list, header, total_length, out_dir) # create feature graphs

    # create Nx graph
    Nx_values = []
    curr_clen = 0
    total_length = total_length / 1000
    N50_clen = None

    for contig in data_list:
        curr_clen += contig[0] / 1000
        Nx_values.append(int((curr_clen / total_length) * 100))
        if curr_clen * 1000 >= N50 and N50_clen is None:
            N50_clen = contig[0]

    log_file.write('N50 contig length: ' + str(N50_clen) + '\n')
    log_file.write('Largest contig: ' + str(max(data_array[:,0])) + '\n')
    log_file.write('Smallest contig: ' + str(min(data_array[:,0])) + '\n')
    plt.figure()
    plt.step(Nx_values, list(np.array(data_list)[:,0] / 1000))
    plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ax = plt.gca()
    ax.xaxis.tick_bottom()                     # and move the X-Axis
    ax.yaxis.tick_left()                    # remove right y-Ticks
    plt.xlabel('Nx')
    plt.ylabel('Contig Length (Kbp)')
    plt.tight_layout()
    plt.show()
    plt.savefig(out_dir + "/Nx_fig.png")
    plt.close()
    


def main():

    parser = argparse.ArgumentParser(description='FRC Curve Generation')
    
    parser.add_argument('-r', '--ref_info', help='Path to folder containing frc info for each reference-guided assembly.', required=False)
    parser.add_argument('-d', '--denovo_info', help='Path to file containing frc info for denovo assembly.', required=False)
    parser.add_argument('-m', '--merged_info', help='Path to file containing frc info for merged assembly.', required=True)
    parser.add_argument('-o', '--out', help='Output path.', required=True)
    parser.add_argument('-l', '--log', help='Log file.', required=False)

    # parse args
    args = parser.parse_args()
    print('Got args: ', args)
    ref_info = args.ref_info
    denovo_info = args.denovo_info
    merged_info = args.merged_info
    out = args.out
    log = args.log
    # TODO: Parameter Checks
    log_file = open(log, 'w')

    # plot refguided
    ref_out = out + '/refguided'
    refguided_dirs = [ref_dir for ref_dir in listdir(ref_info) if isdir(join(ref_info, ref_dir))]
    for ref_dir in refguided_dirs:
        curr_ref_out = ref_out + '/' + ref_dir
        curr_ref_info = ref_info + '/' + ref_dir + '/' + ref_dir + '.tsv'
        graph_all(ref_dir, curr_ref_info, curr_ref_out, log_file)

    # plot denovo
    graph_all('denovo', denovo_info, out + '/denovo', log_file)

    # plot merged
    graph_all('merged', merged_info, out + '/merged', log_file)

    log_file.close()

if __name__ == "__main__":
    main()