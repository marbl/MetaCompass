import argparse
import shutil
import json
import os
from os.path import basename
from zipfile import ZipFile

def open_file(f):

    try:
        open_file = open(f, 'r')
    except:
        return None

    data = open_file.readlines()
    open_file.close()

    return data

# consolidate header info
def get_header_info(forward, reverse, unpaired, n_forward, n_reverse, n_unpaired, db, ref_filter, filter_kmer, depth, breadth, cull_filter):

    header = {}
    header['reads'] = {}
    header['summary'] = {}
    header['parameters'] = {}

    if forward != "NA" and reverse != "NA":
        header['reads']['forward'] = forward
        header['reads']['reverse'] = reverse
        header['summary']['num_forward'] = f'{int(n_forward):,}'
        header['summary']['num_reverse'] = f'{int(n_reverse):,}'
    else:
        header['reads']['forward'] = 'NA'
        header['reads']['reverse'] = 'NA'
        header['summary']['num_forward'] = 'NA'
        header['summary']['num_reverse'] = 'NA'
    
    if unpaired != "NA":
        header['reads']['unpaired'] = unpaired
        header['summary']['num_unpaired'] = f'{int(n_unpaired):,}'
    else:
        header['reads']['unpaired'] = 'NA'
        header['summary']['num_unpaired'] = 'NA'


    header['parameters']['ref_db'] = db
    header['parameters']['ref_filtering'] = ref_filter
    header['parameters']['ref_filtering_kmer_size'] = filter_kmer
    header['parameters']['depth_of_cov'] = depth
    header['parameters']['breadth_of_cov'] = breadth
    header['parameters']['ref_culling_kmer_size'] = cull_filter

    return header

# extract reference info
def get_ref_info(ref_info_file, ref_metrics_file):

    refs = open_file(ref_info_file)
    ref_info = {}

    # extract ref info from json file
    for ref in refs:

        ref = ref.strip()

        # expected file type: jsonl
        # each line is a json object
        ref_data = json.loads(ref)
        ref_id = ref_data['assemblyInfo']['assemblyAccession']

        ref_info[ref_id] = {}
        ref_info[ref_id]['id'] = ref_id
        ref_info[ref_id]['name'] = ref_data['organismName']
        ref_len = int(ref_data['assemblyStats']['totalSequenceLength'])
        ref_info[ref_id]['length'] = f'{ref_len:,}'

    # extract ref metrics from bowtie2 alignment
    ref_metrics = open_file(ref_metrics_file)
        
    for ref in ref_metrics[1:]:
        data = ref.split('\t')

        ref_id = data[0].strip()
        num_seqs = data[1].strip()
        depth = data[2].strip()
        breadth = data[3].strip()

        try:
            ref_info[ref_id]['num_seqs_covered'] = f'{int(num_seqs):,}'
        except:
            ref_info[ref_id]['num_seqs_covered'] = 'NA'
        
        try:
            ref_info[ref_id]['depth'] = f'{float(depth):,}'
        except:
            ref_info[ref_id]['depth'] = 'NA'

        try:
            ref_info[ref_id]['breadth'] = f'{float(breadth):,}'
        except:
            ref_info[ref_id]['breadth'] = 'NA'

    return ref_info

# extract reference info of culled references
def get_min_refs(ref_info, min_refs):

    refs = open_file(min_refs)
    culled_ref_info = {}

    for ref in refs:
        ref = ref.strip()
        culled_ref_info[ref] = ref_info[ref]

    return culled_ref_info

# calculate contig metrics for a given assembly file
def contig_metrics(contig_file):

    contigs = {}

    lines = open_file(contig_file)

    num_contigs = 0
    total_length = 0
    contig_lens = []
    min_contig_len = 0
    max_contig_len = 0
    ng50 = 0

    # default values for empty contigs file
    if lines is None or not lines:
        print('Got empty contig file: ', contig_file)
        return None

    # get contigs
    curr_key = ''
    for line in lines:
        line = line.strip()

        # skip blank lines
        if not line:
            continue

        if line.startswith('>'):
            curr_key = line
            contigs[curr_key] = ''
            continue
        
        contigs[curr_key] += line

    # get contig lengths
    for k,v in contigs.items():
        num_contigs += 1
        
        contig_len = len(v)
        total_length += contig_len
        contig_lens.append(contig_len)

    half_total_len = int(total_length / 2)
    contig_lens.sort(reverse=True)

    # calculate ng50
    curr_len = 0
    for i in contig_lens:
        curr_len += i

        if curr_len >= half_total_len:
            ng50 = curr_len
            break

    if len(contig_lens) > 0:
        min_contig_len = min(contig_lens)
        max_contig_len = max(contig_lens)

    return num_contigs, total_length, min_contig_len, max_contig_len, ng50

# get reference-assembly contig info
def get_ref_contig_info(contigs, min_ref_file, ref_info):

    min_refs = open_file(min_ref_file)

    ref_contigs = {}
    num_refs = 0
    ref_contigs_total_num = 0
    ref_contigs_total_length = 0
    ref_contigs_max_length = 0
    ref_contigs_min_length = -1

    for ref in min_refs:

        ref = ref.strip()
        data = contig_metrics(contigs + '/' + ref + '/contigs.pilon.fasta')

        if data is None:
            continue
            
        num_contigs, total_length, min_contig_len, max_contig_len, ng50 = data

        ref_contigs_total_num += num_contigs
        ref_contigs_total_length += total_length

        if max_contig_len > ref_contigs_max_length:
            ref_contigs_max_length = max_contig_len

        if min_contig_len < ref_contigs_min_length or ref_contigs_min_length == -1:
            ref_contigs_min_length = min_contig_len

        # seq_info = seqInfo(contigs + '/' + ref + '/alignments.sam')

        ref_contigs[ref] = {}
        ref_contigs[ref]['id'] = ref
        ref_contigs[ref]['name'] = ref_info[ref]['name']
        ref_contigs[ref]['length'] = ref_info[ref]['length']
        ref_contigs[ref]['num_contigs'] = f'{num_contigs:,}'
        ref_contigs[ref]['assembly_length'] = f'{total_length:,}'
        ref_contigs[ref]['max_contig_length'] = f'{max_contig_len:,}'
        ref_contigs[ref]['min_contig_length'] = f'{min_contig_len:,}'
        ref_contigs[ref]['ng50'] = f'{ng50:,}'

        num_refs += 1

    # consolidate info
    ref_contig_info = {}
    ref_contig_info['summary'] = {}
    ref_contig_info['summary']['num_refs'] = f'{num_refs:,}'
    ref_contig_info['summary']['total_num_contigs'] = f'{ref_contigs_total_num:,}'
    ref_contig_info['summary']['total_length_of_assembly'] = f'{ref_contigs_total_length:,}'
    ref_contig_info['summary']['max_contig_length'] = f'{ref_contigs_max_length:,}'
    ref_contig_info['summary']['min_contig_length'] = f'{ref_contigs_min_length:,}'
    ref_contig_info['references'] = list(ref_contigs.values())

    return ref_contig_info

def get_contig_info(contigs):

    data = contig_metrics(contigs)

    if data is None:
        return None
    
    num_contigs, total_length, min_contig_len, max_contig_len, ng50 = data

    contig_info = {}
    contig_info['total_num_contigs'] = f'{num_contigs:,}'
    contig_info['total_length_of_assembly'] = f'{total_length:,}'
    contig_info['max_contig_length'] = f'{max_contig_len:,}'
    contig_info['min_contig_length'] = f'{min_contig_len:,}'

    return contig_info

def write_output(header_info, ref_info, culled_ref_info, ref_contig_info, denovo_contig_info, merged_contig_info, out):

    ref_info_summary = {}
    ref_info_summary['num_selected'] = len(ref_info)
    ref_info_summary['references'] = list(ref_info.values())

    culled_info_summary = {}
    culled_info_summary['num_culled'] = len(culled_ref_info)
    culled_info_summary['references'] = list(culled_ref_info.values())

    data = {}
    data['header'] = header_info
    data['reference_selection'] = ref_info_summary
    data['reference_culling'] = culled_info_summary
    data['ref_assembly'] = ref_contig_info
    data['denovo_assembly'] = denovo_contig_info
    data['merged_assembly'] = merged_contig_info

    out_file = open(out, 'w')
    json.dump(data, out_file)
    out_file.close()

# Zip the files from given directory that matches the filter
def zip_output(dir_to_zip, zip_name, filter):
    prefix_to_remove = dir_to_zip.removesuffix('/report_generation')
    # create a ZipFile object
    with ZipFile(zip_name, 'w') as zipObj:
        # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(dir_to_zip):
            for filename in filenames:
                if filter(filename):
                    # create complete filepath of file in directory
                    filePath = os.path.join(folderName, filename)
                    # Add file to zip
                    zipObj.write(filePath, filePath.removeprefix(prefix_to_remove))

    print('All file successfully zipped to ', zip_name)

def seqInfo(alignments):

    read_alignments = open_file(alignments)
    seq_info = {}

    for i in range(len(read_alignments)):
        # remove header info
        if read_alignments[i][0] == '@':
            read_alignments.remove(read_alignments[i])
            continue
        break

    for read_alignment in read_alignments:
        info = read_alignment.split('\t')

        seq_id = info[0]
        contig = info[2]
        
        if seq_info.has_key(contig):
            seq_info[contig] += ',' + seq_id
        else:
            seq_info[contig] = seq_id

    return seq_info


def main():

    parser = argparse.ArgumentParser(description='HTML report generation')
    
    parser.add_argument('-f', '--forward', help='Forward Reads.', required=False)
    parser.add_argument('-r', '--reverse', help='Reverse Reads.', required=False)
    parser.add_argument('-u', '--unpaired', help='Unpaired Reads.', required=False)
    parser.add_argument('-f_reads', '--forward_reads', help='Number of forward reads.', required=False)
    parser.add_argument('-r_reads', '--reverse_reads', help='Number of reverse reads.', required=False)
    parser.add_argument('-u_reads', '--unpaired_reads', help='Number of unpaired reads.', required=False)

    parser.add_argument('-db', '--ref_db', help='Reference database path', required=True)
    parser.add_argument('-filter', '--ref_filter', help='Whether reference filtering was use. True/False.', required=True)
    parser.add_argument('-fil_kmer', '--filter_kmer_size', help='Kmer size used for reference filtering.', required=True)
    parser.add_argument('-cull_kmer', '--culling_kmer_size', help='Kmer size used for culling references.', required=True)
    parser.add_argument('-depth', '--depth_of_cov', help='Depth of coverage cutoff for clustering references.', required=True)
    parser.add_argument('-breadth', '--breadth_of_cov', help='Breadth of coverage cutoff for clustering references.', required=True)

    parser.add_argument('-refInfo', '--references_info', help='Path to "assembly_data_report.jsonl" file downloaded in Reference-Selection step with the datasets tool.', required=True)
    parser.add_argument('-refMetrics', '--references_metrics', help='Path to tsv containing alignment metrics of references from Reference-Selection with input reads.', required=True)
    parser.add_argument('-mrefs', '--min_ref_list', help='List of references used in Reference-Guided Assembly.', required=True)
    parser.add_argument('-refc', '--ref_contigs', help='Path to folder containing folders named by each reference used, which contain contig files for that particular assembly. From reference-guided assembly.', required=True)
    parser.add_argument('-d', '--denovo_contigs', help='Path to contigs file generated from Denovo Assembly.', required=True)
    parser.add_argument('-m', '--merged_contigs', help='Path to contigs file generated from Assembly Merge.', required=True)
    parser.add_argument('-frc', '--frc_curves', help='Path to ref-guided FRC curves.', required=True)
    parser.add_argument('-t', '--templates', help='Path to templates/script for report.', required=True)
    parser.add_argument('-json', '--out_json', help='Output json file. Contains all data to be loaded into html report.', required=True)
    parser.add_argument('-o', '--out', help='Report_Generation output folder path. Used to zip all files necessary for html report.', required=True)

    ######### parse args
    args = parser.parse_args()
    ref_contigs = args.ref_contigs
    denovo_contigs = args.denovo_contigs
    merged_contigs = args.merged_contigs
    templates = args.templates
    # TODO: Parameter Checks

    ######### Header Info #########
    header_info = get_header_info(args.forward, args.reverse, args.unpaired, \
                                  args.forward_reads, args.reverse_reads, args.unpaired_reads, \
                                  args.ref_db, args.ref_filter, args.filter_kmer_size, \
                                  args.depth_of_cov, args.breadth_of_cov, args.culling_kmer_size)

    ######### Reference-Selection Info #########
    ref_info = get_ref_info(args.references_info, args.references_metrics)

    ######### Reference-Culling Info #########
    culled_ref_info = get_min_refs(ref_info, args.min_ref_list)

    ######### Reference-Guided Assembly Metrics #########
    ref_contig_info = get_ref_contig_info(ref_contigs, args.min_ref_list, culled_ref_info)

    ######### Denovo Assembly Metrics #########
    denovo_contig_info = get_contig_info(denovo_contigs)

    ######### Merged Assembly Metrics #########
    merged_contig_info = get_contig_info(merged_contigs)

    ######### Export all info to output file #########
    write_output(header_info, ref_info, culled_ref_info, ref_contig_info, denovo_contig_info, merged_contig_info, args.out_json)

    ######### Zip report file and necessary scripts #########
    shutil.copy(templates + '/template.html', args.out + '/report_metrics/out.html')
    shutil.copy(templates + '/controller.js', args.out + '/report_metrics')
    zip_output(args.out, args.out + '/report.zip', lambda name : '.tsv' in name or '.png' in name or '.js' in name or '.html' in name or '.json' in name)

if __name__ == "__main__":
    main()