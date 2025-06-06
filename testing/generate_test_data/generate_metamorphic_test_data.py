# generates test cases from given input files
# TODO: add optional seed param
# TODO: add output file name param
# TODO: add default output filename
# TODO: check that input files exist
# TODO: check that output files can be written

import sys
import argparse
import argument_parser
from fastq_parser import fastq_parser
from test_definitions import tests

# parse the inputs and call the appropriate functions
def main():

    parser = argument_parser.get_parser()
    args = parser.parse_args()

    file_type = ""
    fastq_files = ""
    fastq_out = [args.test_name + "_1.fq", args.test_name + "_2.fq"]

    # TODO handle this in arg parser
    if args.unpaired:
        file_type = "unpaired"
        fastq_files = args.unpaired

    if args.interleaved:
        file_type = "interleaved"
        fastq_files = args.interleaved

    if args.paired:
        file_type = "paired"
        fastq_files = args.paired

    if args.fastq_out is not None:
        fastq_out = [args.fastq_out + "_1.fq", args.fastq_out + "_2.fq"]

    my_fq_parser = fastq_parser(file_type)
    fastq_data = my_fq_parser.parse(fastq_files)

    my_test = tests(args)

    my_test.apply(fastq_data, args)

    my_fq_parser.write(fastq_data, fastq_out)
    # print(str(fastq_data[3]))
    # print(str(fastq_data))

    return


if __name__ == "__main__":
    main()
