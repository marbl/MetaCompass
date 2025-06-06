import argparse
from test_definitions import tests


# TODO have input args store file type enum automatically in file_type
def get_parser():

    parser = argparse.ArgumentParser(
        description="Create metacompass test cases from given input data. Input files in FASTQ format are required and can be specified with --paired, --interleaved, or --unpaired. Test type is also required."
    )
    input_files = parser.add_mutually_exclusive_group(required=True)
    input_files.add_argument(
        "--paired",
        "-p",
        type=str,
        nargs=2,
        dest="paired",
        action="store",
        help="A pair of fastq files representing mate pair reads",
    )
    input_files.add_argument(
        "--interleaved",
        "-i",
        type=str,
        dest="interleaved",
        action="store",
        help="An interleaved fastq file with paired reads",
    )
    input_files.add_argument(
        "--unpaired",
        "-u",
        type=str,
        dest="unpaired",
        action="store",
        help="A fastq file with unpaired reads.",
    )

    parser.add_argument(
        "--test_name",
        "-t",
        type=str,
        dest="test_name",
        action="store",
        required=True,
        choices=tests.available_test_names,
        help=str(
            "REQUIRED. The type of test data to generate. Options are "
            + str(tests.available_test_names)
        ),
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="fastq_out",
        action="store",
        required=False,
        help=str(
            "The base name of the output file. For example, if you put -o my_output then the output will be my_output_1.fq and my_output_2.fq"
        ),
    )

    parser.add_argument(
        "--seed",
        "-s",
        type=str,
        dest="seed",
        action="store",
        required=False,
        help=str(
            "The seed for any operations that rely on randomness, such as shuffling"
        ),
    )

    parser.add_argument(
        "--substring_prop",
        type=float,
        dest="substr_p",
        action="store",
        required=False,
        default=0.5,
        help=str(
            "For the substring test: The value specifying the length of the substrings. It's a proportion (0.0, 1.0). Substrings will be taken from the center of the reads. Default is 0.5"
        ),
    )

    return parser
