basic usage example:

python generate_metamorphic_test_data.py --paired test1.fq test2.fq --test_name shuffle -o my_shuffle

#AUTO GENERATED USAGE INFORMATION
usage: generate_metamorphic_test_data.py [-h]
                                         (--paired PAIRED PAIRED | --interleaved INTERLEAVED | --unpaired UNPAIRED)
                                         --test_name {shuffle,substring}
                                         [--output FASTQ_OUT] [--seed SEED]
                                         [--substring_prop SUBSTR_P]

Create metacompass test cases from given input data. Input files in FASTQ
format are required and can be specified with --paired, --interleaved, or
--unpaired. Test type is also required.

optional arguments:
  -h, --help            show this help message and exit
  --paired PAIRED PAIRED, -p PAIRED PAIRED
                        A pair of fastq files representing mate pair reads
  --interleaved INTERLEAVED, -i INTERLEAVED
                        An interleaved fastq file with paired reads
  --unpaired UNPAIRED, -u UNPAIRED
                        A fastq file with unpaired reads.
  --test_name {shuffle,substring}, -t {shuffle,substring}
                        REQUIRED. The type of test data to generate. Options
                        are ['shuffle', 'substring']
  --output FASTQ_OUT, -o FASTQ_OUT
                        The base name of the output file. For example, if you
                        put -o my_output then the output will be
                        my_output_1.fq and my_output_2.fq
  --seed SEED, -s SEED  The seed for any operations that rely on randomness,
                        such as shuffling
  --substring_prop SUBSTR_P
                        For the substring test: The value specifying the
                        length of the substrings. It's a proportion (0.0,
                        1.0). Substrings will be taken from the center of the
                        reads. Default is 0.5
