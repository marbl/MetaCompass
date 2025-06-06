import random

# Static class defining how we handle shuffling the order of a data set
class shuffle:

    name = "shuffle"

    def apply(fastq_data, args):
        if args.seed is not None:
            random.seed(args.seed)
        random.shuffle(fastq_data)
        return
