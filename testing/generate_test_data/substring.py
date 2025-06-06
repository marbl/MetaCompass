# Static class defining how to split all the reads in half


class substring:

    name = "substring"

    def substr(read, substr_p):

        rlen = len(read[1]) - 1
        slen = int(rlen * substr_p)
        begin = int((rlen - slen) / 2)
        end = rlen - int((rlen - slen + 1) / 2)
        print(str(rlen) + " " + str(slen) + " " + str(begin) + " " + str(end))
        read[1] = read[1][begin:end] + "\n"
        read[3] = read[3][begin:end] + "\n"

    def apply_paired(fastq_data, substr_p):
        for read_pair in fastq_data:
            substring.substr(read_pair[0], substr_p)
            substring.substr(read_pair[1], substr_p)

        return

    def apply_unpaired(fastq_data, substr_p):
        for read in fastq_data:
            substring.substr(read, substr_p)

        return

    def apply(fastq_data, args):
        # TODO squash these into one arg, file_type

        if args.substr_p <= 0.0 or args.substr_p >= 1.0:
            raise Exception("Substring proportion must be in range (0.0, 1.0)")

        if args.interleaved or args.paired:
            return substring.apply_paired(fastq_data, args.substr_p)

        elif args.unpaired:
            return substring.apply_unpaired(fastq_data, args.substr_p)
