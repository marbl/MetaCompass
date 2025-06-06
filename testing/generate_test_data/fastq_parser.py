# TODO optimize I/O speed if it's a bottleneck
# TODO add output for interleaved files if needed. Right now it outputs two files or one unpaited file


class fastq_parser:
    def __init__(self, file_type):
        self.file_type = file_type

    def parse(self, fastq_files):
        if self.file_type == "paired":
            return fastq_parser.parse_paired(fastq_files)

        elif self.file_type == "unpaired":
            return fastq_parser.parse_unpaired(fastq_files)

        elif self.file_type == "interleaved":
            return fastq_parser.parse_interleaved(fastq_files)

        else:
            return Exception("invalid_fq_file_type")

    def write(self, fastq_data, fastq_out_files):
        if self.file_type == "paired" or self.file_type == "interleaved":
            return fastq_parser.write_paired_reads(fastq_data, fastq_out_files)

        if self.file_type == "unpaired":
            return fastq_parser.write_unpaired_reads(fastq_data, fastq_out_files[0])

        else:
            return Exception("invalid_fq_file_type")

    def parse_unpaired(fastq_filename):
        # print("parsing unpaired fastq file: " + str(fastq_filename))
        fastq_file = open(fastq_filename, "r")

        fastq_data = []
        cur_read = ["", "", "", ""]
        ctr = 0

        for line in fastq_file:
            i = ctr
            cur_read[i] = line

            if i == 3:
                fastq_data.append(cur_read.copy())

            ctr = ctr + 1
            ctr = ctr % 4

        fastq_file.close

        last_line = fastq_data[len(fastq_data) - 1][3]
        if last_line[len(last_line) - 1] != "\n":
            fastq_data[len(fastq_data) - 1][3] = last_line + "\n"

        return fastq_data

    def parse_paired(fastq_filename_pair):

        fastq_file_2 = open(fastq_filename_pair[1], "r")
        fastq_file_1 = open(fastq_filename_pair[0], "r")

        fastq_data = []
        cur_read = ["", "", "", ""]
        ctr = 0

        for line in fastq_file_1:
            i = ctr % 4
            cur_read[i] = line
            if i == 3:
                fastq_data.append([cur_read.copy()])

            ctr = ctr + 1

        cur_read = ["", "", "", ""]
        ctr = 0

        for line in fastq_file_2:
            i = ctr % 4
            cur_read[i] = line
            if i == 3:
                # print(str(ctr // 4) + " " + str(len(fastq_data)))
                fastq_data[ctr // 4].append(cur_read.copy())

            ctr = ctr + 1

        last_line = fastq_data[len(fastq_data) - 1][0][3]
        if last_line[len(last_line) - 1] != "\n":
            fastq_data[len(fastq_data) - 1][0][3] = last_line + "\n"

        last_line = fastq_data[len(fastq_data) - 1][1][3]
        if last_line[len(last_line) - 1] != "\n":
            fastq_data[len(fastq_data) - 1][1][3] = last_line + "\n"

        fastq_file_1.close()
        fastq_file_2.close()

        return fastq_data

    def parse_interleaved(fastq_filename):
        #print("parsing interleaved fastq file: " + str(fastq_filename))

        fastq_file = open(fastq_filename, "r")
        if not fastq_file:
            return fastq_parser.file_not_found(fastq_filename)

        fastq_data = []
        cur_read_pair = [["", "", "", ""], ["", "", "", ""]]
        ctr = 0

        for line in fastq_file:
            i = ctr // 4
            j = ctr % 4
            cur_read_pair[i][j] = line

            if i == 1 and j == 3:
                fastq_data.append([cur_read_pair[0].copy(), cur_read_pair[1].copy()])

            ctr = ctr + 1
            ctr = ctr % 8

        last_line = fastq_data[len(fastq_data) - 1][1][3]
        if last_line[len(last_line) - 1] != "\n":
            fastq_data[len(fastq_data) - 1][1][3] = last_line + "\n"

        fastq_file.close()

        return fastq_data

    def write_paired_reads(fastq_data, out_filenames):
        out_file_1 = open(out_filenames[0], "w")
        out_file_2 = open(out_filenames[1], "w")

        for read_pair in fastq_data:
            out_file_1.write(
                "{0}{1}{2}{3}".format(
                    read_pair[0][0], read_pair[0][1], read_pair[0][2], read_pair[0][3]
                )
            )

        for read_pair in fastq_data:
            out_file_2.write(
                "{0}{1}{2}{3}".format(
                    read_pair[1][0], read_pair[1][1], read_pair[1][2], read_pair[1][3]
                )
            )

        out_file_1.close()
        out_file_2.close()
        return

    def write_unpaired_reads(fastq_data, out_filename):
        out_file = open(out_filename, "w")

        for read in fastq_data:
            out_file.write("{0}{1}{2}{3}".format(read[0], read[1], read[2], read[3]))

        out_file.close()
        return
