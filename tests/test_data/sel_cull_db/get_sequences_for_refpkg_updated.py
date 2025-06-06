import shutil, os
import sys
import re
from itertools import groupby
from pathlib import Path


def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fh.close()


import json
import csv


def detectContaminated(input_file, output_file):
    # Open the JSONL file for reading and the CSV file for writing
    with open(input_file, 'r') as jsonl_file, open(output_file, 'a', newline='') as csv_file:
        # Create a CSV writer object
        csv_writer = csv.writer(csv_file)
        # Check if the CSV file is empty. If it's empty, write the header row.
        if csv_file.tell() == 0:
            csv_writer.writerow(["assemblyAccession", "organismName", "isAtypical", "hasWarning", "warnings"])
        # Process each line in the JSONL file
        for line in jsonl_file:
            # Load the JSON data from the current line
            data = json.loads(line)
            # Extract the required data
            assemblyAccession = data['assemblyInfo']['assemblyAccession']
            # Extracting organismName can be tricky since it can be nested under different keys
            # We'll use a try-except block to handle potential KeyErrors
            try:
                organismName = data['description']['organism']['organismName']
            except KeyError:
                organismName = data.get('organismName', '')

            # Extract the 'atypical' information
            atypical_info = data['assemblyInfo'].get('atypical', {})

            # Get the 'isAtypical' value, default to False if it's not present or blank
            isAtypical = atypical_info.get('isAtypical', False)

            # # Convert 'isAtypical' to a boolean if it's a string
            # if isinstance(isAtypical, str):
            # 	isAtypical = isAtypical.lower() == "true"

            # isAtypical = data['assemblyInfo'].get('atypical', {}).get('isAtypical', '')
            warnings = ';'.join(data['assemblyInfo'].get('atypical', {}).get('warnings', []))
            hasWarning = bool(warnings.strip())
            # Write the extracted data to the CSV file
            csv_writer.writerow([assemblyAccession, organismName, isAtypical, hasWarning, warnings])
            return isAtypical, hasWarning
    print(f"Data has been successfully parsed and saved to {output_file}.")


def main():
    # wrkdir = "/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2/RefSeq_V2_db/"
    wrkdir = sys.argv[1] + "/RefSeq_V2_db/"
    print("Working directory set as: " + wrkdir)

    # version = 'TEST'
    if len(sys.argv) == 3:
        version = ''
    else:
        version = sys.argv[3]
    cognames = {}
    with open(wrkdir + 'data/FetchMG_COGid2gene_name.tsv') as f:
        for num, line in enumerate(f):
            if num == 0:
                continue
            val = line.strip().split('\t')
            cognames[val[0]] = val[1]
    filehandle_cogs_fna = {}
    filehandle_cogs_faa = {}
    filehandle_cogs_genome = {}
    for gene in cognames:
        os.system("mkdir -p " + wrkdir + 'marker_index' + version + '/' + cognames[gene] + '_' + gene)
        filehandle_cogs_fna[gene] = open(
            wrkdir + 'marker_index' + version + '/' + cognames[gene] + '_' + gene + '/' + cognames[
                gene] + '_' + gene + '.fna', 'w')
        filehandle_cogs_faa[gene] = open(
            wrkdir + 'marker_index' + version + '/' + cognames[gene] + '_' + gene + '/' + cognames[
                gene] + '_' + gene + '.faa', 'w')
        filehandle_cogs_genome[gene] = open(
            wrkdir + 'marker_index' + version + '/' + cognames[gene] + '_' + gene + '/' + cognames[
                gene] + '_' + gene + '_genome.tsv', 'w')
        filehandle_cogs_genome[gene].write('seq_id\tmarker_gene\taccession\tlen_gene_nuc\n')

    # file_accession=[wrkdir + 'data/prokaryotes_'+ version + '.txt']
    # file_accession=[wrkdir + 'data/prokaryotes.txt']
    file_accession = [sys.argv[2]]
    # file_accession = '/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2/RefSeq_V2_db/data/prokaryotes_unique_new.txt'
    # This file is from https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
    for bac_arc_ind, refseq_file in enumerate(file_accession):
        with open(refseq_file) as f:
            for num, line in enumerate(f):
                print("Parsing line number " + str(num))
                val = line.strip().split('\t')
                if line.startswith('#'):
                    continue
                accession = val[18]
                taxid = val[1]
                ftp = val[20]
                assemblyname = re.sub("^.*" + accession + "_", "", ftp)

                try:
                    # datasets_command = "/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2/datasets download genome accession " + accession + " --exclude-gff3 --exclude-rna --filename "+ wrkdir+ "marker_index_"+version+"/ncbi_dataset.zip"
                    datasets_command = "datasets download genome accession " + accession + " --exclude-gff3 --exclude-rna --filename " + wrkdir + "marker_index" + version + "/ncbi_dataset.zip"
                    os.system(datasets_command)
                    os.system(
                        "unzip " + wrkdir + "marker_index" + version + "/ncbi_dataset.zip -d " + wrkdir + "marker_index" + version)
                except Exception as e:
                    print("Download genome data FAILED for ", accession)
                    continue

                trailname = wrkdir + "marker_index" + version + "/ncbi_dataset/data/"
                # Look at genome sequence (nucleotide), protein sequences (aa), and cds sequences (nucleotide) for each genome accession
                accession_folder = Path(trailname + accession + "/")
                accession_files = [file.name for file in accession_folder.iterdir() if file.is_file()]
                genomefilename = trailname + accession + "/" + accession + "_" + assemblyname + "_genomic.fna"
                for fname in accession_files:
                    if accession in fname and ".fna" in fname:
                        genomefilename = trailname + accession + "/" + fname

                cdsfilename = trailname + accession + "/cds_from_genomic.fna"
                proteinfilename = trailname + accession + "/protein.faa"

                assemblyfilename = trailname + "/assembly_data_report.jsonl"

                try:
                    f1 = open(genomefilename)
                    f1.close()
                except Exception as e:
                    print("Missing genome file from ", accession)
                    print(genomefilename)
                    continue

                try:
                    f2 = open(cdsfilename)
                    f2.close()
                except Exception as e:
                    print("Missing cds file from ", accession)
                    continue

                try:
                    f3 = open(proteinfilename)
                    f3.close()
                except Exception as e:
                    print("Missing protein file from ", accession)
                    continue

                try:
                    f1 = open(genomefilename)
                    f1.close()
                except Exception as e:
                    print("Missing genome file from ", accession)
                    os.system("rm " + cdsfilename)
                    os.system("rm " + proteinfilename)
                    continue

                try:
                    f2 = open(cdsfilename)
                    f2.close()
                except Exception as e:
                    print("Missing cds file from ", accession)
                    os.system("rm " + genomefilename)
                    os.system("rm " + proteinfilename)
                    continue

                try:
                    f3 = open(proteinfilename)
                    f3.close()
                except Exception as e:
                    print("Missing protein file from ", accession)
                    os.system("rm " + genomefilename)
                    os.system("rm " + cdsfilename)
                    continue

                assembly_output_stats = wrkdir + "parsed_assembly_data_report" + version + ".csv"
                isAtypical, hasWarning = detectContaminated(assemblyfilename, assembly_output_stats)
                if isAtypical or hasWarning:
                    print("Record: " + accession + " is contaminated (is atypical or has warning) and was skipped.")
                    continue

                proteinseqmap = {}
                fiter = fasta_iter(proteinfilename)
                for ff in fiter:
                    header = ff[0].strip().split(' ')[0]
                    proteinseqmap[header] = ff[1]

                gene_nuc_filename = accession + "_gene.fna"
                gene_aa_filename = accession + "_gene.faa"
                report_filename = accession + "_gene_summary.txt"
                fgn = open(gene_nuc_filename, 'w')
                fga = open(gene_aa_filename, 'w')
                fgr = open(report_filename, 'w')

                # Keep only those genes that have protein id and are not pseudo genes
                fiternuc = fasta_iter(cdsfilename)
                for ff in fiternuc:
                    header = ff[0]
                    seq_id = accession + '_' + header.split(' ')[0].split('|')[-1]
                    nucseq = ff[1]
                    if "protein_id=" in header:
                        protein_name = header.split('protein_id=')[1].split(']')[0]
                    elif "pseudo=true" not in header:
                        print("#", accession, header, "\tcds with no protein assigned and not pseudo gene")
                        continue
                    else:
                        continue
                    if protein_name in proteinseqmap:
                        corr_protein = proteinseqmap[protein_name]
                    else:
                        print("#", accession, header, "\tprotein missing")
                        continue

                    nuclen = len(nucseq)
                    aalen = len(corr_protein)

                    if nuclen - 3 == 3 * aalen:
                        matchflag = 1
                    else:
                        matchflag = 0
                        print("#", accession, header, "\tnucleotide and aa seq length don't match\t", nuclen, aalen,
                              str(nuclen - 3 - (aalen * 3)), str(bac_arc_ind))
                    ## Writes nucleotide gene sequences, aa gene sequences and the metadata information for all genes for the accession. Later, we will keep records for only those genes that belong to 40 marker genes.
                    fgr.write('\t'.join(
                        [seq_id, accession, protein_name, str(nuclen), str(aalen), str(matchflag), taxid, header,
                         str(bac_arc_ind)]) + '\n')
                    fgn.write('>' + seq_id + '\n' + nucseq + '\n')
                    fga.write('>' + seq_id + '\n' + corr_protein + '\n')
                fgr.close()
                fgn.close()
                fga.close()

                ## Fetchmg retrieves all genes that align well to one of the 40 marker genes. -v parameter retrieves very best hit for each marker gene.
                output_folder = wrkdir + "marker_index" + version + "/" + accession + "_cogoutput"
                try:
                    fetchmg_command = sys.argv[
                                          1] + "/fetchMG/fetchMG.pl -m extraction -v " + gene_aa_filename + " -o " + output_folder
                    os.system(fetchmg_command)

                except Exception as e:
                    print("FetchMG failed", accession)
                    continue

                ## write to the final files
                genesofinterest = {}
                for gene in cognames:
                    ## reading nuc sequence first
                    if os.path.exists(output_folder + '/' + gene + '.fna'):
                        fitergene = fasta_iter(output_folder + '/' + gene + '.fna')
                        for ff in fitergene:
                            genesofinterest[ff[0]] = gene
                            filehandle_cogs_fna[gene].write('>' + ff[0] + '\n' + ff[1] + '\n')

                    ## reading aa sequence second
                    if os.path.exists(output_folder + '/' + gene + '.faa'):
                        fitergene = fasta_iter(output_folder + '/' + gene + '.faa')
                        for ff in fitergene:
                            filehandle_cogs_faa[gene].write('>' + ff[0] + '\n' + ff[1] + '\n')

                ## write relevant enteries from report file
                with open(report_filename) as f:
                    for line in f:
                        val = line.strip().split('\t')
                        if val[0] in genesofinterest:
                            gene = genesofinterest[val[0]]
                            filehandle_cogs_genome[gene].write(
                                val[0] + '\t' + cognames[gene] + '_' + gene + '\t' + val[1] + '\t' + val[3] + '\n')

                ##cleanup unnecessary files
                try:
                    os.system("rm " + gene_nuc_filename + "*")
                    os.system("rm " + gene_aa_filename + "*")
                    os.system("rm " + report_filename)
                    os.system("rm " + wrkdir + "marker_index" + version + "/README.md")
                    os.system("rm " + wrkdir + "marker_index" + version + "/ncbi_dataset.zip")
                    shutil.rmtree(output_folder)
                # shutil.rmtree(wrkdir+"marker_index"+version+"/ncbi_dataset")
                except:
                    print("Error while cleaning up", accession)
                    continue

                print("#completed", accession)


if __name__ == '__main__':
    main()
