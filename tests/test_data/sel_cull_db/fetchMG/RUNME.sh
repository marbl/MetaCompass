# Print help text
echo -e "\nfetchMG retrieves sequences homologous to 40 universal marker genes"
echo -e "from a fasta formatted file with protein coding sequences.\n"

# Run fethMG on metagenome examples
echo -e "Running fetchMG on metagenomic genes"
./fetchMG.pl -m extraction -o example_output_metagenome example_datasets/example_data.faa -t 8

# Run fethMG on genome examples
echo -e "Running fetchMG on genomes"
./fetchMG.pl -m extraction -v -o example_output_genome example_datasets/example_data_genomes.faa -t 8

# Print help text
echo -e "============================================================================="
echo -e "\nfetchMG has finished. The sequences of the 40 marker genes"
echo -e "are now saved in the ./example_output_metagenome and ./example_output_genome"
echo -e "============================================================================="
