# MetaCompass v2.0-beta

## Publication

Victoria Cepeda, Bo Liu, Mathieu Almeida, Christopher M. Hill, Sergey Koren, Todd J. Treangen, Mihai Pop.
bioRxiv 212506; doi: https://doi.org/10.1101/212506

## Installation

1. Check if you have all the [requirements](docs/installation_and_requirements.md#requirements-).
2. Clone the repository
    ```shell
    git clone https://gitlab.umiacs.umd.edu/mpop/metacompass.git
    ```

3. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) and
   create a conda environment using the `metacompass_environment.yml` file

    ```shell
    conda env create -f metacompass_environment.yml
    ```

4. Setup reference database.

    You have the option to either generate a reference database on your own or utilize the one provided by us. Please refer to the instructions [here](docs/installation_and_requirements.md#3-setup-reference-database) for guidance on how to proceed.


Detailed installation instructions are [here](docs/installation_and_requirements.md).

## Usage 

This provides a quick overview of running the Metacompass.

1. Set up your input data and reference database paths.

2. Create an output directory.

3. Run Metacompass Nextflow script using appropriate parameters.

```bash
    
    NXF_OPTS="-Dleveldb.mmap=false -Xmx500g" nextflow run metacompass.nf \
    --reference_db "${ref_db_path}" \ # [required]
    --forward "$forward_read" \ # [required]
    --reverse "$reverse_read" \ # [required]
    --output "$output_folder" \ # [required]
    --threads 8 \ # [optional] by default it is 16
    --trace_file_name "$output_folder/trace.txt" \ # [optional] 
    -with-timeline "$output_folder/timeline.html" \ # [optional]
```

#### Parameter description: <br/>

- **forward** : Path to the file with forward reads, providing sequence information from one end of the DNA fragment in
  a specific direction.
- **reverse** : Path to the file with reverse reads, representing sequence information from the opposite end of the DNA
  fragment in the reverse direction.
- **output** : The output directory or file path where the results will be stored.
- **threads**: The number of threads to use during reference-guided assembly with Pilon.
  <br/> <br/>
  Full parameter configuration can be found [here](./docs/parameter_configuration.md)

Detailed usage instructions are [here](docs/usage_guide.md)