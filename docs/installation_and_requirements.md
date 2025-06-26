## Requirements 

### Hardware requirements

* 90GB or more hard disk space to perform a normal installation.
* 8GB or more memory to allocate to the JVM (needed for pilon error correction step (https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage ). The amount of memory required depends on the genome, the read data, and how many fixes Pilon needs to make. Generally, bacterial genomes with ~200x of Illumina coverage will require at least 8GB, though 16GB is recommended.

### Software requirements

The file metacompass_environment.yml contains the list of packages needed to run MetaCompass. Currently all these packages are available through conda/bioconda.

## Installation

Follow these steps to get started with running the MetaCompass software:

### 1. Clone the MetaCompass Repository

    ```console
    git clone https://github.com/marbl/MetaCompass.git
    ```

### 2. Install Conda and Create Conda Environment

   - Install Conda, a package and environment manager, if you haven't already.
   Follow the official Conda installation [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
   
   - Create a Conda environment from the `metacompass_environment.yml` file provided in the cloned repository. Use the following command in your terminal:

     ```console
        conda env create -f metacompass_environment.yml
     ```

   - Metacompass  makes use of Nextflow, an open-source workflow framework for scientific 
   and data-intensive computing.<br/><br/>
      
      To ensure Nextflow is properly installed and configured on your system, follow these steps:

      1.  Activate the metacompass environment.
      2. Type `nextflow -v` in the command line. If everything is setup correctly you should see an output similar
      to this:

        ```console
        ()$ conda activate metacompass
        (metacompass)$ nextflow help -v 
        nextflow version 21.10.6.5660
        ```

### 3. Set up reference database.

A reference database is required to run Metacompass. You have the option to either generate a reference database on your own or utilize the one provided by us.

### Option A: Use Pre-built Database

1. Download the reference database (approximately 16GB):
   ```bash
   wget https://obj.umiacs.umd.edu/metacompass-db/RefSeq_V2_db.tar.gz
   ```

2. Extract the database:
   ```bash
   tar -xzf RefSeq_V2_db.tar.gz
   ```

3. Use in pipeline:
   ```bash
   nextflow run metacompass.nf 
    --reference_db /path/to/RefSeq_V2_db/RefSeq_V2_db 
    --forward "$forward_read" \ # [required]
    --reverse "$reverse_read" \ # [required]
    --output "$output_folder" \ # [required]
    --threads 8 \ # [optional] by default it is 16
    --trace_file_name "$output_folder/trace.txt" \ # [optional] 
    -with-timeline "$output_folder/timeline.html" \ # [optional]

   ```
   ```

### Option B: Generate the database yourself 

The reference database can also be set up manually with an input file formatted according to the NCBI RefSeq assembly_summary.txt file (documented at [link](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt)). An example is included in the repository: ref_db/RefSeq_V2_db/data/filtered.txt. Only the columns labeled "accession", "taxid", and "ftp" are currently used and the rest of the columns can be left blank.

The manual reference database setup has two options, depending on the size of the input file:

* If this input file is small (less than ~10,000 lines; estimate, will depend on computational resources):
  * Place your input text file in the following location: `/metacompass/ref_db/RefSeq_V2_db/data/`
  * Activate the metacompass conda environment
  * Identify the full local path of the directory "ref_db" within the cloned metacompass repository. e.g.:
    ```bash 
    ./metacompass/ref_db
    ```
  * Run the script "setup_ref_db_small.sh" located in this directory with the ref_db path as the first argument and the input file as the second argument:
    ```bash 
    ./setup_ref_db_small.sh ./metacompass/ref_db filtered.txt
    ```
  
  Note: This script will take a long time to complete

* If this input file is large (greater than ~10,000 lines; estimate, will depend on computational resources):
  * Place your input text file in the following location: `./metacompass/ref_db/RefSeq_V2_db/data/`
  * Identify the full local path of the directory "ref_db" within the cloned metacompass repository.  Set repository path (e.g., if repository is in the current directory): 
    ```bash
    repository_path=./metacompass/ref_db
    ```
  * Navigate to the data directory:
    ```bash
    cd ${repository_path}/RefSeq_V2_db/data
    ```
  * Split the accession text file:
    ```bash
    split -l 5000 ${repository_path}/RefSeq_V2_db/data/prokaryotes.txt ${repository_path}/RefSeq_V2_db/data/prokaryotes_
    ```
  * Activate the metacompass conda environment
  * Run submit_build_ref.sh with the ref_db as the first argument. This will submit several jobs that will process the split up input file in parallel:
    ```bash 
    cd $repository_path
    ./submit_build_ref.sh ${repository_path}
    ```
  * Once ALL jobs complete successfully, combine the outputs:
    ```bash
    cd $repository_path
    ./combineOutputs.sh ${repository_path}
    ```
  * (Optional) Remove split directories now that everything is combined:
    ```bash
    rm -r marker_index_[a-z][a-z]
    ```
  * Finally, run the script "setup_ref_db.sh" located in this directory with the ref_db path as the first argument, and the filename of the input text file as the second argument:
    ```bash 
    ./setup_ref_db.sh ./metacompass/ref_db prokaryotes.txt
    ```
