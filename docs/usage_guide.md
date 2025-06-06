# Running Metacompass

This guide provides instructions on how to run the Metacompass.

## Requirements

Before running the Metacompass, ensure that you have the following:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (setting it up using our conda environment is recommended)
- Reference database (`ref_db_path`) (see build instructions [here](./installation_and_requirements.md#3-setup-reference-database))
- Input forward and reverse read files in FASTQ format
- Other [requirements](./installation_and_requirements.md#requirements-)

## Usage

To run Metacompass effectively, we recommend creating a script that incorporates the following key elements:

1. **Input Data Paths**: Specify the paths to the input forward and reverse read files.

2. **Parameter Configuration**: Define various parameters, including:
   - The reference database
   - The output directory
   - The number of threads to use
   

   ([Full parameter configuration](./parameter_configuration.md))
   

Creating such a script will simplify the execution of Metacompass for your specific analysis.

## Example: 

#### example.sh
   ```bash
   # Set the paths to your input data and reference database by modifying
   # the following variables in your shell script:

   forward_read="forwad_read.fastq.gz" # will also work with .fastq file
   reverse_read="reverse_read.fastq.gz" # will also work with .fastq file
   ref_db_path="/path/to/your/reference/database"
   output_folder="/path/to/your/output/directory"
   
   # Run metacompass on these variables using the following command:
   
   nextflow run metacompass2.nf \
    --reference_db "${ref_db_path}" \ # [required]
    --forward "$forward_read" \ # [required]
    --reverse "$reverse_read" \ # [required]
    --output "$output_folder" \ # [required]
    --threads 8 \ # [optional] by default it is 16
    --trace_file_name "$output_folder/trace.txt" \ # [optional] 
    -with-timeline "$output_folder/timeline.html" \ # [optional]
    -with-dag "$output_folder/${read}_dag.png" # [optional]
    
    # --trace_file_name: Path to a nextflow trace file
    # -with-timeline: Generates a timeline HTML report. [optional]
    # -with-dag: Generates a Directed Acyclic Graph (DAG) visualization. [optional]
  ```

```shell
./example.sh
```

Monitor the progress and view the results in the specified output directory.

### Toggling *de novo* assembly on and off
By default, MetaCompass also generates a *de novo* assembly using the reads that were not used by the reference-guided step. To turn off this feature, you can use the parameter ```--de_novo 0``` . 
