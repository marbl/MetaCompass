#!/usr/bin/env nextflow

// -----------------------------------------------------------------------------
// Default parameters
// -----------------------------------------------------------------------------
params.output = params.output ?: ''

// -----------------------------------------------------------------------------
// Helper values
// -----------------------------------------------------------------------------
final String REFERENCE_ASSEMBLY_DIR = params.output ?
    "${params.output}/reference_assembly" :
    "reference_assembly"

// -----------------------------------------------------------------------------
// Interleave reads so downstream steps can consume a single FASTQ
// -----------------------------------------------------------------------------
process interleaveReads {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    val reads

    output:
    path "interleaved.fq"

    script:
    """
    reformat.sh -t=${params.threads} ${reads} out=interleaved.fq
    """
}

// -----------------------------------------------------------------------------
// Combine cluster FASTA files, align reads to the combined references,
// and retain only reads that map to at least one cluster
// -----------------------------------------------------------------------------
process reduceClusters {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    path cluster_files
    path interleaved_reads

    output:
    path "concat_refs_mapped.fq", emit: mapped_reads
    path "read_to_genome.tsv", emit: read_to_genome
    path "concat_refs_mapped.sam", emit: global_sam

    script:
    """
    python ${projectDir}/scripts/combine_clusters.py ${cluster_files} concat_refs.fna

    # Get file size in GB.
    # The +1 is needed because the exit code is 1 for arithmetic operations
    # that yield a result of 0.
    size_bytes=\$(stat -c%s concat_refs.fna)
    let sz=size_bytes/1024**3+1

    if [[ \$sz > 5 ]]; then
        let batchs=sz+4
        batch="-I\${batchs}g"
    else
        batch=""
    fi

    # Added to remove extra /1 and /2 suffixes from read names that minimap2
    # may produce when aligning paired-end reads
    minimap2 -a -t ${params.threads} --heap-sort=yes -x sr \$batch \
        concat_refs.fna ${interleaved_reads} > concat_refs_mapped.sam

    awk 'BEGIN{OFS="\t"} !/^@/ {
        read=\$1;
        flag=\$2;
        gsub(/\\/1\\/1\$/, "/1", read);
        gsub(/\\/2\\/2\$/, "/2", read);
        ref=\$3;

        if (read !~ /\\/[12]\$/) {
            if (and(flag, 64)) {
                read = read "/1";
            } else if (and(flag, 128)) {
                read = read "/2";
            }
        }

        if (ref != "*") {
            print read >> "aligned_reads.txt";
            print read, ref >> "read_to_genome.tsv";
        }
    }' concat_refs_mapped.sam

    seqkit grep -I -j ${params.threads} -f aligned_reads.txt ${interleaved_reads} > concat_refs_mapped.fq
    """
}


// MOUMI: added for parallelization
// -----------------------------------------------------------------------------
// Split clusters.csv into one cluster file per row
// -----------------------------------------------------------------------------
process splitClusters {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    path cluster_files

    output:
    path "split_clusters/cluster_*.csv", emit: cluster_refs

    script:
    """
    mkdir -p split_clusters

    python - <<'PY'
import csv
from pathlib import Path

cluster_file = Path("${cluster_files}")
outdir = Path("split_clusters")
outdir.mkdir(exist_ok=True)

with cluster_file.open() as fh:
    rows = list(csv.reader(fh))

for i, row in enumerate(rows, start=1):
    out = outdir / f"cluster_{i}.csv"
    with out.open("w", newline="") as ofh:
        writer = csv.writer(ofh)
        writer.writerow(row)
PY
    """
}


// MOUMI: added for building cluster read subsets for parallelization
// -----------------------------------------------------------------------------
// Build FASTQ files containing only reads that map to each cluster.
// -----------------------------------------------------------------------------
process buildClusterReadSubsets {
    input:
    path cluster_files
    path mapped_fq
    path ref_to_cluster
    path global_sam

    output:
    path "cluster_reads/cluster_*.fq", emit: cluster_reads
    path "cluster_maps/cluster_*.read_to_genome.tsv", emit: cluster_maps
    path "cluster_sams/cluster_*.sam", emit: cluster_sams
    
    script:
    """
    mkdir -p cluster_reads cluster_maps cluster_sams

    python ${projectDir}/scripts/build_cluster_read_subsets.py \
        --cluster-list ${cluster_files} \
        --fq ${mapped_fq} \
        --ref-map ${ref_to_cluster} \
        --outdir cluster_reads \
        --map-outdir cluster_maps \
        --sam ${global_sam} \
        --sam-outdir cluster_sams
    """
}

// MOUMI: added for parallelization
// -----------------------------------------------------------------------------
// Run one reference-guided assembly task per cluster
// -----------------------------------------------------------------------------
process refAssemblyCluster {
    tag { cluster_name }
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    maxForks 4

    // MOUMI: change made for sending subset reads to corresponding cluster from initial global mapping 
    input:
    // path cluster_list
    // path reduced_reads
    tuple val(cluster_name), path(cluster_list), path(cluster_reads), path(cluster_map), path(cluster_sam)
    path "*"
    path "*"

    output:
    path "*.refctgs.fna", emit: genomes

    script:
    """
    PYTHONPATH=${projectDir}/scripts
    export PYTHONPATH

    # MOUMI: change made for sending subset reads to corresponding cluster from initial global mapping

    cluster_id=\$(basename "${cluster_list}" .csv | sed 's/^cluster_//')
    
    python -u -m align_reads \
        -ir ${cluster_reads} \
        -cl ${cluster_list} \
        --read-to-genome ${cluster_map} \
        --global-sam ${cluster_sam} \
        -rs . \
        -as reads.base_2.kmers \
        -o . \
        -debug \
        -mcl 2000 \
        -t ${params.threads} \
        --single-cluster \
        --cluster-id \${cluster_id}

    python ${projectDir}/scripts/collate_genomes.py
    """
}

// -----------------------------------------------------------------------------
// Perform reference-guided assembly
// -----------------------------------------------------------------------------
process refAssembly {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    path reads
    path cluster_list
    path reduced_reads
    path "*"
    path "*"

    output:
    path "unmapped.fq", emit: unmapped_reads
    path "*.refctgs.fna", emit: genomes

    script:
    """
    PYTHONPATH=${projectDir}/scripts
    export PYTHONPATH

    python -u -m align_reads \
        -ir ${reads} \
        -cl ${cluster_list} \
        -rs . \
        -as reads.base_2.kmers \
        -o . \
        -debug \
        -mcl 2000 \
        -t ${params.threads}

    python ${projectDir}/scripts/collate_genomes.py
    """
}
