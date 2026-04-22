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
    path "concat_refs_mapped.fq"

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
    minimap2 -t ${params.threads} --heap-sort=yes -x sr \$batch \
        concat_refs.fna ${interleaved_reads} | cut -f 1 | sed -E 's#/(1|2)\$##' > aligned_reads.txt

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

// -----------------------------------------------------------------------------
// Run one reference-guided assembly task per cluster
// -----------------------------------------------------------------------------
process refAssemblyCluster {
    tag { cluster_list.baseName }
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    /*
     * Be conservative at first. Increase later if memory allows.
     */
    maxForks 4

    input:
    path cluster_list
    path reduced_reads
    path "*"
    path "*"

    output:
    path "*.refctgs.fna", emit: genomes

    script:
    """
    PYTHONPATH=${projectDir}/scripts
    export PYTHONPATH

    cluster_id=\$(basename "${cluster_list}" .csv | sed 's/^cluster_//')

    python -u -m align_reads \
        -ir ${reduced_reads} \
        -cl ${cluster_list} \
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
