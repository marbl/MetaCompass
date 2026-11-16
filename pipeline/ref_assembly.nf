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
// process interleaveReads {
//     publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

//     input:
//     val reads

//     output:
//     path "interleaved.fq"

//     script:
//     """
//     reformat.sh -t=1 ${reads} out=interleaved.fq
//     """
// }

// -----------------------------------------------------------------------------
// Combine cluster FASTA files, align reads to the combined references,
// and retain only reads that map to at least one cluster
// -----------------------------------------------------------------------------
process reduceClusters {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    path cluster_files
    path reads_r1
    path reads_r2

    output:
    path "concat_refs_mapped_R1.fq", emit: mapped_r1
    path "concat_refs_mapped_R2.fq", emit: mapped_r2
    path "concat_refs_unmapped_R1.fq", emit: unmapped_r1
    path "concat_refs_unmapped_R2.fq", emit: unmapped_r2
    path "read_to_genome.tsv", emit: read_to_genome

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

    touch read_to_genome.tsv
    minimap2 -t ${params.threads} --heap-sort=yes -x sr \$batch \
        concat_refs.fna ${reads_r1} ${reads_r2} | \
        awk 'BEGIN{OFS="\\t"} {
            base=\$1;
            sub(/\\/[12]\$/, "", base);

            print base;
            print base "/1";
            print base "/2";

            print base, \$6 >> "read_to_genome.tsv";
        }' > aligned_reads.txt

    seqkit grep -I -j ${params.threads} \
        -f aligned_reads.txt ${reads_r1} > concat_refs_mapped_R1.fq

    seqkit grep -I -j ${params.threads} \
        -f aligned_reads.txt ${reads_r2} > concat_refs_mapped_R2.fq

    seqkit grep -I -j ${params.threads} -v \
        -f aligned_reads.txt ${reads_r1} > concat_refs_unmapped_R1.fq

    seqkit grep -I -j ${params.threads} -v \
        -f aligned_reads.txt ${reads_r2} > concat_refs_unmapped_R2.fq
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
    path mapped_r1
    path mapped_r2
    path ref_to_cluster

    output:
    path "cluster_reads/cluster_*_R1.fq", emit: reads_r1
    path "cluster_reads/cluster_*_R2.fq", emit: reads_r2

    script:
    """
    mkdir -p cluster_reads

    python ${projectDir}/scripts/build_cluster_read_subsets.py \
        --cluster-list ${cluster_files} \
        --fq1 ${mapped_r1} \
        --fq2 ${mapped_r2} \
        --ref-map ${ref_to_cluster} \
        --outdir cluster_reads

    mkdir -p cluster_reads/paired

    for r1 in cluster_reads/cluster_*_R1.fq; do
        r2="\${r1/_R1.fq/_R2.fq}"

        seqkit pair \
            --id-regexp '^([^/[:space:]]+)' \
            -1 "\$r1" \
            -2 "\$r2" \
            -O cluster_reads/paired \
            -f

        mv "cluster_reads/paired/\$(basename "\$r1")" "\$r1"
        mv "cluster_reads/paired/\$(basename "\$r2")" "\$r2"
    done

    rmdir cluster_reads/paired
    """
}


// -----------------------------------------------------------------------------
// Split each cluster into:
//   1. relatively contiguous genomes for Pangraph
//   2. fragmented genomes to be mapped after Pangraph residual blocks
// -----------------------------------------------------------------------------
process preparePangraphRefs {

    tag { cluster_list.baseName }

    publishDir "${REFERENCE_ASSEMBLY_DIR}/pangraph", mode: 'copy'

    input:
    path cluster_list

    output:
    path "${cluster_list.baseName}.pangraph.fna", emit: pangraph_fasta
    path "${cluster_list.baseName}.pangraph_refs.csv", emit: pangraph_refs
    path "${cluster_list.baseName}.fragmented_refs.csv", emit: fragmented_refs

    script:
    """
    python - <<'PY'
    import csv
    import shutil
    from pathlib import Path

    cluster_file = Path("${cluster_list}")
    cluster_name = cluster_file.stem

    max_contigs = ${params.pangraph_max_contigs}

    with cluster_file.open() as fh:
        rows = list(csv.reader(fh))

    if not rows:
        raise ValueError(f"Empty cluster file: {cluster_file}")

    refs = [
        Path(x)
        for x in rows[0]
        if x.strip()
    ]


    def count_contigs(fasta):
        count = 0

        with fasta.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1

        return count


    stats = []

    for ref in refs:

        if not ref.exists():
            raise FileNotFoundError(ref)

        n_contigs = count_contigs(ref)

        stats.append(
            (ref, n_contigs)
        )


    pangraph_refs = [
        ref
        for ref, n_contigs in stats
        if n_contigs <= max_contigs
    ]

    fragmented_refs = [
        ref
        for ref, n_contigs in stats
        if n_contigs > max_contigs
    ]

    # ---------------------------------------------------------
    # Write list of genomes participating in Pangraph
    # ---------------------------------------------------------

    with open(
        f"{cluster_name}.pangraph_refs.csv",
        "w",
        newline=""
    ) as fh:

        writer = csv.writer(fh)

        writer.writerow(
            [str(x) for x in pangraph_refs]
        )


    # ---------------------------------------------------------
    # Write fragmented genomes for later Level-2 mapping
    # ---------------------------------------------------------

    with open(
        f"{cluster_name}.fragmented_refs.csv",
        "w",
        newline=""
    ) as fh:

        writer = csv.writer(fh)

        writer.writerow(
            [str(x) for x in fragmented_refs]
        )


    # ---------------------------------------------------------
    # Concatenate ONLY Pangraph-eligible genome FASTAs
    #
    # IMPORTANT:
    # preserve original FASTA headers.
    # This allows read_aligner.py to match backbone contig IDs
    # directly against Pangraph path IDs.
    # ---------------------------------------------------------

    pangraph_fasta = Path(
        f"{cluster_name}.pangraph.fna"
    )

    with pangraph_fasta.open("wb") as out:

        for ref in pangraph_refs:

            with ref.open("rb") as src:

                shutil.copyfileobj(
                    src,
                    out
                )


    print(
        f"{cluster_name}: "
        f"total={len(refs)}, "
        f"pangraph={len(pangraph_refs)}, "
        f"fragmented={len(fragmented_refs)}"
    )

    for ref, n in sorted(
        stats,
        key=lambda x: x[1]
    ):

        category = (
            "PANGRAPH"
            if ref in pangraph_refs
            else "LEVEL2"
        )

        print(
            f"{category}\\t"
            f"{n}\\t"
            f"{ref}"
        )
    PY
        """
    }


// MOUMI: added for generating pangraph
// -----------------------------------------------------------------------------
// Build one pangenome graph per reference cluster
// -----------------------------------------------------------------------------
process buildPangraph {

    tag { pangraph_fasta.baseName }

    publishDir "${REFERENCE_ASSEMBLY_DIR}/pangraph", mode: 'copy'

    cpus 4

    input:
    path pangraph_fasta

    output:
    tuple val("${pangraph_fasta.baseName.replace('.pangraph', '')}"),
          path("*_graph.json"),
          emit: graphs

    script:
    
    def clusterName = pangraph_fasta.baseName.replace(".pangraph", "")

    """
    num_sequences=\$(awk '/^>/{n++} END{print n+0}' ${pangraph_fasta})

    echo "${clusterName}: Pangraph input contains \$num_sequences sequence(s)"

    if [[ \$num_sequences -ge 2 ]]; then

        pangraph \
            --verbosity info \
            build \
            -j ${task.cpus} \
            ${pangraph_fasta} \
            -o ${clusterName}_graph.json

    else

        echo "${clusterName}: fewer than 2 sequences; skipping Pangraph."
        touch ${clusterName}_graph.json

    fi
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
    tuple val(cluster_name), path(cluster_list), path(cluster_r1), path(cluster_r2), path(pangraph_graph), path(pangraph_refs), path(fragmented_refs)
    path "*"
    path "*"

    output:
    path "*.refctgs.fna", emit: genomes
    tuple path("${cluster_name}.unmapped_R1.fq"),
        path("${cluster_name}.unmapped_R2.fq"),
        emit: unmapped_reads
    tuple path("${cluster_name}.mapped_R1.fq"),
        path("${cluster_name}.mapped_R2.fq"),
        emit: mapped_reads

    script:
    """
    PYTHONPATH=${projectDir}/scripts
    export PYTHONPATH

    # MOUMI: change made for sending subset reads to corresponding cluster from initial global mapping

    cluster_id=\$(basename "${cluster_list}" .csv | sed 's/^cluster_//')
    
    python -u -m align_reads \
        -f ${cluster_r1} \
        -r ${cluster_r2} \
        -cl ${cluster_list} \
        -rs . \
        -as reads.base_2.kmers \
        -o . \
        -debug \
        -mcl 2000 \
        -t ${params.threads} \
        --single-cluster \
        --cluster-id \${cluster_id} \
        --pangraph-graph ${pangraph_graph} \
        --pangraph-refs ${pangraph_refs} \
        --fragmented-refs ${fragmented_refs}
    
    for mate in R1 R2; do
        if [[ -f unmapped_\${mate}.fq ]]; then
            mv unmapped_\${mate}.fq "${cluster_name}.unmapped_\${mate}.fq"
        else
            touch "${cluster_name}.unmapped_\${mate}.fq"
        fi

        if [[ -f cluster_mapped_\${mate}.fq ]]; then
            mv cluster_mapped_\${mate}.fq "${cluster_name}.mapped_\${mate}.fq"
        else
            touch "${cluster_name}.mapped_\${mate}.fq"
        fi
    done

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

process combineClusterUnmappedReads {
    publishDir "${REFERENCE_ASSEMBLY_DIR}", mode: 'copy'

    input:
    path original_r1
    path original_r2
    path cluster_mapped_reads

    output:
    tuple path("unmapped_R1.fq"),
          path("unmapped_R2.fq"),
          emit: unmapped_reads

    script:
    """
    cat ${cluster_mapped_reads.join(' ')} | \
        seqkit seq -n -i | \
        sed -E 's#/[12]\$##' | \
        sort -u > mapped_ids.txt

    awk '{
        print;
        print \$0 "/1";
        print \$0 "/2"
    }' mapped_ids.txt > exclude_ids.txt

    seqkit grep -v -f exclude_ids.txt \
        ${original_r1} > unmapped_R1.fq

    seqkit grep -v -f exclude_ids.txt \
        ${original_r2} > unmapped_R2.fq
    """
}