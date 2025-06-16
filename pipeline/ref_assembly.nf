#!/usr/bin/env nextflow

/*
* Interleave reads so that the next steps have an easier time
*/
process interleaveReads {
    publishDir {
        path: file("$workflow.outputDir/reference_assembly"),
        mode: 'copy'
    }

    input:
    val reads

    output:
    path "interleaved.fq"

    script:
    """
    reformat.sh -t=${params.threads} ${reads} out=interleaved.fq
    """
}

/*
*  Combines all clusters into one single FASTA file
*  aligns all reads to them and retains just the reads
*  that align to clusters
*/ 
process reduceClusters {
    publishDir {
        path: file("$workflow.outputDir/reference_assembly"),
        mode: 'copy'
    }
    input:
    path cluster_files
    path interleaved_reads

    output:
    path "concat_refs_mapped.fq"

    script:
    """
    python ${projectDir}/scripts/combine_clusters.py ${cluster_files} concat_refs.fna
   # get file size 
    let sz=`stat -c%s concat_refs.fna`/1024**3
    if [[ \$sz > 4 ]] ;then
       let batchs=sz+4
       batch="-I\${batchs}g"
    else
       batch=""
    fi
       
    minimap2 -t ${params.threads} --heap-sort=yes -x sr \$batch \
        concat_refs.fna ${interleaved_reads}| cut -d / -f 1 > aligned_reads.txt
    seqkit grep -I -j ${params.threads} -f aligned_reads.txt ${interleaved_reads} > concat_refs_mapped.fq
    """
}

/*
* 
*/
process refAssembly {
    publishDir {
        path: file("$workflow.outputDir/reference_assembly"),
        mode: 'copy'
    }
    
    input:
        path reads // input reads
        path cluster_list // clusters.txt from Cluster process
//        val readIndex // kmer index of reads
        path reduced_reads // just the reads that align to clusters
        path "*" // read index
        path "*" // cluster index

    output:
        path "unmapped.fq", emit: unmapped_reads
        path "*.refctgs.fna", emit: genomes


    script:
    """
    PYTHONPATH=${projectDir}/scripts
    export PYTHONPATH
    python -u -m align_reads \
        -ir ${reads} -cl ${cluster_list} \
        -rs . \
        -as reads.base_2.kmers  \
        -o . \
        -debug \
        -mcl 2000 \
        -t ${params.threads}
    python ${projectDir}/scripts/collate_genomes.py

    """
}
