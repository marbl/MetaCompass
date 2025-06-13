#!/usr/bin/env nextflow

// finalize pipeline
// copy directory for each genome into output and just what we need
// copy assembly from de novo if needed
// combine all contigs into a single ouptut

process createOutputs {
    publishDir {
        path: file("${params.output}/output")
        mode: 'move'
    }

    input:
    val denovo    
    path "*"  

    output:
    path "output/*"

    script:
    """
    mkdir output
    cat *.refctgs.fna > output/all_contigs.fna
    cp -L *.refctgs.fna output
    if [ $denovo == 1 ] ;then
        cat $PWD/${params.output}/denovo_assembly/megahit_out/final.contigs.fa >> output/all_contigs.fna
        cp -L $PWD/${params.output}/denovo_assembly/megahit_out/final.contigs.fa output/denovo.contigs.fna
    fi
    """
     
}
