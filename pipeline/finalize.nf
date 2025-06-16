#!/usr/bin/env nextflow

// finalize pipeline
// copy directory for each genome into output and just what we need
// copy assembly from de novo if needed
// combine all contigs into a single ouptut

process createOutputs {
    publishDir {
        path: file("${workflow.outputDir}/output")
        mode: 'link'
    }

    input:
    val flag
    path "*"  

    output:
    path "output/*"

    script:
    """
    mkdir output
    cat *.refctgs.fna > output/all_contigs.fna
    cp -L *.refctgs.fna output
    if [ $flag == 1 ] ;then
       cat ${workflow.outputDir}/denovo_assembly/megahit_out/final.contigs.fa >> output/all_contigs.fna
       cp -L ${workflow.outputDir}/denovo_assembly/megahit_out/final.contigs.fa output/denovo.contigs.fna 
    fi
    """
}
/*
process createDeNovoOutputs {
    input:
    path "megahit_out"  

    script:
    """
    cp -L megahit_out/final.contigs.fa ${workflow.outputDir}/output/denovo.contigs.fna
    cat megathit_out/final.contigs.fa >> ${workflow.outputDir}/output/all_contigs.fna
    """
}
*/
