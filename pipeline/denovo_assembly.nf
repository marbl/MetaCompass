#!/usr/bin/env nextflow

/*
* Runs megahit on a set of reads
* expects either format -1 forward -2 reverse, or -r unpaired in the reads
* variable
*/

params.minlen=500  // default minimum contig size

process deNovoAssembly {
    publishDir {
        path: file("$workflow.outputDir/denovo_assembly"),
        mode: 'copy'
    }

    input:
        path reads

    output:
        val 1, emit: flag
        path "megahit_out"

    script:
    """
    megahit --min-count 3 --min-contig-len ${params.minlen} -t ${params.threads} --12 $reads
    """

}
