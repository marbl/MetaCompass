#!/usr/bin/env nextflow

/*
* Runs megahit on a set of reads
* expects either format -1 forward -2 reverse, or -r unpaired in the reads
* variable
*/
process deNovoAssembly {
    publishDir {
       path: file("params.output/denovo"),
       mode: 'copy'
    }

    when:
        params.denovo == 1

    input:
        val reads

    output:
        path "megahit_out"

    script:
    """
    megahit --min-count 3 --min-contig-len ${params.minlen} -t ${params.threads} $reads
    """
}
