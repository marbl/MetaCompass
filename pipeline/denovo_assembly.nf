#!/usr/bin/env nextflow

/*
* Runs megahit on a set of reads
* expects either format -1 forward -2 reverse, or -r unpaired in the reads
* variable
*/

params.minlen=500  // default minimum contig size

process deNovoAssembly {
    publishDir {
        path: file("$params.output/denovo_assembly"),
        mode: 'copy'
    }

    when:
        params.de_novo == 1

    input:
        path reads

    output:
        path "megahit_out"
        val flag

    script:
    """
    megahit --min-count 3 --min-contig-len ${params.minlen} -t ${params.threads} --12 $reads
    flag=1
    """

}
