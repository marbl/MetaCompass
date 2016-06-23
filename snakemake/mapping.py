configfile: "config.json"
rule build_index:
    input:
        reffile = '{sample}.assembly.out/mc.refseq.fna'
    output:
       index= '{sample}.assembly.out/{sample}.index',
       log= '{sample}.step2.log',
       pref='{sample}.assembly.out/{sample}.index',
       sam= '{sample}.assembly.out/{sample}.sam'
    message: """---Build index ."""
    shell:"bowtie2-build -q {input.reffile} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -k 30 -p 1 -x {output.pref} -f reads.fasta -S {output.sam} > {output.log} 2>&1"
