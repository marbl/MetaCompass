configfile: "config.json"
rule build_index:
    input:
       '{sample}.assembly.out/mc.refseq.fna'
    output:
       index= '{sample}.assembly.out/{sample}.index',
       pref='{sample}.assembly.out/{sample}.index',
       sam= '{sample}.assembly.out/{sample}.sam'
    log: '{sample}.step2.log'
    threads:1
    message: """---Build index ."""
    shell:"bowtie2-build -q {input} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -k 30 -p 1 -x {output.pref} -f reads.fasta -S {output.sam} > {log} 2>&1"
