configfile: "config.json"
rule bowtie2_map:
    input:
       ref='{sample}.assembly.out/mc.refseq.fna',
       reads='{sample}.fasta'
    output:
       index= '{sample}.assembly.out/{sample}.index',
       pref='{sample}.assembly.out/{sample}.index',
       sam= '{sample}.assembly.out/{sample}.sam'
    log: '{sample}.step2.log'
    threads:config["nthreads"]
    message: """---Build index ."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -D 15 -R 2 -N 0 -L 31 -i S,1,1.15 -p {threads} -x {output.pref} -f {input.reads} -S {output.sam} > {log} 2>&1"
