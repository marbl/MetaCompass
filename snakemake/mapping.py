configfile: "config.json"
<<<<<<< HEAD
rule bowtie2_map:
=======
rule build_index:
>>>>>>> c3e0b58d0fbb423c53fd9ec8ec39f7b468aa9a6c
    input:
       '{sample}.assembly.out/mc.refseq.fna'
    output:
       index= '{sample}.assembly.out/{sample}.index',
       pref='{sample}.assembly.out/{sample}.index',
       sam= '{sample}.assembly.out/{sample}.sam'
    log: '{sample}.step2.log'
    threads:1
    message: """---Build index ."""
<<<<<<< HEAD
    shell:"bowtie2-build -q {input} {output.pref} 1>> {output.index} 2>&1;bowtie2 -L 31 -p {threads} -x {output.pref} -f Sample1.fasta -S {output.sam} > {log} 2>&1"
=======
    shell:"bowtie2-build -q {input} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -k 30 -p 1 -x {output.pref} -f reads.fasta -S {output.sam} > {log} 2>&1"
>>>>>>> c3e0b58d0fbb423c53fd9ec8ec39f7b468aa9a6c
