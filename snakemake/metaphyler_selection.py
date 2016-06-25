configfile: "config.json"
rule fastq2fasta:
    input: '{sample}.fastq'
    output:'{sample}.fasta'
    message: """---Converting fastq to fasta."""
<<<<<<< HEAD
    shell : "perl ./MetaCompass/bin/fq2fa.pl -i {input} -o {output}"
=======
    shell : "perl ./MetaCompass-beta/bin/fq2fa.pl -i {input} -o {output}"
>>>>>>> c3e0b58d0fbb423c53fd9ec8ec39f7b468aa9a6c


rule reference_recruitment:
    input:
        rules.fastq2fasta.output#'{sample}.fasta'
    output:
        out = '{sample}.assembly.out',
	reffile = '{sample}.assembly.out/mc.refseq.fna'
<<<<<<< HEAD
    log:'{sample}.step1.log'
    message: """---reference recruitment."""
    shell:"mkdir -p {output.out}; ./MetaCompass/bin/pickrefseqs.pl {input} {output.out} 1  1> {log} 2>&1"
=======
    log:
	'{sample}.step1.log'
    message: """---reference recruitment."""
    shell:"mkdir -p {output.out}; ./MetaCompass-beta/bin/pickrefseqs.pl {input} {output.out} 1   1> {log} 2>&1"

>>>>>>> c3e0b58d0fbb423c53fd9ec8ec39f7b468aa9a6c
