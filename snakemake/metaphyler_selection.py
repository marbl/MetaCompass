configfile: "config.json"
rule fastq2fasta:
    input: '{sample}.fastq'
    output:'{sample}.fasta'
    message: """---Converting fastq to fasta."""
    shell : "perl ./MetaCompass/bin/fq2fa.pl -i {input} -o {output}"


rule reference_recruitment:
    input:
        rules.fastq2fasta.output#'{sample}.fasta'
    output:
        out = '{sample}.assembly.out',
	reffile = '{sample}.assembly.out/mc.refseq.fna'
    log:'{sample}.step1.log'
    message: """---reference recruitment."""
    shell:"mkdir -p {output.out}; ./MetaCompass/bin/pickrefseqs.pl {input} {output.out} 1  1> {log} 2>&1"
