configfile: "config.json"
rule fastq2fasta:
    input: '{sample}.fastq'
    output:'{sample}.fasta'
    message: """---Converting fastq to fasta."""
    shell : "perl ./MetaCompass-beta/bin/fq2fa.pl -i {input} -o {output}"


rule reference_recruitment:
    input:
        fasta = '{sample}.fasta'
    output:
        out = '{sample}.assembly.out',
	log='{sample}.step1.log',
	reffile = '{sample}.assembly.out/mc.refseq.fna'
    message: """---reference recruitment."""
    shell:"mkdir -p {output.out}; ./MetaCompass-beta/bin/pickrefseqs.pl {input.fasta} {output.out} 1   1> {output.log} 2>&1"

