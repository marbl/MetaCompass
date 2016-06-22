configfile: "config.json"
#sample= expand('{sample}.fastq', sample=config["reads"]["S1"]),#newref.fasta  selected_maps.sam
rule fastq2fasta:
    input: '{sample}.fastq'#[wildcards.unit]
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
