configfile: "config.json"

rule kmer_mask:
    input: '{prefix}/{sample}.fastq'
    output:
        fastq='{prefix}/{sample}.marker.match.1.fastq'
    message: """---kmer-mask fastq"""
    params:'{prefix}/{sample}.marker'
    threads:config["nthreads"]
    shell : "kmer-mask -ms 28 -mdb ../refseq/kmer-mask_db/markers.mdb -1 {input} -clean 0.0 -match 0.01 -nomasking -t {threads} -l 103 -o {params}"

rule fastq2fasta:
    input: rules.kmer_mask.output
    output:'{prefix}/{sample}.fasta'
    message: """---Converting fastq to fasta."""
    shell : "perl ../bin/fq2fa.pl -i {input} -o {output}"

rule reference_recruitment:
    input:
        rules.fastq2fasta.output
    output:
        out ='{prefix}/{sample}.assembly.out',
	reffile = '{prefix}/{sample}.assembly.out/mc.refseq.fna'
    log:'{prefix}/{sample}.step1.log'
    message: """---reference recruitment."""
    threads:config["nthreads"]
    shell:"mkdir -p {output.out}; ../bin/pickrefseqs.pl {input} {output.out} {threads}  1> {log} 2>&1"
