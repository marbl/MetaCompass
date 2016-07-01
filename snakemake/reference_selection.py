configfile: "config.json"

rule kmer_mask:
    input: '{sample}.fastq'
    output:
        fastq='{sample}.marker.match.1.fastq'
    message: """---kmer-mask fastq"""
    params:'{sample}.marker'
    threads:config["nthreads"]
    shell : "kmer-mask -ms 28 -mdb ./MetaCompass/refseq/kmer-mask_db/markers.mdb -1 {input} -clean 0.0 -match 0.01 -nomasking -t {threads} -l 103 -o {params}"

rule reference_recruitment:
    input:
        rules.kmer_mask.output
    output:
        out = '{sample}.assembly.out',
	reffile = '{sample}.assembly.out/mc.refseq.fna'
    log:'{sample}.step1.log'
    message: """---reference recruitment."""
    threads:config["nthreads"]
    shell:"mkdir -p {output.out}; ./MetaCompass/bin/pickrefseqs.pl {input} {output.out} {threads}  1> {log} 2>&1"
