configfile: "config.json"
rule mapping_bowtie2:
    input:'{sample}.assembly.out/contigs.fasta'
    output:
       index= '{sample}.pilon/{sample}.index',
       log= '{sample}.pilon/btw.log',
       pref='{sample}.pilon/{sample}.index',
       sam= '{sample}.pilon/{sample}.sam'
    message: """---Build index to pilon."""
    shell: "mkdir -p Sample1.pilon; bowtie2-build -q {input} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -k 30 -p 1 -x {output.pref} -f reads.fasta -S {output.sam} > {output.log} 2>&1"

rule samtools_sort:
    input:'{sample}.pilon/{sample}.sam'
    output:
        log='{sample}.pilon/{sample}.samtoolssort.log',
	bam='{sample}.pilon/{sample}.sorted.bam'
    message: """---Samtools sort ."""
    shell:"samtools view -bS {input}| samtools sort - Sample1.pilon/Sample1.sorted 1>> {output.log} 2>&1"# ; touch {output.bam}"

rule samtools_index:
    input:'{sample}.pilon/{sample}.sorted.bam'
    output:
        log='{sample}.pilon/{sample}.samtoolsindex.log',
        bam='{sample}.pilon/{sample}.sorted.bam.bai'
    message: """---Samtools index ."""
    shell:"samtools index {input} 1>> {output.log} 2>&1"

