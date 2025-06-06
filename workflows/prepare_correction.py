configfile: "config.json"

rule indexing_bowtie2:
    input:
       ref='{sample}.assembly.out/contigs.fasta'
    output:#index='{sample}.pilon/{sample}.index',
       outdir='{sample}.pilon',
       fwd='{sample}.pilon/{sample}.index.1.bt2'
#       fwd=expand('{sample}.pilon/{sample}.index.{index}.bt2', sample=config["reads"]["S1"], index=range(1,5) ),
#       rev=expand("{sample}.pilon/{sample}.index.rev.{index}.bt2",sample=config["reads"]["S1"], index=range(1,3))
    log:'{sample}.pilon/bowtie-build.log'
    threads:config["nthreads"]
    params:'{sample}.pilon/{sample}.index'
    message: """---Build index to pilon."""
    shell: "mkdir -p {output.outdir}; bowtie2-build -q {input.ref} {params} 1>> {log} 2>&1"

rule mapping_bowtie2:
    input:
       reads='{sample}.fasta',
       fwd='{sample}.pilon/{sample}.index.1.bt2'
#=expand('{reads}.pilon/{reads}.index.{index}.bt2', reads=config["reads"]["S1"], index=range(1,5) ),
#       rev=expand("{reads}.pilon/{reads}.index.rev.{index}.bt2",reads=config["reads"]["S1"], index=range(1,3))
       #iindex=rules.indexing_bowtie2.output.index
    output:'{sample}.pilon/{sample}.sam'
    log:'{sample}.pilon/bowtie2.log'
    threads:config["nthreads"]
    params:"{sample}.pilon/{sample}.index"
    message: """---Build index to pilon."""
    shell: "bowtie2 --end-to-end -k 30 -p {threads} -x {params} -f {input.reads} -S {output} > {log} 2>&1"

rule samtools_sort:
    input:'{sample}.pilon/{sample}.sam'
    output:'{sample}.pilon/{sample}.sorted.bam'
    log:'{sample}.pilon/{sample}.samtoolssort.log'
    message: """---Samtools sort ."""
    threads:config["nthreads"]
    params:'{sample}.pilon/{sample}.sorted'
    shell:"samtools view -bS {input}| samtools sort - {params} 1>> {log} 2>&1"# ; touch {output.bam}"

rule samtools_index:
    input:'{sample}.pilon/{sample}.sorted.bam'
    output:'{sample}.pilon/{sample}.sorted.bam.bai'
    log:'{sample}.pilon/{sample}.samtoolsindex.log'
    message: """---Samtools index ."""
    threads:config["nthreads"]
    shell:"samtools index {input} 1>> {log} 2>&1"

