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

#rule build_contigs:
#    input:
#        genome = '{sample}.assembly.out/mc.refseq.fna',
#        sam=  '{sample}.assembly.out/{sample}.sam'
#    output:
#        log='{sample}.assembly.out/{sample}.assembly.log',
#        out='{sample}.assembly.out',
#	contigs='{sample}.contigs.fasta'
#    message: """---Build contigs ."""
#    shell:'./buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 2 -l 300 -n T -b T -u T -k breadth  1>> {output.log} 2>&1;cp {output.out}/contigs.fasta {output.contigs}'


#rule align:
#    input:'{sample}.assembly.out/contigs.fasta'
#    output:
#       index= '{sample}.pilon/{sample}.index',
#       log= '{sample}.pilon/btw.log',
#       pref='{sample}.pilon/{sample}.index',
#       sam= '{sample}.pilon/{sample}.sam'
#    message: """---Build index to pilon."""
#    shell: "mkdir -p Sample1.pilon; bowtie2-build -q {input} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -k 30 -p 1 -x {output.pref} -f reads.fasta -S {output.sam} > {output.log} 2>&1"

#rule correct_contigs_pilon:
#    input:'{sample}.pilon/Sample1.sam'
#    output:
#        log='{sample}.samtoolssort.log',
#	bam='{sample}.sorted.bam'
#    message: """---Samtools sort ."""
#    shell:"samtools view -bS {input}| samtools sort - Sample1.sorted 1>> {output.log} 2>&1 ; touch {output.bam}"

#rule correct_contigs_pilon2:
#    input:'{sample}.sorted.bam'
#    output:
#        log='{sample}.samtoolsindex.log',
#        bam='{sample}.sorted.bam.bai'
#    message: """---Samtools index ."""
#    shell:"samtools index {input} 1>> {output.log} 2>&1"

