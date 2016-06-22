configfile: "config.json"

rule build_contigs:
    input:
        genome = '{sample}.assembly.out/mc.refseq.fna',
        sam=  '{sample}.assembly.out/{sample}.sam'
    output:
        log='{sample}.assembly.out/{sample}.assembly.log',
        out='{sample}.assembly.out',
	contigs='{sample}.contigs.fasta'
    message: """---Build contigs ."""
    shell:'./buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 2 -l 300 -n T -b T -u T -k breadth  1>> {output.log} 2>&1;cp {output.out}/contigs.fasta {output.contigs}'


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

