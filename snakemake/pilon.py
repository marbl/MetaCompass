#from os.path import join
import os#csv
#FILE ="5.fasta"#"Sample1.contigs.fasta"#,i 'Samplei2']
configfile: "config.json"
rule extract_single_fasta: 
    input:'{sample}.assembly.out/contigs.fasta'# config["references"]["genome"]#"wildcards.reference]
    output:'{sample}.genomes/NC_000915.1_0.fasta'
    message: """---single fasta genomes.""" 
#    shell:"echo"
    run:shell("echo ")
	references='Sample1.assembly.out/mc.refseq.fna'#config["references"]["genome"]
	f = open(references)
	shell("mkdir -p Sample1.genomes")
	contigs=[]
	for line in f:
            if line[0] == ">" :
                contigs=(line[1:].rstrip().split(' ',1)[0])	
		shell("samtools faidx {references} {contigs}> Sample1.genomes/{contigs}.fasta")#shell("mkdir -p Sample1.g;samtools faidx {FILE} {contigs}> Sample1.g/{contigs}.fasta")

rule correct_pilon:
    input:'{sample}.genomes/NC_000915.1_0.fasta'
    output:'{sample}.output.log'
    shell:"java -jar pilon --genome {input} --frags Sample1.pilon/Sample1.sorted.bam --output Sample1.pilon/NC_000915 --fix bases --tracks --vcf --changes  1> {output} 2>&1"

