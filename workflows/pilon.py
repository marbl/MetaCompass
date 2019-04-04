import os
import re
configfile: "config.json"
rule extract_single_fasta: 
    input:'{sample}.assembly.out/contigs.fasta'
    output:'{sample}.genomes/NC_000915.1_0.fasta'
    message: """---single fasta genomes.""" 
    run:
        references='Sample1.assembly.out/contigs.fasta'
	f = open(references)
	shell("mkdir -p Sample1.genomes")
	contigs=[]
	for line in f:
            if line[0] == ">" :
                contigs=(line[1:].rstrip().split(' ',1)[0])	
		shell("samtools faidx {references} {contigs}> Sample1.genomes/{contigs}.fasta")
rule correct_pilon:
    input:'{sample}.genomes'
    output:'{sample}.pilon/{sample}.output.log'
    run:
        rootdir = os.getcwd()+'/'+'Sample1'+'.genomes'
    	for root, dirs, filenames in os.walk(rootdir):
            for f in filenames:
                name=re.sub(r'\.fasta$', '', f)
                shell("echo {f} >> Sample1.pilon/Sample1.output.log;java -jar pilon --genome {rootdir}/{f} --frags Sample1.pilon/Sample1.sorted.bam --output Sample1.pilon/{name} --fix bases --tracks --vcf --changes  1>> Sample1.pilon/Sample1.output.log 2>&1")
