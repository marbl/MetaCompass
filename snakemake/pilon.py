import os
configfile: "config.json"
rule extract_single_fasta: 
    input:'{sample}.assembly.out/contigs.fasta'# config["references"]["genome"]#"wildcards.reference]
    output:'{sample}.genomes/NC_000915.1_0.fasta'
    message: """---single fasta genomes.""" 
#    shell:"echo"

    run:shell("echo ")
	references='Sample1.assembly.out/contigs.fasta'#config["references"]["genome"]
	f = open(references)
	shell("mkdir -p Sample1.genomes")
	contigs=[]
	for line in f:
            if line[0] == ">" :
                contigs=(line[1:].rstrip().split(' ',1)[0])	
		shell("samtools faidx {references} {contigs}> Sample1.genomes/{contigs}.fasta")#shell("mkdir -p Sample1.g;samtools faidx {FILE} {contigs}> Sample1.g/{contigs}.fasta")

#def _run_pilon(genome,bam,contig_id):
#    shell("java -jar pilon --genome {input} --frags Sample1.pilon/Sample1.sorted.bam --output Sample1.pilon/NC_000915 --fix bases --tracks --vcf --changes  1> {output} 2>&1")
import re
rule correct_pilon:
    input:'{sample}.genomes'#/NC_000915.1_0.fasta'
    output:'{sample}.pilon/NC_000915.1_0.fasta'#{sample}.output.log'
    run:x='Sample1'
    rootdir = os.getcwd()+'/'+'Sample1'+'.genomes'
    for root, dirs, filenames in os.walk(rootdir):
    	for f in filenames:    
	     name=re.sub(r'\.fasta$', '', f)
             shell("echo {f} >> Sample1.pilon/Sample1.output.log;java -jar pilon --genome {rootdir}/{f} --frags Sample1.pilon/Sample1.sorted.bam --output Sample1.pilon/{name} --fix bases --tracks --vcf --changes  1>> Sample1.pilon/Sample1.output.log 2>&1")
