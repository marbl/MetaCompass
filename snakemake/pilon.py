#from os.path import join
#import csv
#FILE ="5.fasta"#"Sample1.contigs.fasta"#,i 'Samplei2']
configfile: "config.json"
rule extract_single_fasta: 
    #input: config["references"]["genome"]#"wildcards.reference]
    message: """---single fasta genomes."""
    run:shell("echo ")
	#references={input}
	references=config["references"]["genome"]
	f = open(references)#open(FILE, 'r')
	shell("mkdir -p Sample1.g")
	contigs=[]
	for line in f:
            if line[0] == ">" :
                contigs=(line[1:].rstrip().split(' ',1)[0])	
		shell("samtools faidx {references} {contigs}> Sample1.g/{contigs}.fasta")#shell("mkdir -p Sample1.g;samtools faidx {FILE} {contigs}> Sample1.g/{contigs}.fasta")

rule correct_pilon:
    input:'{sample}.sorted.bam'
    output:
        log='{sample}.piloncorr.log'
    message:"""---Running pilon."""#memory optional
    shell:"java -jar pilon --genome Sample1.genomes/NC_000915_0.fasta --frags {input} --output Sample1.genomes/NC_000915_0. --fix bases --tracks --vcf --changes  1> {output.log} 2>&1"
