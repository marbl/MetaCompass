"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json",
ruleorder: pilon_contigs

rule all:
     input:
         sixth=expand('{outdir}/error_correction/contigs.pilon.fasta',outdir=config['outdir'])
   
rule pilon_contigs:
    input:
        contigs='%s/assembly/contigs.fasta'%(config['outdir']),
        sam = "%s/error_correction/sorted.bam"%(config['outdir'])
    output:
        pilonctg='%s/error_correction/contigs.pilon.fasta'%(config['outdir'])
    params:
        memory="%d"%(int(config['memory'])),
        retry='%s/error_correction/run4.ok'%(config['outdir'])
    log:'%s/logs/pilon.log'%(config['outdir'])
    threads:int(config['nthreads'])
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx{params.memory}G -jar %s/bin/pilon-1.23.jar --flank 5 --threads {threads} --mindepth 3 --genome {input.contigs} --frags {input.sam} --output %s/error_correction/contigs.pilon --fix bases,amb --tracks --changes --vcf 1>> {log} 2>&1  && touch {params.retry}"%(config["mcdir"],config['outdir'])
