"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: "config.json"

ruleorder: merge_reads > bowtie2_map > build_contigs 

rule all:   
     input:expand('{prefix}/{sample}.{iter}.assembly.out/contigs.fasta',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])

rule merge_reads:
    input:
        reads=config['reads'].split(",")[0]
    message: """---merge fastq reads"""
    output:
        merged='%s/%s.merged.fq'%(config['prefix'],config['sample'])
    run:

        for read in config['reads'].split(','):
            if read != "" and len(read) != 0:
                os.system("ln -s %s  %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
                #if len(read) == 1:
                #    os.system("ln -s %s  %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
                #else:
                #    os.system("cat %s >> %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
        if config['reads'] != "" and config['reference'] != "%s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',sample=config['sample'],prefix=config["prefix"])[0]:
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
             os.system("touch %s"%expand('{prefix}/{sample}.{iter}.assembly.out/contigs.fasta',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
#if len(read) == 1:
                #    os.system("ln -s %s  %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
                #elif len(read) > 1:
                #    os.system("cat %s >> %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
rule bowtie2_map:
    input:
       ref=config['reference'],
       r1=rules.merge_reads.output.merged
    output:
       index=expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
    log: '%s/%s.%s.bowtie2map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Build index ."""
    shell:"bowtie2-build -o 3 --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --end-to-end --sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1} -S {output.sam}.all > {log} 2>&1;python3 %s/bin/best_strata.py {output.sam}.all {output.sam}"%(config["mcdir"])

rule build_contigs:
    input:
        genome = '%s'%(config['reference']),
        sam=  rules.bowtie2_map.output.sam
    params:
        pickref="%s"%(config['pickref']),
        mincov="%d"%(int(config['mincov'])),
        minlen="%d"%(int(config['minlen']))
    output:
        out='%s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']),
        contigs='%s/%s.%s.assembly.out/contigs.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.buildcontigs.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Build contigs ."""
    shell:"%s/bin/buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c {params.mincov} -l {params.minlen} -n T -b F -u F -k {params.pickref}  1>> {log} 2>&1"%(config["mcdir"])

