"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json",

#include_prefix="https://"
#include_prefix + "rules"

#os.system("touch %s"%config['reads'][0])


ruleorder: pilon_map > sam_to_bam > bam_sort > pilon_contigs 
#if config['reads'] != "" and config['reference'] != "%s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',sample=config['sample'],prefix=config["prefix"])[0]:
     #os.system("touch %s"%expand('{prefix}/{sample}.marker.match.1.fastq',prefix=config['prefix'],sample=config['sample'])[0])
     #os.system("touch %s"%expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])[0])
     #os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
     #os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
     #os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])[0])

     #os.system("touch %s/%s.%s.assembly.out/%s.sam"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     #os.system("touch %s/%s.%s.assembly.out/%s.sam.unmapped.fq"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     #os.system('touch %s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']))
     #os.system("touch %s/%s.%s.assembly.out/%s.index"%(config['prefix'],config['sample'],config['iter'],config['sample']))


rule pilon_map:
    input:
       ref='%s/%s.%d.assembly.out/contigs.final.fasta'%(config['prefix'],config['sample'],int(config['iter'])-1),
       r1='%s/%s.merged.fq'%(config['prefix'],config['sample'])
    output:
       index= expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.mc.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.mc.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       unmapped='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.fq'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.pilon.map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Map reads for pilon polishing."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --very-sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1} -S {output.sam} --un {output.sam}.unmapped.fq > {log} 2>&1"


rule sam_to_bam:
    input:
        sam=  rules.pilon_map.output.sam
    output:
        bam = "%s.bam"%(rules.pilon_map.output.sam)
    log:'%s/%s.%s.assembly.out/%s.assembly.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config["nthreads"])
    message: """---Convert sam to bam ."""
    shell:"samtools view -bS {input.sam} -o {output.bam} 1>> {log} 2>&1"

rule bam_sort:
    input:
        bam = rules.sam_to_bam.output.bam 
    output:
        bam_sorted = "%s/%s.%s.assembly.out/sorted.bam"%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.assembly.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config["nthreads"])
    message: """---Sort bam ."""
    shell: "samtools sort -@ {threads} {input.bam} -o %s/%s.%s.assembly.out/sorted.bam -O bam -T tmp 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1"%(config['prefix'],config['sample'],config['iter'])


#samtools view -bS $samfile | samtools sort - $pilon_dir/$ref.sorted"
#$cmd = "samtools index $pilon_dir/$ref.sorted.bam";

rule pilon_contigs:
    input:
        contigs='%s/%s.%d.assembly.out/contigs.final.fasta'%(config['prefix'],config['sample'],int(config['iter'])-1),
        sam = rules.bam_sort.output.bam_sorted
    benchmark:
       "%s/benchmarks/pilon_contigs/%s.txt"%(config['prefix'],config['sample'])
    params:
        threads="%s"%(int(config['nthreads'])),
        mincov="%d"%(int(config['mincov']))
    output:
        pilonctg='%s/%s.%s.assembly.out/contigs.final.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.pilon.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config["nthreads"])
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx64G -jar %s/bin/pilon-1.18.jar --flank 5 --threads {threads} --mindepth {params.mincov} --genome {input.contigs} --unpaired {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases,local 1>> {log} 2>&1;cp %s/%s.%s.assembly.out/contigs.pilon.fasta %s/%s.%s.assembly.out/contigs.final.fasta"%(config["mcdir"],config['prefix'],config['sample'],config['iter'],config['prefix'],config['sample'],config['iter'],config['prefix'],config['sample'],config['iter'])
