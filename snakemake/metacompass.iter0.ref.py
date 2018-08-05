"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json")

#include_prefix="https://"
#include_prefix + "rules"

#os.system("touch %s"%config['reads'][0])


ruleorder: merge_reads > bowtie2_map > build_contigs > pilon_map > sam_to_bam > bam_sort > pilon_contigs > assemble_unmapped > join_contigs > create_tsv

rule all:
     input:expand('{prefix}/metacompass.tsv',prefix=config["prefix"])

rule merge_reads:
    input:
        reads=config['reads'].split(",")[0]
    message: """---merge fastq reads"""
    output:
        merged='%s/%s.merged.fq'%(config['prefix'],config['sample'])
    run:

        for read in config['reads'].split(','):
            if read != "" and len(read) != 0:
                os.system("cat %s >> %s/%s.merged.fq"%(read,config['prefix'],config['sample']))
        if config['reads'] != "" and config['reference'] != "%s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',sample=config['sample'],prefix=config["prefix"])[0]:
             os.system("touch %s"%expand('{prefix}/{sample}.marker.match.1.fastq',prefix=config['prefix'],sample=config['sample'])[0])
             os.system("touch %s"%expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])[0])
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
             os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])[0])

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
    #shell:"bowtie2-build -o 3 --threads {threads} -q %s {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1}  -S {output.sam} > {log} 2>&1"%(config['reference'])
    shell:"bowtie2-build -o 3 --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --end-to-end --sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1} -S {output.sam}.all > {log} 2>&1; %s/bin/best_strata.py {output.sam}.all {output.sam}; rm {output.sam}.all"%(config["mcdir"])

rule build_contigs:
    input:
        genome = '%s'%(config['reference']),
        sam=  rules.bowtie2_map.output.sam
    params:
        pickref="%s"%(config['pickref']),
        mincov="%d"%(int(config['mincov']))
    output:
        out='%s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']),
        contigs='%s/%s.%s.assembly.out/contigs.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.buildcontigs.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Build contigs ."""
    shell:"%s/bin/buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c {params.mincov} -l 300 -n T -b F -u F -k {params.pickref}  1>> {log} 2>&1"%(config["mcdir"])

#/cbcb/project2-scratch/treangen/test78/c_rudii.0.assembly.out/c_rudii.megahit/final.contigs.fa 

rule pilon_map:
    input:
       ref=rules.build_contigs.output.contigs,
       #r1=rules.merge_reads.output.merged
       r1=config['reads'].split(",")[0],
       r2=config['reads'].split(",")[1]
    output:
       index=expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.mc.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.mc.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       #unmapped='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.fq'%(config['prefix'],config['sample'],config['iter'],config['sample'])
       unmappedr1='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.1.fq'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       unmappedr2='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.2.fq'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.pilon.map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Map reads for pilon polishing."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --no-mixed --sensitive --no-unal -p {threads} -x {output.pref} -q -1 {input.r1} -2 {input.r2} -S {output.sam} --un-conc {output.sam}.unmapped.fq > {log} 2>&1"


rule sam_to_bam:
    input:
        sam=rules.pilon_map.output.sam
    output:
        bam = "%s.bam"%(rules.pilon_map.output.sam)
    log:'%s/%s.%s.assembly.out/%s.assembly.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Convert sam to bam ."""
    shell:"samtools view -bS {input.sam} -o {output.bam} 1>> {log} 2>&1"

rule bam_sort:
    input:
        bam = rules.sam_to_bam.output.bam 
    output:
        bam_sorted = "%s/%s.%s.assembly.out/sorted.bam"%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.assembly.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Sort bam ."""
    shell: "samtools sort -@ {threads} {input.bam} -o %s/%s.%s.assembly.out/sorted.bam -O bam -T tmp 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1"%(config['prefix'],config['sample'],config['iter'])


#samtools view -bS $samfile | samtools sort - $pilon_dir/$ref.sorted"
#$cmd = "samtools index $pilon_dir/$ref.sorted.bam";

rule pilon_contigs:
    input:
        contigs=rules.build_contigs.output.contigs,
        sam = rules.bam_sort.output.bam_sorted
    #benchmark:
    #   "%s/benchmarks/pilon_contigs/%s.txt"%(config['prefix'],config['sample'])
    output:
        pilonctg='%s/%s.%s.assembly.out/contigs.pilon.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.pilon.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx19G -jar %s/bin/pilon-1.22.jar --flank 5 --threads {threads} --mindepth 0.5 --genome {input.contigs} --frags {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases 1>> {log} 2>&1"%(config["mcdir"],config['prefix'],config['sample'],config['iter'])
   
#rule remove_zerocov:
#    input:
#        contigs=rules.pilon_contigs.output.pilonctg,
#        bam=rules.bam_sort.output.bam_sorted
#    output:
#        filtctg='%s/%s.%s.assembly.out/contigs.pilon.fasta.fixed'%(config['prefix'],config['sample'],config['iter'])
#    log:'%s/%s.%s.assembly.out/%s.remove_zerocov.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
#    threads:int(config['nthreads'])
#    message: """---Remove zero coverage regions."""
#    shell:"python %s/bin/cut_zeros.py {input.contigs} {input.bam} 1>> {log} 2>&1"%(config["mcdir"])
#concatenate this output with buildcontigs for pilon improvement
#reads=rules.pilon_map.output.unmapped

rule assemble_unmapped:
    input:
        r1=rules.pilon_map.output.unmappedr1,
        r2=rules.pilon_map.output.unmappedr2
    output:
        megahit_contigs='%s/%s.0.assembly.out/%s.megahit/final.contigs.fa'%(config['prefix'],config['sample'],config['sample'])
    threads:int(config["nthreads"])
    log: '%s/%s.%s.megahit.log'%(config['prefix'],config['sample'],config['iter'])
    message: """---Assemble unmapped reads ."""
    shell:"if [[ -s %s/%s.0.assembly.out/%s.mc.sam.unmapped.1.fq || -s %s/%s.0.assembly.out/%s.mc.sam.unmapped.2.fq ]]; then rm -rf %s/%s.0.assembly.out/%s.megahit; megahit -o %s/%s.0.assembly.out/%s.megahit --min-count 3 --min-contig-len %d --presets meta-sensitive -t {threads} -1 {input.r1} -2 {input.r2}  1>> {log} 2>&1; else touch {output.megahit_contigs} {log}; echo 'No unmapped reads to run de novo assembly' >{log} ;fi"%(config['prefix'],config['sample'],config['sample'],config['prefix'],config['sample'],config['sample'],config['prefix'],config['sample'],config['sample'],config['prefix'],config['sample'],config['sample'],int(config['minlen']))
    #then touch {output.megahit_contigs} {log}; else
    #shell:"rm -rf %s/%s.0.assembly.out/%s.megahit; if  [[ -s tutorial_thao_mash/thao2000.0.assembly.out/thao2000.mc.sam.unmapped.1.fq ]]; then empty=0; else echo empty=1; fi; if [[ $empty ]]; then touch {output.megahit_contigs} {log}; fi;  megahit -o %s/%s.0.assembly.out/%s.megahit --min-count 3 --min-contig-len %d --presets meta-sensitive -t {threads} -1 {input.r1} -2 {input.r2}  1>> {log} 2>&1"%(config['prefix'],config['sample'],config['sample'],config['prefix'],config['sample'],config['sample'],int(config['minlen']))

#        mc_contigs=rules.remove_zerocov.output.filtctg,
rule join_contigs:
    input:
        mc_contigs=rules.pilon_contigs.output.pilonctg,
        mh_contigs=rules.assemble_unmapped.output.megahit_contigs
    message: """---concanenate reference-guided and de novo contigs"""
    output:
        final_contigs="%s/%s.0.assembly.out/contigs.final.fasta"%(config['prefix'],config['sample'])
    shell:"cat {input.mh_contigs} {input.mc_contigs} > {output.final_contigs}"

#shell:"python %s/bin/fixhdr.py {input.contigs} ;java -Xmx16G -jar %s/bin/pilon-1.18.jar --threads {threads} --mindepth 0.75 --genome {input.contigs}.fna --frags {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases 1>> {log} 2>&1"%(config["mcdir"],config["mcdir"],config['prefix'],config['sample'],config['iter'])

rule create_tsv:
    input:
        all_contigs=rules.join_contigs.output.final_contigs,
        mc_contigs=rules.build_contigs.output.contigs
    message: """---information reference-guided and de novo contigs"""
    output:"%s/metacompass.tsv"%(config['prefix'])
    shell:"sh %s/bin/create_tsv.sh {input.all_contigs} {input.mc_contigs} NA {output}"%(config["mcdir"])


