"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json",

#include_prefix="https://"
#include_prefix + "rules"

#os.system("touch %s"%config['reads'][0])


ruleorder: merge_reads > kmer_mask > fastq2fasta > reference_recruitment > bowtie2_map > build_contigs > pilon_map > sam_to_bam > bam_sort > pilon_contigs > assemble_unmapped > join_contigs

#ruleorder: bowtie2_map > assemble_unmapped > build_contigs > pilon_map > sam_to_bam > bam_sort > pilon_contigs
#code to skip initial steps if reference genomes provided 
if config['reads'] != "" and config['reference'] != "%s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',sample=config['sample'],prefix=config["prefix"])[0]:
     os.system("touch %s"%expand('{prefix}/{sample}.marker.match.1.fastq',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("touch %s"%expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
     os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])[0])

     os.system("touch %s/%s.%s.assembly.out/%s.sam"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     os.system("touch %s/%s.%s.assembly.out/%s.sam.unmapped.fq"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     os.system('touch %s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']))
     os.system("touch %s/%s.%s.assembly.out/%s.index"%(config['prefix'],config['sample'],config['iter'],config['sample']))

     if int(config['iter']) >= 1:
         os.system('cp %s/%s.%d.assembly.out/contigs.final.fasta %s/%s.%s.assembly.out/contigs.fasta'%(config['prefix'],config['sample'],int(config['iter'])-1, config['prefix'],config['sample'],int(config['iter'])))


#rule all:
#    input:expand('{prefix}/{sample}.{iter}.assembly.out/contigs.fasta',sample=config["sample"],prefix=config["prefix"],iter=config["iter"])


#expand('{reads}',reads=config['r1'])[0]

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


rule kmer_mask:
    input:
        r1=rules.merge_reads.output.merged
    output:
        fastq1=expand('{prefix}/{sample}.marker.match.1.fastq',prefix=config['prefix'],sample=config['sample'])[0],
    message: """---kmer-mask fastq"""
    params:
        out=expand('{prefix}/{sample}.marker',prefix=config['prefix'],sample=config['sample'])[0],
        len=str(config["length"]+3)
    threads:config["nthreads"]
    log:'%s/%s.%s.kmermask.log'%(config['prefix'],config['sample'],config['iter'])
    shell : "kmer-mask -ms 28 -mdb %s/refseq/kmer-mask_db/markers.mdb -1 {input.r1} -clean 0.0 -match 0.01 -nomasking -t {threads} -l {params.len} -o {params.out} 1>> {log} 2>&1"%(config["mcdir"])


rule fastq2fasta:
    input: rules.kmer_mask.output.fastq1
    output:expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])
    message: """---Converting fastq to fasta."""
    shell : "perl %s/bin/fq2fa.pl -i {input} -o {output}"%(config["mcdir"])
    log:'%s/%s.%s.fastq2fasta.log'%(config['prefix'],config['sample'],config['iter'])

#    benchmark:
#       "%s/benchmarks/reference_recruitment/%s.txt"%(config['prefix'],config['sample'])

#fqfile =expand('{prefix}/{sample}.0.assembly.out/{sample}.fq',prefix=config['prefix'],sample=config['sample'])
#'{prefix}/{sample}.{iter}.assembly.out/mc.refseq.fna'
rule reference_recruitment:
    input:
        rules.fastq2fasta.output
    output:
        out =expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter']),
	reffile =expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])
    message: """---reference recruitment."""
    threads:config["nthreads"]
    log:'%s/%s.%s.reference_recruitement.log'%(config['prefix'],config['sample'],config['iter'])
    shell:"mkdir -p {output.out}; %s/bin/pickrefseqs.pl {input} {output.out} {threads}  1>> {log} 2>&1"%(config["mcdir"])


#reads=rules.reference_recruitment.output.fqfile
#    benchmark:
#       "%s/benchmarks/bowtie2_map/%s.txt"%(config['prefix'],config['sample'])

rule bowtie2_map:
    input:
       ref=rules.reference_recruitment.output.reffile,
       r1=rules.merge_reads.output.merged
    output:
       index= '%s/%s.%s.assembly.out/%s.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       pref='%s/%s.%s.assembly.out/%s.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),

    log: '%s/%s.%s.bowtie2map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:config["nthreads"]
    message: """---Build index ."""
    shell:"bowtie2-build -o 3 --threads {threads} -q %s {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --end-to-end --very-sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1}  -S {output.sam} > {log} 2>&1"%(config['reference'])



#shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --end-to-end -D 15 -R 2 -N 0 -L 31 -i S,1,1.15 -p {threads} -x {output.pref} -f {input.reads} -S {output.sam} > {log} 2>&1"

#    benchmark:
#       "%s/benchmarks/build_contigs/%s.txt"%(config['prefix'],config['sample'])
#
rule build_contigs:
    input:
        genome = '%s'%(config['reference']),
        sam=  rules.bowtie2_map.output.sam
    params:
        pickref="%s"%(config['pickref'])
    output:
        out='%s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']),
        contigs='%s/%s.%s.assembly.out/contigs.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.buildcontigs.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Build contigs ."""
    shell:"%s/bin/buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 5 -l 300 -n T -b F -u F -k {params.pickref}  1>> {log} 2>&1"%(config["mcdir"])

#/cbcb/project2-scratch/treangen/test78/c_rudii.0.assembly.out/c_rudii.megahit/final.contigs.fa 

rule pilon_map:
    input:
       ref=rules.build_contigs.output.contigs,
       r1=rules.merge_reads.output.merged
    output:
       index= '%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       pref='%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.mc.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       unmapped='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.fq'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.pilon.map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:config["nthreads"]
    message: """---Map reads for pilon polishing."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --very-sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1} -S {output.sam} --un {output.sam}.unmapped.fq > {log} 2>&1"


rule sam_to_bam:
    input:
        sam=  rules.pilon_map.output.sam
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
    threads:1
    message: """---Sort bam ."""
    shell: "samtools sort {input.bam} %s/%s.%s.assembly.out/sorted 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1"%(config['prefix'],config['sample'],config['iter'])


#samtools view -bS $samfile | samtools sort - $pilon_dir/$ref.sorted"
#$cmd = "samtools index $pilon_dir/$ref.sorted.bam";

rule pilon_contigs:
    input:
        contigs=rules.build_contigs.output.contigs,
        sam = rules.bam_sort.output.bam_sorted
    benchmark:
       "%s/benchmarks/pilon_contigs/%s.txt"%(config['prefix'],config['sample'])
    params:
        threads="%s"%(config['nthreads'])
    output:
        pilonctg='%s/%s.%s.assembly.out/contigs.pilon.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.pilon.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx16G -jar %s/bin/pilon-1.18.jar --flank 5 --threads {threads} --mindepth 4 --genome {input.contigs} --frags {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases 1>> {log} 2>&1"%(config["mcdir"],config['prefix'],config['sample'],config['iter'])
   
#concatenate this output with buildcontigs for pilon improvement
rule assemble_unmapped:
    input:
        reads=rules.pilon_map.output.unmapped
    output:
        megahit_contigs='%s/%s.0.assembly.out/%s.megahit/final.contigs.fa'%(config['prefix'],config['sample'],config['sample'])
    threads:config["nthreads"]
    log: '%s/%s.%s.megahit.log'%(config['prefix'],config['sample'],config['iter'])
    message: """---Assemble unmapped reads ."""
    shell:"rm -rf %s/%s.0.assembly.out/%s.megahit; megahit -o %s/%s.0.assembly.out/%s.megahit -t {threads} -r {input.reads}  1>> {log} 2>&1"%(config['prefix'],config['sample'],config['sample'],config['prefix'],config['sample'],config['sample'])


rule join_contigs:
    input:
        mc_contigs=rules.pilon_contigs.output.pilonctg,
        mh_contigs=rules.assemble_unmapped.output.megahit_contigs
    message: """---concanenate reference-guided and de novo contigs"""
    output:
        final_contigs="%s/%s.0.assembly.out/contigs.final.fasta"%(config['prefix'],config['sample'])
    shell:"cat {input.mh_contigs} {input.mc_contigs} > {output.final_contigs}"


#shell:"python %s/bin/fixhdr.py {input.contigs} ;java -Xmx16G -jar %s/bin/pilon-1.18.jar --threads {threads} --mindepth 0.75 --genome {input.contigs}.fna --frags {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases 1>> {log} 2>&1"%(config["mcdir"],config["mcdir"],config['prefix'],config['sample'],config['iter'])

