"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: "config.json"

ruleorder: merge_reads > pilon_map > sam_to_bam > bam_sort > pilon_contigs > assemble_unmapped > join_contigs > create_tsv

rule all:   
     input:expand('{prefix}/metacompass_summary.tsv',prefix=config["prefix"])

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
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
             os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
             os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])[0])

rule pilon_map:
    input:
       ref=config['assembly'],
       r1=config['r1'],#['reads'].split(",")[0],
       r2=config['r2'],#['reads'].split(",")[1],
       ru=config['ru']#['reads'].split(",")[2],
    output:
       index=expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.mc.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.mc.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam2='%s/%s.%s.assembly.out/%s.mc_unpaired.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       unmappedr1='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.1.fq'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       unmappedr2='%s/%s.%s.assembly.out/%s.mc.sam.unmapped.2.fq'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.pilon.map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Map reads for pilon polishing."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --no-mixed --sensitive --no-unal -p {threads} -x {output.pref} -q -1 {input.r1} -2 {input.r2} -S {output.sam} --un-conc {output.sam}.unmapped.fq > {log} 2>&1;bowtie2 --no-mixed --sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.ru}  -S {output.sam2} >> {log} 2>&1"

rule sam_to_bam:
    input:
        sam=rules.pilon_map.output.sam,
        sam2=rules.pilon_map.output.sam2
    output:
        bam = "%s.bam"%(rules.pilon_map.output.sam),
        bam2= "%s.bam"%(rules.pilon_map.output.sam2)
    log:'%s/%s.%s.assembly.out/%s.samtools.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Convert sam to bam ."""
    shell:"samtools view -bS {input.sam} -o {output.bam} 1>> {log} 2>&1;samtools view -bS {input.sam2} -o {output.bam2} 1>> {log} 2>&1;"

rule bam_sort:
    input:
        bam = rules.sam_to_bam.output.bam,
        bam2 = rules.sam_to_bam.output.bam
    output:
        bam_sorted = "%s/%s.%s.assembly.out/sorted.bam"%(config['prefix'],config['sample'],config['iter']),
        bam_sorted2 = "%s/%s.%s.assembly.out/sorted2.bam"%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.samtools.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Sort bam ."""
    shell: "samtools sort -@ {threads} {input.bam} -o %s/%s.%s.assembly.out/sorted.bam -O bam -T tmp 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1;samtools sort -@ {threads} {input.bam2} -o %s/%s.%s.assembly.out/sorted2.bam -O bam -T tmp 1>> {log} 2>&1; samtools index {output.bam_sorted2} 1>> {log} 2>&1"%(config['prefix'],config['sample'],config['iter'],config['prefix'],config['sample'],config['iter'])

rule pilon_contigs:
    input:
        contigs=config['assembly'],
        sam = rules.bam_sort.output.bam_sorted,
        sam2 = rules.bam_sort.output.bam_sorted2
    output:
        pilonctg='%s/%s.%s.assembly.out/contigs.pilon.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.pilon.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx19G -jar %s/bin/pilon-1.22.jar --flank 5 --threads {threads} --mindepth 3 --genome {input.contigs} --frags {input.sam} --unpaired {input.sam2} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases,amb --tracks --changes 1>> {log} 2>&1"%(config["mcdir"],config['prefix'],config['sample'],config['iter'])

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

rule join_contigs:
    input:
        mc_contigs=rules.pilon_contigs.output.pilonctg,
        mh_contigs=rules.assemble_unmapped.output.megahit_contigs
    message: """---concanenate reference-guided and de novo contigs"""
    output:
        final_contigs="%s/%s.0.assembly.out/contigs.final.fasta"%(config['prefix'],config['sample'])
    shell:"cat {input.mh_contigs} {input.mc_contigs} > {output.final_contigs}"

rule create_tsv:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs=config['assembly'],
        mc_contigs_pilon=rules.pilon_contigs.output.pilonctg,
        mg_contigs=rules.assemble_unmapped.output.megahit_contigs,
        ref=config['reference']
    params:    
        minlen="%d"%(int(config['minlen']))
    message: """---information reference-guided and de novo contigs"""
    output:
        summary="%s/metacompass_summary.tsv"%(config['prefix']),
        stats="%s/metacompass_assembly_stats.tsv"%(config['prefix'])
    #change3.6 to 4 for github
    shell:"sh %s/bin/create_tsv.sh {input.mc_contigs_pilon} {input.mc_contigs} {input.mg_contigs} {input.ref} {output.summary};python3.6 %s/bin/assembly_stats.py {input.contigs} {params.minlen} > {output.stats}"%(config["mcdir"],config["mcdir"])

