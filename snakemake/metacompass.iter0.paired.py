"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: "config.json"

ruleorder: merge_reads > kmer_mask > fastq2fasta > reference_recruitment > mash_filter > bowtie2_map > build_contigs > pilon_map > sam_to_bam > bam_sort > pilon_contigs >  assemble_unmapped > join_contigs > create_tsv

if config['reads'] != "" and config['reference'] != "%s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',sample=config['sample'],prefix=config["prefix"])[0]:
     os.system("touch %s"%expand('{prefix}/{sample}.marker.match.1.fastq',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("touch %s"%expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("mkdir -p %s"%expand('{prefix}/{sample}.0.assembly.out',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("mkdir -p %s"%expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter'])[0])
     os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("touch %s"%expand('{prefix}/{sample}.0.assembly.out/mc.refseq.ids',prefix=config['prefix'],sample=config['sample'])[0])
     os.system("touch %s/%s.%s.assembly.out/%s.sam"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     os.system("touch %s/%s.%s.assembly.out/%s.sam.unmapped.fq"%(config['prefix'],config['sample'],config['iter'],config['sample']))
     os.system('touch %s/%s.%s.assembly.out'%(config['prefix'],config['sample'],config['iter']))
     os.system("touch %s/%s.%s.assembly.out/%s.index"%(config['prefix'],config['sample'],config['iter'],config['sample']))

     if int(config['iter']) >= 1:
         os.system('cp %s/%s.%d.assembly.out/contigs.final.fasta %s/%s.%s.assembly.out/contigs.fasta'%(config['prefix'],config['sample'],int(config['iter'])-1, config['prefix'],config['sample'],int(config['iter'])))

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
        len=str(int(config["length"])+3)
    threads:int(config['nthreads'])
    log:'%s/%s.%s.kmermask.log'%(config['prefix'],config['sample'],config['iter'])
    shell:"kmer-mask -ms 28 -mdb %s/refseq/kmer-mask_db/markers.mdb -1 {input.r1} -clean 0.0 -match 0.01 -nomasking -t {threads} -l {params.len} -o {params.out} 1>> {log} 2>&1"%(config["mcdir"])


rule fastq2fasta:
    input: rules.kmer_mask.output.fastq1
    output:expand('{prefix}/{sample}.fasta',prefix=config['prefix'],sample=config['sample'])
    message: """---Converting fastq to fasta."""
    shell : "perl %s/bin/fq2fa.pl -i {input} -o {output}"%(config["mcdir"])
    log:'%s/%s.%s.fastq2fasta.log'%(config['prefix'],config['sample'],config['iter'])

rule reference_recruitment:
    input:
        fasta = rules.fastq2fasta.output
        fastq = rules.merge_reads.output
    params:
        cogcov = "%d"%(int(config['cogcov'])),
        readlen = "%d"%(int(config['length']))
    output:
        out =expand('{prefix}/{sample}.{iter}.assembly.out',prefix=config['prefix'],sample=config['sample'],iter=config['iter']),
	    reffile =expand('{prefix}/{sample}.0.assembly.out/mc.refseq.fna',prefix=config['prefix'],sample=config['sample']),
        refids=expand('{prefix}/{sample}.0.assembly.out/mc.refseq.ids',prefix=config['prefix'],sample=config['sample'])
    message: """---reference recruitment."""
    threads:int(config['nthreads'])
    log:'%s/%s.%s.reference_recruitement.log'%(config['prefix'],config['sample'],config['iter'])
    shell:"mkdir -p {output.out}; python3 %s/bin/select_references.py {input.fasta} {input.fastq} {output.out} {threads} {params.cogcov}  1>> {log} 2>&1"%(config["mcdir"])

rule mash_filter:
    input:
        r1=rules.merge_reads.output.merged,
        g1=rules.reference_recruitment.output.reffile
    output:
        reffile=expand('{prefix}/{sample}.0.assembly.out/mc.refseq.filt.fna',prefix=config['prefix'],sample=config['sample'])
    params:
        mfilter= "%f"%(float(config['mfilter']))
    message: """---mash filter recruited references"""
    threads:int(config["nthreads"])
    log:'%s/%s.%s.mash.log'%(config['prefix'],config['sample'],config['iter'])
    shell:"echo {params.mfilter};python3 %s/bin/mash_filter.py {input.r1} {input.g1} {output.reffile} {params.mfilter} 1>> {log} 2>&1"%(config["mcdir"])

rule mrsfast_map:
    input:
       ref=rules.reference_recruitment.output.reffile,
       r1=rules.merge_reads.output.merged
    output:
       index= '%s.index'%(config['reference']),
       sam='%s/%s.%s.assembly.out/%s.sam'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.mrsfastmap.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Build index and map ."""
    shell:"mrsfast --index {input.ref} --ws 14 > {log} 2>&1;mrsfast --search {input.ref} --crop 100 --mem 8 --seq {input.r1} --threads {threads} -o {output.sam} --disable-nohits >> {log} 2>&1"

rule bwa_map:
    input:
       ref=rules.reference_recruitment.output.reffile,
       r1=rules.merge_reads.output.merged
    output:
       index= '%s.bwt'%(config['reference']),
       sam='%s/%s.%s.assembly.out/%s.sam'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    log: '%s/%s.%s.bwamap.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Build index and map ."""
    shell:"bwa index {input.ref} > {log} 2>&1;bwa aln -R 110 -N -o 0 -t {threads} -f {output.sam}.sai {input.ref} {input.r1} >> {log} 2>&1;bwa samse -n 1000 {input.ref} {output.sam}.sai {input.r1} > {output.sam}.full 2>&1;samtools view -F4 -@ {threads} {output.sam}.full -o {output.sam} >> {log} 2>&1"

rule bowtie2_map:
    input:
       ref=rules.reference_recruitment.output.reffile,#rules.mash_filter.output.reffile,
       r1=rules.merge_reads.output.merged
    output:
       index=expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
    log: '%s/%s.%s.bowtie2map.log'%(config['prefix'],config['sample'],config['iter'])
    threads:int(config["nthreads"])
    message: """---Build index ."""
    shell:"bowtie2-build -o 3 --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --end-to-end --sensitive --no-unal -p {threads} -x {output.pref} -q -U {input.r1} -S {output.sam}.all > {log} 2>&1; python3 %s/bin/best_strata.py {output.sam}.all {output.sam}; rm {output.sam}.all"%(config["mcdir"])

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

rule pilon_map:
    input:
       ref=rules.build_contigs.output.contigs,
       r1=config['r1'],#['reads']#.split(",")[0],
       r2=config['r2']#['reads'].split(",")[1]
    output:
       index=expand('{prefix}/{sample}.{itera}.assembly.out/{sample}.mc.index',prefix=config['prefix'],sample=config['sample'],itera=config['iter']),
       pref='%s/%s.%s.assembly.out/%s.mc.index'%(config['prefix'],config['sample'],config['iter'],config['sample']),
       sam='%s/%s.%s.assembly.out/%s.mc.sam'%(config['prefix'],config['sample'],config['iter'],config['sample']),
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
    log:'%s/%s.%s.assembly.out/%s.samtools.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:1
    message: """---Convert sam to bam ."""
    shell:"samtools view -bS {input.sam} -o {output.bam} 1>> {log} 2>&1"

rule bam_sort:
    input:
        bam = rules.sam_to_bam.output.bam 
    output:
        bam_sorted = "%s/%s.%s.assembly.out/sorted.bam"%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.samtools.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Sort bam ."""
    shell: "samtools sort -@ {threads} {input.bam} -o %s/%s.%s.assembly.out/sorted.bam -O bam -T tmp 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1"%(config['prefix'],config['sample'],config['iter'])

rule pilon_contigs:
    input:
        contigs=rules.build_contigs.output.contigs,
        sam = rules.bam_sort.output.bam_sorted
    output:
        pilonctg='%s/%s.%s.assembly.out/contigs.pilon.fasta'%(config['prefix'],config['sample'],config['iter'])
    log:'%s/%s.%s.assembly.out/%s.pilon.log'%(config['prefix'],config['sample'],config['iter'],config['sample'])
    threads:int(config['nthreads'])
    message: """---Pilon polish contigs ."""
    shell:"java -Xmx19G -jar %s/bin/pilon-1.22.jar --flank 5 --threads {threads} --mindepth 3 --genome {input.contigs} --frags {input.sam} --output %s/%s.%s.assembly.out/contigs.pilon --fix bases,amb --tracks --changes 1>> {log} 2>&1"%(config["mcdir"],config['prefix'],config['sample'],config['iter'])

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
        mc_contigs=rules.build_contigs.output.contigs,
        mc_contigs_pilon=rules.pilon_contigs.output.pilonctg,
        mg_contigs=rules.assemble_unmapped.output.megahit_contigs,
        ref=rules.reference_recruitment.output.reffile
    params:    
        minlen="%d"%(int(config['minlen']))
    message: """---information reference-guided and de novo contigs"""
    output:
        summary="%s/metacompass_summary.tsv"%(config['prefix']),
        stats="%s/metacompass_assembly_stats.tsv"%(config['prefix'])
    shell:"sh %s/bin/create_tsv.sh {input.mc_contigs_pilon} {input.mc_contigs} {input.mg_contigs} {input.ref} {output.summary};python %s/bin/assembly_stats.py {input.contigs} {params.minlen} > {output.stats}"%(config["mcdir"],config["mcdir"])
