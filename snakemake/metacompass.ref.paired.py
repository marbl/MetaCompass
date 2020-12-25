"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

import os
ruleorder: bowtie2_map > build_contigs > assembled_references >polish_map > polish_contigs >  assemble_unmapped > join_contigs > create_tsv >stats_all >stats_genome >mapping_stats

if config['reads'] != "" and config['reference'] != "%s"%expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config["outdir"][0]):
     #print("%s"%(config["outdir"]))
     #print("%s"%(config["reads"]))
     #print("%s"%(config["reference"]))
     os.system("mkdir -p %s"%expand('{outdir}/mapped_reads',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/unmapped_reads',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/assembly',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/error_correction',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/logs',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/metacompass_output',outdir=config['outdir'])[0])

rule all:
     input:
         first=expand('{outdir}/metacompass_output/metacompass_summary.tsv',outdir=config["outdir"]),
         second=expand('{outdir}/metacompass_output/metacompass_assembly_stats.tsv',outdir=config["outdir"]),
         third=expand('{outdir}/metacompass_output/metacompass_assembly_pergenome_stats.tsv',outdir=config["outdir"]),
         fourth=expand('{outdir}/metacompass_output/metacompass_mapping_stats.tsv',outdir=config["outdir"]),
         fifth=expand('{outdir}/metacompass_output/metacompass_mapping_pergenome_stats.tsv',outdir=config["outdir"])

rule bowtie2_map:
    input:
       ref=config['reference'],#rules.mash_filter.output.reffile,
       r1=config['reads'].split(",")[0]
    output:
       index=expand('{outdir}/assembly/mc.index',outdir=config['outdir']),
       pref='%s/assembly/mc.index'%(config['outdir']),
       sam='%s/assembly/mc.sam'%(config['outdir']),
       log= '%s/logs/bowtie2map.log'%(config['outdir'])
    params:
        r1=config['reads'],
        retry='%s/assembly/.run1.ok'%(config['outdir']) 
    threads:int(config["nthreads"])
    message: """---Build index ."""
    shell:"bowtie2-build -o 3 --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 -a --end-to-end --sensitive --no-unal -p {threads} -x {output.pref} -q -U {params.r1} -S {output.sam}.all > {output.log} 2>&1; python3 %s/bin/best_strata.py {output.sam}.all {output.sam}; rm {output.sam}.all && touch {params.retry}"%(config["mcdir"])

rule build_contigs:
    input:
        genome= config['reference'],
        sam=  rules.bowtie2_map.output.sam
    params:
        pickref="%s"%(config['pickref']),
        mincov="%d"%(int(config['mincov'])),
        minlen="%d"%(int(config['minlen'])),
        outputdir='%s/assembly'%(config['outdir']),
        retry='%s/assembly/.run2.ok'%(config['outdir'])   
    output:
        contigs='%s/assembly/contigs.fasta'%(config['outdir']),
        mapped_reads='%s/assembly/selected_maps.sam'%(config['outdir'])
    log:'%s/logs/buildcontigs.log'%(config['outdir'])
    threads:1
    message: """---Build contigs ."""
    shell:"%s/bin/buildcontig -r {input.genome} -s {input.sam} -o {params.outputdir} -c {params.mincov} -l {params.minlen} -n F -b F -u F -k {params.pickref}  1>> {log} 2>&1 && touch {params.retry}"%(config["mcdir"])

rule assembled_references:
    input:
        genomes=config['reference'],
        assembly=rules.build_contigs.output.contigs
    output:
        fna='%s/assembly/metacompass.assembled.fna'%(config['outdir']),
        ids='%s/assembly/metacompass.assembled.ids'%(config['outdir'])
    params:
        retry='%s/assembly/.run3.ok'%(config['outdir'])
    #log:'%s/logs/assembled_references.log'%(config['outdir'])
    threads:1
    message: """---Assembled references ."""
    shell:"grep '>' {input.assembly} |rev| cut -f2- -d '_'|rev|tr -d '>'|uniq > {output.ids};%s/bin/extractSeq {input.genomes} {output.ids} > {output.fna} && touch {params.retry}"%(config["mcdir"])

rule polish_map:
    input:
       ref=rules.build_contigs.output.contigs,
       r1=config['r1'],
       r2=config['r2']
    output:
       index=expand('{outdir}/error_correction/mc.index',outdir=config['outdir']),
       pref='%s/error_correction/mc.index'%(config['outdir']),
       sam='%s/error_correction/mc.sam'%(config['outdir']),
       unmappedr1='%s/error_correction/mc.sam.unmapped.1.fq'%(config['outdir']),
       unmappedr2='%s/error_correction/mc.sam.unmapped.2.fq'%(config['outdir']),
       log='%s/logs/polishmap.log'%(config['outdir'])
    params:
        retry='%s/error_correction/.run1.ok'%(config['outdir'])
    threads:int(config["nthreads"])
    message: """---Map reads for polish polishing."""
    shell:"bowtie2-build --threads {threads} -q {input.ref} {output.pref} 1>> {output.index} 2>&1;bowtie2 --no-mixed --sensitive --no-unal -p {threads} -x {output.pref} -q -1 {input.r1} -2 {input.r2} -S {output.sam} --un-conc {output.sam}.unmapped.fq > {output.log} 2>&1 && touch {params.retry}"
    
rule sam_to_bam:
    input:
        sam=rules.polish_map.output.sam,
    output:
        bam = "%s.bam"%(rules.polish_map.output.sam),
    params:
        retry='%s/error_correction/.run2.ok'%(config['outdir'])
    log:'%s/logs/samtools_sam2bam.log'%(config['outdir'])
    threads:1
    message: """---Convert sam to bam ."""
    shell:"samtools view -bS {input.sam} -o {output.bam} 1>> {log} 2>&1 && touch {params.retry}"

rule bam_sort:
    input:
        bam = rules.sam_to_bam.output.bam,
    output:
        bam_sorted = "%s/error_correction/sorted.bam"%(config['outdir']),
    params:
        retry='%s/error_correction/.run3.ok'%(config['outdir'])
    log:'%s/logs/samtools_bamsort.log'%(config['outdir'])
    threads:int(config['nthreads'])
    message: """---Sort bam ."""
    shell: "samtools sort -@ {threads} {input.bam} -o %s/error_correction/sorted.bam -O bam -T $RANDOM 1>> {log} 2>&1; samtools index {output.bam_sorted} 1>> {log} 2>&1 && touch {params.retry} "%(config['outdir'])

rule polish_contigs:
    input:
        contigs=rules.build_contigs.output.contigs,
        r1=config['r1'],
        r2=config['r2']
    output:
        polishctg='%s/error_correction/contigs_edited.fa'%(config['outdir'])
    params:
        memory="%d"%(int(config['memory'])),
        bf='%s/error_correction/solidBF_c1'%(config['outdir']),
        fld='%s/error_correction/contigs'%(config['outdir'])
    log:'%s/logs/polish.log'%(config['outdir'])
    threads:int(config['nthreads'])
    message: """---ntEDit polish contigs ."""
    run:
        shell("touch {output}")
        shell("/usr/bin/time -v -o {params.bf}.time nthits -c1 -b 36 -k 25 -t16 --outbloom {input.r1} {input.r2} -p {params.bf}> {log} 2>&1")
        shell("/usr/bin/time -v -o {params.fld}.time ntedit -f {input.contigs} -r {params.bf}_k25.bf -b {params.fld} -m 1 >> {log} 2>&1")
  
rule assemble_unmapped:
    input:
        r1=rules.polish_map.output.unmappedr1,
        r2=rules.polish_map.output.unmappedr2,
    output:
        megahit_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir'])
    params:
        retry='%s/assembly/.run4.ok'%(config['outdir'])
    threads:int(config["nthreads"])
    log: '%s/logs/megahit.log'%(config['outdir'])
    message: """---Assemble unmapped reads ."""
    shell:"if [[ -s {input.r1} || -s {input.r2} ]]; then rm -rf %s/assembly/megahit; megahit -o %s/assembly/megahit --min-count 3 --min-contig-len %d --presets meta-sensitive -t {threads} -1 {input.r1} -2 {input.r2}  1>> {log} 2>&1; else touch {output.megahit_contigs} {log}; echo 'No unmapped reads to run de novo assembly' >{log} ;fi && touch {params.retry}"%(config['outdir'],config['outdir'],int(config['minlen']))

rule join_contigs:
    input:
        mc_contigs=rules.polish_contigs.output.polishctg,
        mh_contigs=rules.assemble_unmapped.output.megahit_contigs
    message: """---concanenate reference-guided and de novo contigs"""
    output:
        final_contigs="%s/metacompass_output/metacompass.final.ctg.fa"%(config['outdir'])
    params:
        retry='%s/assembly/.run5.ok'%(config['outdir'])
    shell:"cat {input.mh_contigs} {input.mc_contigs} > {output.final_contigs} && touch {params.retry}"

rule create_tsv:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs=rules.build_contigs.output.contigs,
        mc_contigs_polish=rules.polish_contigs.output.polishctg,
        mg_contigs=rules.assemble_unmapped.output.megahit_contigs,
        asm_contigs=rules.assembled_references.output.fna,
        ref=config['reference']
    params:    
        minlen="%d"%(int(config['minlen'])),
        retry='%s/metacompass_output/.run1.ok'%(config['outdir'])
    message: """---information reference-guided and de novo contigs"""
    output:
        summary="%s/metacompass_output/metacompass_summary.tsv"%(config['outdir']),
    shell:"sh %s/bin/create_tsv.sh {input.mc_contigs_polish} {input.mc_contigs} {input.mg_contigs} {input.asm_contigs} {input.ref} {output.summary} && touch {params.retry}"%(config["mcdir"])

rule stats_all:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs=rules.build_contigs.output.contigs,
        mc_contigs_polish=rules.polish_contigs.output.polishctg,
        mg_contigs=rules.assemble_unmapped.output.megahit_contigs,
        ref=config['reference']
    params:    
        minlen="%d"%(int(config['minlen'])),
        out='%s/assembly'%(config['outdir']),
        retry='%s/metacompass_output/.run2.ok'%(config['outdir'])
    message: """---assembly stats reference-guided contigs"""
    output:
        stats="%s/metacompass_output/metacompass_assembly_stats.tsv"%(config['outdir']),
    shell:"python %s/bin/assembly_stats.py {input.contigs} {params.minlen} > {output.stats} && touch {params.retry}"%(config["mcdir"])
     
rule stats_genome:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs=rules.build_contigs.output.contigs,
        mc_contigs_polish=rules.polish_contigs.output.polishctg,
        mg_contigs=rules.assemble_unmapped.output.megahit_contigs
    params:
        path="%s"%(config["mcdir"]),
        out="%s"%(config['outdir']),
        assembly="assembly",
        minlen="%d"%(int(config['minlen'])),
        log='%s/logs/statspercontig.log'%(config['outdir']),
        retry='%s/metacompass_output/.run3.ok'%(config['outdir'])
    message: """---assembly stats per genome in reference-guided contigs"""
    output:
        statspercontig="%s/metacompass_output/metacompass_assembly_pergenome_stats.tsv"%(config['outdir'])
    shell:"sh %s/bin/assembly_percontig_stats.sh {params.path} {params.out} {params.assembly} {input.mc_contigs_polish} {params.minlen}  > {output.statspercontig} && touch {params.retry}"%(config["mcdir"])

rule mapping_stats:
    input:
        assembled_references=rules.assembled_references.output.ids,
        bowtie2reads= rules.bowtie2_map.output.log,# '%s/bowtie2map.log'%(config['outdir']),
        bowtie2contigs=rules.polish_map.output.log,#'%s/polishmap.log'%(config['outdir']),
        bowtie2sam=rules.build_contigs.output.mapped_reads#'%s/assembly/selected_maps.sam'%(config['outdir'])
    message: """---mapping stats per genome in reference-guided contigs"""
    log: '%s/mapping_stats.log'%(config['outdir'])
    output:
        mapping="%s/metacompass_output/metacompass_mapping_stats.tsv"%(config['outdir']),
        mappingpergenome="%s/metacompass_output/metacompass_mapping_pergenome_stats.tsv"%(config['outdir'])
    params:
        retry='%s/metacompass_output/.run4.ok'%(config['outdir']),
        references='%s/mc.refseq.ids'%(config['outdir'])
    shell:"grep '>' %s|tr -d '>' > {params.references} ;sh %s/bin/mapping_stats.sh {params.references} {input.assembled_references} {input.bowtie2reads} {input.bowtie2contigs} {input.bowtie2sam} {output.mapping} {output.mappingpergenome} && touch {params.retry}"%(config['reference'],config["mcdir"])

onsuccess:
    print("MetaCompass finished successfully!")
    os.system("touch %s/.run.ok"%(config['outdir']))

#onerror:
#    print("One or more errors occurred. See MetaCompass Log files for more info")
#    sys.exit(1)