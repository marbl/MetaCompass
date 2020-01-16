"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json",
ruleorder: pilon_contigs > join_contigs > create_tsv >stats_all >stats_genome >mapping_stats

rule all:
     input:
         first=expand('{outdir}/metacompass_output/metacompass_summary.tsv',outdir=config["outdir"]),
         second=expand('{outdir}/metacompass_output/metacompass_assembly_stats.tsv',outdir=config["outdir"]),
         third=expand('{outdir}/metacompass_output/metacompass_assembly_pergenome_stats.tsv',outdir=config["outdir"]),
         fourth=expand('{outdir}/metacompass_output/metacompass_mapping_stats.tsv',outdir=config["outdir"]),
         fifth=expand('{outdir}/metacompass_output/metacompass_mapping_pergenome_stats.tsv',outdir=config["outdir"]),
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

rule assemble_unmapped:
    input:
        r1='%s/error_correction/mc.sam.unmapped.1.fq'%(config['outdir']),
        r2='%s/error_correction/mc.sam.unmapped.2.fq'%(config['outdir'])
    output:
        megahit_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir'])
    params:
        retry='%s/assembly/run4.ok'%(config['outdir'])
    threads:int(config["nthreads"])
    log: '%s/logs/megahit.log'%(config['outdir'])
    message: """---Assemble unmapped reads ."""
    shell:"if [[ -s {input.r1} || -s  {input.r2} ]]; then rm -rf %s/assembly/megahit; megahit -o %s/assembly/megahit --min-count 3 --min-contig-len %d --presets meta-sensitive -t {threads} -1 {input.r1} -2 {input.r2}  1>> {log} 2>&1; else touch {output.megahit_contigs} {log}; echo 'No unmapped reads to run de novo assembly' >{log} ;fi && touch {params.retry}"%(config['outdir'],config['outdir'],int(config['minlen']))

rule join_contigs:
    input:
        mc_contigs=rules.pilon_contigs.output.pilonctg,
        mh_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir'])
    message: """---concanenate reference-guided and de novo contigs"""
    output:
        final_contigs="%s/metacompass_output/metacompass.final.ctg.fa"%(config['outdir'])
    params:
        retry='%s/assembly/run5.ok'%(config['outdir'])
    shell:"cat {input.mh_contigs} {input.mc_contigs} > {output.final_contigs} && touch {params.retry}"

rule create_tsv:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs='%s/assembly/contigs.fasta'%(config['outdir']),
        mc_contigs_pilon=rules.pilon_contigs.output.pilonctg,
        mg_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir']),
        asm_contigs='%s/assembly/metacompass.assembled.ids'%(config['outdir']),
        ref=expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config['outdir'])
    params:    
        minlen="%d"%(int(config['minlen'])),
        retry='%s/metacompass_output/run1.ok'%(config['outdir'])
    message: """---information reference-guided and de novo contigs"""
    output:
        summary="%s/metacompass_output/metacompass_summary.tsv"%(config['outdir']),
    shell:"echo bin/create_tsv.sh {input.mc_contigs_pilon} {input.mc_contigs} {input.mg_contigs} {input.asm_contigs} {input.ref} {output.summary};sh %s/bin/create_tsv.sh {input.mc_contigs_pilon} {input.mc_contigs} {input.mg_contigs} {input.asm_contigs} {input.ref} {output.summary} && touch {params.retry}"%(config["mcdir"])
        
rule stats_all:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs='%s/assembly/contigs.fasta'%(config['outdir']),
        mc_contigs_pilon=rules.pilon_contigs.output.pilonctg,
        mg_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir']),
        ref=expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config['outdir'])
    params:    
        minlen="%d"%(int(config['minlen'])),
        out='%s/assembly'%(config['outdir']),
        retry='%s/metacompass_output/run2.ok'%(config['outdir'])
    message: """---assembly stats reference-guided contigs"""
    output:
        stats="%s/metacompass_output/metacompass_assembly_stats.tsv"%(config['outdir']),
    shell:"python %s/bin/assembly_stats.py {input.contigs} {params.minlen} > {output.stats} && touch {params.retry}"%(config["mcdir"])
    
rule stats_genome:
    input:
        contigs=rules.join_contigs.output.final_contigs,
        mc_contigs='%s/assembly/contigs.fasta'%(config['outdir']),
        mc_contigs_pilon=rules.pilon_contigs.output.pilonctg,
        mg_contigs='%s/assembly/megahit/final.contigs.fa'%(config['outdir'])
    params:
        path="%s"%(config["mcdir"]),
        out="%s"%(config['outdir']),
        assembly="assembly",
        minlen="%d"%(int(config['minlen'])),
        log='%s/logs/statspercontig.log'%(config['outdir']),
        retry='%s/metacompass_output/run3.ok'%(config['outdir'])
    message: """---assembly stats per genome in reference-guided contigs"""
    output:
        statspercontig="%s/metacompass_output/metacompass_assembly_pergenome_stats.tsv"%(config['outdir'])
    shell:"echo %s/bin/assembly_percontig_stats.sh {params.path} {params.out} {params.assembly} {input.mc_contigs_pilon} {params.minlen} >{params.log};sh %s/bin/assembly_percontig_stats.sh {params.path} {params.out} {params.assembly} {input.mc_contigs_pilon} {params.minlen}  > {output.statspercontig} && touch {params.retry}"%(config["mcdir"],config["mcdir"])

rule mapping_stats:
    input:
        references=expand('{outdir}/reference_selection/mc.refseq.ids',outdir=config['outdir']),
        assembled_references='%s/assembly/metacompass.assembled.ids'%(config['outdir']),
        bowtie2reads= '%s/logs/bowtie2map.log'%(config['outdir']),# '%s/bowtie2map.log'%(config['outdir']),
        bowtie2contigs='%s/logs/pilonmap.log'%(config['outdir']),#'%s/pilonmap.log'%(config['outdir']),
        bowtie2sam='%s/assembly/selected_maps.sam'%(config['outdir'])#'%s/assembly/selected_maps.sam'%(config['outdir'])
    message: """---assembly stats per genome in reference-guided contigs"""
    log: '%s/mapping_stats.log'%(config['outdir'])
    output:
        mapping="%s/metacompass_output/metacompass_mapping_stats.tsv"%(config['outdir']),
        mappingpergenome="%s/metacompass_output/metacompass_mapping_pergenome_stats.tsv"%(config['outdir'])
    params:
        retry='%s/metacompass_output/run4.ok'%(config['outdir'])
    shell:"sh %s/bin/mapping_stats.sh {input.references} {input.assembled_references} {input.bowtie2reads} {input.bowtie2contigs} {input.bowtie2sam} {output.mapping} {output.mappingpergenome} && touch {params.retry}"%(config["mcdir"])


onsuccess:
    print("MetaCompass finished succesfully!")
    os.system("touch %s/run.ok"%(config['outdir']))

#onerror:
#    print("One or more errors occurred. See MetaCompass Log files for more info")
#    sys.exit(1)