"""
DESCRIPTION
"""
#__author__ = "Victoria Cepeda

#configfile: expand("config.json",
ruleorder: merge_reads > kmer_mask > fastq2fasta > reference_selection

if config['reads'] != "" and config['reference'] != "%s"%expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config["outdir"][0]):
     #print("%s"%(config["outdir"]))
     #print("%s"%(config["reads"]))
     #print("%s"%(config["reference"]))
     os.system("mkdir -p %s"%expand('{outdir}/intermediate_files',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/reference_selection',outdir=config['outdir'])[0])
     os.system("mkdir -p %s"%expand('{outdir}/logs',outdir=config['outdir'])[0])

rule all:
     input:
         first=expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config["outdir"])

rule merge_reads:
    input:
        reads=config['reads'].split(",")[0]
    message: """---merge fastq reads"""
    output:
        merged='%s/intermediate_files/merged.fq'%(config['outdir'])
    params:
        retry='%s/intermediate_files/run.ok'%(config['outdir'])
    run:

        for read in config['reads'].split(','):
            if read != "" and len(read) != 0:
                os.system("cat %s >> %s/intermediate_files/merged.fq && touch %s/intermediate_files/run.ok"%(read,config['outdir'],config['outdir']))


rule kmer_mask:
    input:
        r1=rules.merge_reads.output.merged
    output:
        fastq1=expand('{outdir}/reference_selection/marker.match.1.fastq',outdir=config['outdir'])[0],
    message: """---kmer-mask fastq"""
    params:
        out=expand('{outdir}/reference_selection/marker',outdir=config['outdir'])[0],
        len=str(int(config["length"])+3),
        retry='%s/reference_selection/run1.ok'%(config['outdir'])
    threads:int(config['nthreads'])
    log:'%s/logs/kmermask.log'%(config['outdir'])
    shell:"kmer-mask -ms 28 -mdb %s/refseq/kmer-mask_db/markers.mdb -1 {input.r1} -clean 0.0 -match 0.01 -nomasking -t {threads} -l {params.len} -o {params.out} 1>> {log} 2>&1  && touch {params.retry}"%(config["mcdir"])


rule fastq2fasta:
    input: rules.kmer_mask.output.fastq1
    output:expand('{outdir}/reference_selection/masked_reads.fasta',outdir=config['outdir'],sample=config['sample'])
    message: """---Converting fastq to fasta."""
    params:
        retry='%s/reference_selection/run2.ok'%(config['outdir'])
    log:'%s/logs/fastq2fasta.log'%(config['outdir'])
    shell : "%s/bin/fq2fa -i {input} -o {output} && touch {params.retry}"%(config["mcdir"])


rule reference_selection:
    input:
        fasta = rules.fastq2fasta.output,
        fastq = rules.merge_reads.output
    params:
        cogcov = "%d"%(int(config['cogcov'])),
        identity = "%s"%(config['ani']),
        readlen = "%d"%(int(config['length'])),
        refsel = "%s"%(config['refsel']),
        out =expand('{outdir}/reference_selection',outdir=config['outdir']),
        retry='%s/reference_selection/run3.ok'%(config['outdir'])
    output:
        reffile =expand('{outdir}/reference_selection/mc.refseq.fna',outdir=config['outdir']),
        refids=expand('{outdir}/reference_selection/mc.refseq.ids',outdir=config['outdir'])
    message: """---reference selection."""
    threads:int(config['nthreads'])
    log:'%s/logs/reference_selection.log'%(config['outdir'])
    shell:"python3 %s/bin/select_references.py {params.refsel} {input.fasta} {input.fastq} {params.out} {threads} {params.cogcov} {params.identity} 1>> {log} 2>&1 && touch {params.retry}"%(config["mcdir"])
    
onsuccess:
    print("MetaCompass finished successfully!")
    os.system("touch %s/run.ok"%(config['outdir']))
