configfile: "config.json"

rule build_contigs:
    input:
        genome = '{prefix}/{sample}.assembly.out/{reference}',#mc.refseq.fna',
        sam=  '{prefix}/{sample}.assembly.out/{sample}.{iter}.sam'
    output:
        out='{prefix}/{sample}.{iter}.assembly.out',
	#contigs='{prefix}/{sample}.assembly.out/contigs.fasta'#'{sample}.contigs.fasta'
        contigs='{prefix}/{sample}.{iter}.assembly.out/contigs.fasta'
    log:'{prefix}/{sample}.assembly.out/{sample}.assembly.log'
    threads:1
    message: """---Build contigs ."""
    shell:'../bin/buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 2 -l 300 -n T -b T -u T -k breadth  1>> {log} 2>&1'


