configfile: "config.json"

rule build_contigs:
    input:
        genome = '{sample}.assembly.out/mc.refseq.fna',
        sam=  '{sample}.assembly.out/{sample}.sam'
    output:
        log='{sample}.assembly.out/{sample}.assembly.log',
        out='{sample}.assembly.out',
	contigs='{sample}.assembly.out/contigs.fasta'#'{sample}.contigs.fasta'
    message: """---Build contigs ."""
    shell:'./buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 2 -l 300 -n T -b T -u T -k breadth  1>> {output.log} 2>&1'#;cp {output.out}/contigs.fasta {output.contigs}'


