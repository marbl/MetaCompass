configfile: "config2.json"

rule build_contigs:
    input:
        genome = config["reference"],
        sam=  '{prefix}/{sample}.assembly.out/{sample}.sam'
    output:
        out='{prefix}/{sample}.assembly.out',
	contigs='{prefix}/{sample}.assembly.out/contigs.fasta'
    log:'{prefix}/{sample}.assembly.out/{sample}.assembly.log'
    threads:config["nthreads"]
    message: """---Build contigs ."""
    shell:'../bin/buildcontig -r {input.genome} -s {input.sam} -o {output.out} -c 2 -l 300 -n T -b T -u T -k breadth  1>> {log} 2>&1'


