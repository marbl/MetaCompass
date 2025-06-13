#!/usr/bin/env nextflow

// setting variables

nextflow.preview.output = 1
params.forward = "" // input forward reads
params.reverse = "" // input reverse reads
forward_gz = ""  // uncompressed forward reads
reverse_gz = ""  // uncompressed reverse reads

params.unpaired = "" // input unpaired reads
unpaired_gz = ""  // uncompressed unpaired reads

reads = "" // string containing command line parameters for passing reads along

params.reference_db = "" // location of reference database
params.output = "" // location for output files
params.skip_rs = false // skip reference selection
params.skip_rc = false // skip reference culling
params.clean_uf = false // clean up as we go
params.de_novo = 1 // execute de novo assembly

gzip_flag = false  // input files are gziped
paired = false // input files are paired
unpaired = false // input files are unpaired

// some modularization
include {filter_reads} from './pipeline/ref_selection.nf'
include {map_to_gene} from './pipeline/ref_selection.nf'
include {select_genomes} from './pipeline/ref_selection.nf'
include {collect_refs} from './pipeline/ref_culling.nf'
include {SkaniTriangle} from './pipeline/ref_culling.nf'
include {Cluster} from './pipeline/ref_culling.nf'
include {ConcatFasta} from './pipeline/ref_culling.nf'
include {IndexReads} from './pipeline/ref_culling.nf'
include {ClusterIndex} from './pipeline/ref_culling.nf'
include {interleaveReads} from './pipeline/ref_assembly.nf'
include {reduceClusters} from './pipeline/ref_assembly.nf'
include {refAssembly} from './pipeline/ref_assembly.nf'
include {deNovoAssembly} from './pipeline/denovo_assembly.nf'
include {createOutputs} from './pipeline/finalize.nf'

// Usage information
def usage(status) {
    log.info "Usage: \nnextflow run metacompass.nf \n\
            [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] | \n\
            [--forward /path/to/forwardReads --reverse /path/to/reverseReads] | \n\
            [--unpaired /path/to/unpairedReads]] \n\
            --output /path/to/outputDir"
    log.info " --reference_db   Path to reference marker gene database."
    
    log.info ""
    log.info "Required:"
    log.info ""

    log.info " --forward        Path to forward paired-end read."
    log.info " --reverse        Path to reverse paried-end read."
    log.info " --unpaired       Path to unpaired read fasta file(s)."
    log.info " --output         Path to output folder."
    
    log.info ""
    log.info "Optional:"
    log.info ""
    
    log.info " --ref_sel        Reference selection method. Default set to 'tax'."
    log.info " --ref_pick       Pick reference selection method. Default set to 'breadth'."
    log.info " --readlen        Read length for filtering. Default set to 200."
    log.info " --mincov         Minimum coverage. Default set to 1."
    log.info " --minctglen      Minimum contig length. Default set to 1."
    log.info " --run_valet      Whether or not run VALET during reference-guided assembly."
    log.info " --de_novo        Set to 0 to skip de_novo assembly. Default set to 1."
    log.info " --skip_rs        Set to true to skip ref selection process. Default set to false."
    log.info " --skip_rc        Set to true to skip ref culling process. Default set to false."
    log.info " --tracks         Tracks. Default set to false."
    log.info " --threads        Number of threads to use. Default set to 12."
    log.info " --memory         Amount of memory to use."
    log.info " --clean_uf       Remove unnecessary files while running the job. Default is false."
    log.info " --executor       Executor to use."
    log.info " --help           Print help message."

    exit status
}

if (params.help){
    usage(0)
}
// if the path is invalid, throw an error
if (params.output == "" || ! file(params.output).mkdirs() ){
    println "ERROR: Need proper output directory path!"
    usage(1)
}

// check if reference file exists
if (params.reference_db != "" && !file(params.reference_db).isDirectory() ){
    println "ERROR: Reference genome files not found!"
    usage(1)
}

// check input reads
if (params.forward != "" && params.reverse != ""){
        if(!file(params.forward).isFile() ||
       !file(params.reverse).isFile()){
        println "ERROR: Incorrect filepath to paired reads!"
        usage(1)
    }

    paired = true
}

// Dealing with unpaired reads
if (params.unpaired != ""){
    if(!file(params.unpaired).isFile()){
        println "ERROR: Incorrect filepath to unpaired reads!"
        usage(1)
    }

    unpaired = true
}

if (paired == false && unpaired == false){
    println "ERROR: Incorrect combination of reads given!"
    usage(1)
}

// Here we define the actual workflow

workflow {
  main:

// filter the reads aligned to marker genes
  if (paired == true) {
    reads = file(params.forward) + " " + file(params.reverse)
  } else {
    reads = file(params.unpaired)
  }

// keep only reads aligned to marker genes
  mapped_reads = filter_reads(reads)

  def toProcess = []

  file("${params.reference_db}/marker_index").eachFile {item ->
     if (item.isDirectory()) { 
       // println "DEBUG: got ${item.getName()}"
        toProcess = toProcess + item.getName()
     } //else {
	//println "not dir ${item.getName()}"
     //}
  }


  // mapped_reads.view() // for debugging

  // map reads to each gene and collect info
  marker_covs = map_to_gene(mapped_reads, Channel.fromList(toProcess))
  
  // generate list of genomes selected
  select_genomes(marker_covs.collect())  

  // download all the genomes needed
  collect_refs(select_genomes.out.candidates)
  // cluster the genomes based on sequence similarity
  SkaniTriangle(collect_refs.out, select_genomes.out.candidates.countLines())
  Cluster(SkaniTriangle.out)
  // create per-cluster FASTA files
  ConcatFasta(Cluster.out.clusters)
  
  // build k-mer index for reads
  if (paired == true) {
    toindex = file(params.forward) + "\n" + file(params.reverse)
  } else {
    toindex = file(params.unpaired)
  }
  IndexReads(toindex)

  // build k-mer index from clusters
  ClusterIndex(ConcatFasta.out.cluster_files.flatten()) 

  // run alignment/assembly

  // first interleave
  if (paired == true) {
    toassemble = "-in1=" + file(params.forward) + " -in2=" + file(params.reverse)
  } else {
    toassemble = "-in=" + file(params.unpaired)
  }
  interleaveReads(toassemble)


  // combine all clusters in a single fasta file
  // and identify just the reads that map to some of the clusters
  reduceClusters(Cluster.out.clusters, interleaveReads.out)

  // run the cluster-by-cluster assembly
  refAssembly(interleaveReads.out, Cluster.out.clusters, reduceClusters.out, IndexReads.out.collect(), ClusterIndex.out.collect()) 

  // de novo assembly (if needed)
  deNovoAssembly(refAssembly.out.unmapped_reads)
 
  // stage results
  if (params.de_novo > 0) {
     createOutputs(1, refAssembly.out.genomes)   
  } else {
     createOutputs(0, refAssembly.out.genomes)
  }

  publish:
  contigs = createOutputs.out

}

output {
  contigs { 
    mode 'move'
    path "."
  }
}
