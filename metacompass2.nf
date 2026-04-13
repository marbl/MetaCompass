#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * MetaCompass main pipeline
 *
 * Responsibilities:
 *   1. Validate user inputs
 *   2. Build read arguments for downstream modules
 *   3. Orchestrate reference selection, culling, assembly, and final output staging
 *
 * Best-practice note:
 *   Parameter defaults are also in `nextflow.config`.
 */

// -----------------------------------------------------------------------------
// Parameter defaults
// -----------------------------------------------------------------------------

params.help         = params.help ?: false
params.forward      = params.forward ?: ''
params.reverse      = params.reverse ?: ''
params.unpaired     = params.unpaired ?: ''
params.output = params.output ?: "${launchDir}/results"
workflow.outputDir = file(params.output)
params.reference_db = params.reference_db ?: ''
params.ref_cache    = params.ref_cache ?: params.output

if (params.skip_rs  == null) params.skip_rs  = false
if (params.skip_rc  == null) params.skip_rc  = false
if (params.clean_uf == null) params.clean_uf = false
if (params.denovo  == null) params.denovo  = 1

// -----------------------------------------------------------------------------
// Module imports
// -----------------------------------------------------------------------------

include { filter_reads  } from './pipeline/ref_selection.nf'
include { map_to_gene   } from './pipeline/ref_selection.nf'
include { select_genomes } from './pipeline/ref_selection.nf'

include { collect_refs  } from './pipeline/ref_culling.nf'
include { SkaniTriangle } from './pipeline/ref_culling.nf'
include { Cluster       } from './pipeline/ref_culling.nf'
include { ConcatFasta   } from './pipeline/ref_culling.nf'
include { IndexReads    } from './pipeline/ref_culling.nf'
include { ClusterIndex  } from './pipeline/ref_culling.nf'

include { interleaveReads } from './pipeline/ref_assembly.nf'
include { reduceClusters  } from './pipeline/ref_assembly.nf'
include { refAssembly     } from './pipeline/ref_assembly.nf'

include { deNovoAssembly } from './pipeline/denovo_assembly.nf'
include { createOutputs  } from './pipeline/finalize.nf'

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

/**
 * Returns CLI usage text.
 */
String usageText() {
    return """
    Usage:
      nextflow run metacompass.nf \\
        [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] |
         [--forward /path/to/forwardReads --reverse /path/to/reverseReads] |
         [--unpaired /path/to/unpairedReads]] \\
        --reference_db /path/to/reference_db \\
        --output /path/to/outputDir

    Required:
      --reference_db   Path to reference marker gene database
      --output         Path to output folder
      --forward        Path to forward paired-end read
      --reverse        Path to reverse paired-end read
      --unpaired       Path to unpaired read file(s)

    Optional:
      --ref_cache      Path to cache reference genomes (default: output directory)
      --ref_sel        Reference selection method (default: tax)
      --ref_pick       Reference picking method (default: breadth)
      --readlen        Read length for filtering (default: 200)
      --mincov         Minimum coverage (default: 1)
      --minctglen      Minimum contig length (default: 1)
      --run_valet      Run VALET during reference-guided assembly
      --denovo         Set to 0 to skip de novo assembly (default: 1)
      --skip_rs        Skip reference selection (default: false)
      --skip_rc        Skip reference culling (default: false)
      --tracks         Tracks output (default: false)
      --threads        Number of threads to use (default: 12)
      --memory         Amount of memory to use
      --clean_uf       Remove unnecessary files while running (default: false)
      --executor       Executor to use
      --help           Print this help message
    """.stripIndent().trim()
}

/**
 * Print usage text and exit.
 */
void printUsageAndExit(int status = 0) {
    log.info usageText()
    System.exit(status)
}

/**
 * Log an error, show usage, and exit.
 */
void failWithUsage(String message) {
    log.error message
    printUsageAndExit(1)
}

/**
 * Ensure a file exists and is a regular file.
 */
void requireFile(String path, String flagName) {
    def f = file(path)
    if (!f.exists() || !f.isFile()) {
        failWithUsage("Invalid path for ${flagName}: ${path}")
    }
}

/**
 * Ensure a directory exists.
 */
void requireDirectory(String path, String flagName) {
    def d = file(path)
    if (!d.exists() || !d.isDirectory()) {
        failWithUsage("Invalid directory for ${flagName}: ${path}")
    }
}

/**
 * Validate CLI inputs and infer read mode.
 *
 * Returns:
 *   [paired: boolean, unpaired: boolean]
 */
Map validateInputs() {
    if (params.help) {
        printUsageAndExit(0)
    }

    final boolean hasForward  = params.forward as boolean
    final boolean hasReverse  = params.reverse as boolean
    final boolean hasUnpaired = params.unpaired as boolean
    final boolean hasPaired   = hasForward && hasReverse
    final boolean hasPartialPair = (hasForward && !hasReverse) || (!hasForward && hasReverse)

    if (!(params.reference_db as boolean)) {
        failWithUsage("Missing required parameter: --reference_db")
    }

    requireDirectory(params.reference_db, '--reference_db')

    if (hasPartialPair) {
        failWithUsage("Paired-end input requires both --forward and --reverse.")
    }

    if (!hasPaired && !hasUnpaired) {
        failWithUsage("No reads supplied. Provide paired-end reads or unpaired reads.")
    }

    if (hasPaired) {
        requireFile(params.forward, '--forward')
        requireFile(params.reverse, '--reverse')
    }

    if (hasUnpaired) {
        requireFile(params.unpaired, '--unpaired')
    }

    return [
        paired  : hasPaired,
        unpaired: hasUnpaired
    ]
}

/**
 * Discover marker-index subdirectories under the reference database.
 */
List<String> discoverMarkerTargets(String referenceDb) {
    def markerIndexDir = file("${referenceDb}/marker_index")

    if (!markerIndexDir.exists() || !markerIndexDir.isDirectory()) {
        failWithUsage("Missing marker index directory: ${markerIndexDir}")
    }

    def targets = []
    markerIndexDir.eachFile { item ->
        if (item.isDirectory()) {
            targets << item.name
        }
    }

    if (targets.isEmpty()) {
        failWithUsage("No marker-index subdirectories found in: ${markerIndexDir}")
    }

    return targets.sort()
}

/**
 * Build read arguments for downstream modules based on input mode.
 *
 * Note:
 *   This preserves the original behavior:
 *   - If paired reads are provided, downstream steps use paired input
 *   - Otherwise, downstream steps use unpaired input
 */
Map buildReadArguments(boolean pairedInput) {
    if (pairedInput) {
        return [
            filterArg    : "${file(params.forward)} ${file(params.reverse)}",
            indexArg     : "${file(params.forward)}\n${file(params.reverse)}",
            interleaveArg: "-in1=${file(params.forward)} -in2=${file(params.reverse)}"
        ]
    }

    return [
        filterArg    : "${file(params.unpaired)}",
        indexArg     : "${file(params.unpaired)}",
        interleaveArg: "-in=${file(params.unpaired)}"
    ]
}

// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------

workflow {
    final Map inputMode      = validateInputs()
    final Map readArgs       = buildReadArguments(inputMode.paired)
    final List<String> markerTargets = discoverMarkerTargets(params.reference_db)

    log.info "Output directory: ${params.output}"
    log.info "Input mode: ${inputMode.paired ? 'paired-end' : 'unpaired'}"
    log.info "Discovered ${markerTargets.size()} marker-index partitions"

    /*
     * Step 1: Filter reads aligned to marker genes
     */
    def mappedReads = filter_reads(readArgs.filterArg)

    /*
     * Step 2: Map reads to each marker gene partition
     */
    def markerCoverages = map_to_gene(mappedReads, Channel.fromList(markerTargets))

    /*
     * Step 3: Select candidate genomes
     */
    select_genomes(markerCoverages.collect())

    /*
     * Step 4: Collect and cluster reference genomes
     */
    collect_refs(select_genomes.out.candidates)
    SkaniTriangle(collect_refs.out, select_genomes.out.candidates.countLines())
    Cluster(SkaniTriangle.out)
    ConcatFasta(Cluster.out.clusters)

    /*
     * Step 5: Build read and cluster indexes
     */
    IndexReads(readArgs.indexArg)
    ClusterIndex(ConcatFasta.out.cluster_files.flatten())

    /*
     * Step 6: Interleave reads and reduce candidate clusters
     */
    interleaveReads(readArgs.interleaveArg)
    reduceClusters(Cluster.out.clusters, interleaveReads.out)

    /*
     * Step 7: Reference-guided assembly
     */
    refAssembly(
        interleaveReads.out,
        Cluster.out.clusters,
        reduceClusters.out,
        IndexReads.out.collect(),
        ClusterIndex.out.collect()
    )

    /*
     * Step 8: Optional de novo assembly and final output staging
     */
    if ((params.denovo as Integer) > 0) {
        deNovoAssembly(refAssembly.out.unmapped_reads)
        createOutputs(deNovoAssembly.out.flag, refAssembly.out.genomes)
    } else {
        createOutputs(0, refAssembly.out.genomes)
    }
}
