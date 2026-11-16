#!/usr/bin/env python3

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


STATUS_SCORE = {
    "Complete Genome": 5,
    "Chromosome": 4,
    "Scaffold": 2,
    "Contig": 1,
    "NA": 0,
}

REFSEQ_CATEGORY_SCORE = {
    "reference genome": 6,
    "representative genome": 5,
    "REFR": 6,
    "REPR": 5,
    "-": 0,
    "NA": 0,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Quality-aware reference genome filtering using NCBI metadata, "
            "CheckM/ANI summaries when available, and FASTA/sequence-report fallback."
        )
    )

    parser.add_argument("--matching-files", required=True, type=Path)
    parser.add_argument("--taxonomy-info", required=True, type=Path)
    parser.add_argument("--prokaryotes", required=True, type=Path)

    parser.add_argument(
        "--assembly-data-report",
        required=False,
        type=Path,
        default=None,
        help="Local NCBI assembly_data_report.jsonl. Contains assembly stats for many genomes.",
    )

    parser.add_argument(
        "--summary-metadata",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional richer metadata from `datasets summary genome accession`, "
            "with checkm_info and average_nucleotide_identity."
        ),
    )

    parser.add_argument(
        "--fetch-missing-summary",
        action="store_true",
        help=(
            "If enabled, use `datasets summary genome accession` once in batch "
            "to retrieve missing CheckM/ANI metadata."
        ),
    )

    parser.add_argument(
        "--summary-cache",
        type=Path,
        default=Path("datasets_summary_cache.json"),
        help="Where to write/read fetched datasets summary metadata.",
    )

    parser.add_argument(
        "--datasets-bin",
        default="datasets",
        help="Path to NCBI datasets executable.",
    )

    parser.add_argument("--output-matching-files", required=True, type=Path)
    parser.add_argument("--output-metrics", required=True, type=Path)
    parser.add_argument("--output-decisions", required=True, type=Path)

    parser.add_argument("--top-n-per-group", type=int, default=5)

    parser.add_argument("--min-normalized-n50", type=float, default=0.05)
    parser.add_argument("--max-contigs-good", type=int, default=500)

    parser.add_argument("--high-completeness", type=float, default=95.0)
    parser.add_argument("--medium-completeness", type=float, default=90.0)
    parser.add_argument("--max-high-contamination", type=float, default=5.0)
    parser.add_argument("--max-medium-contamination", type=float, default=10.0)

    parser.add_argument(
        "--prefer-refseq-representatives",
        action="store_true",
        help="Give strong priority to REFR/REPR or reference/representative genomes.",
    )

    parser.add_argument(
        "--drop-extreme-low-quality",
        action="store_true",
        help=(
            "Drop extremely fragmented or contaminated genomes even if they are the only representative. "
            "By default, the best available genome is kept for each taxonomic group."
        ),
    )

    return parser.parse_args()


def accession_from_path(path):
    path = Path(path)
    text = str(path)
    m = re.search(r"(GC[AF]_\d+\.\d+)", text)
    if m:
        return m.group(1)
    return None


def gca_from_accession(accession):
    if accession is None:
        return None
    if accession.startswith("GCA_"):
        return accession
    return None


def normalize_value(value, default="NA"):
    if value is None:
        return default
    if isinstance(value, str):
        value = value.strip()
        return value if value else default
    return value


def safe_int(value):
    try:
        if value in [None, "", "NA"]:
            return None
        return int(float(value))
    except Exception:
        return None


def safe_float(value):
    try:
        if value in [None, "", "NA"]:
            return None
        return float(value)
    except Exception:
        return None


def read_matching_files(path):
    fasta_paths = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                fasta_paths.append(line)
    return fasta_paths


def fasta_lengths(fasta_path):
    lengths = []
    current_len = 0

    try:
        with open(fasta_path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_len > 0:
                        lengths.append(current_len)
                    current_len = 0
                else:
                    current_len += len(line)

            if current_len > 0:
                lengths.append(current_len)
    except FileNotFoundError:
        pass

    return lengths


def sequence_report_lengths(fasta_path, accession):
    genome_dir = Path(fasta_path).parent
    seq_report = genome_dir / "sequence_report.jsonl"

    if not seq_report.exists():
        return []

    lengths = []
    with open(seq_report) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except Exception:
                continue

            rec_acc = record.get("assemblyAccession") or record.get("assembly_accession")
            if rec_acc and accession and rec_acc != accession:
                continue

            length = safe_int(record.get("length"))
            if length is not None:
                lengths.append(length)

    return lengths


def compute_n50(lengths):
    if not lengths:
        return 0

    total = sum(lengths)
    half = total / 2.0
    running = 0

    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length

    return 0


def compute_fasta_stats(fasta_path, accession):
    lengths = sequence_report_lengths(fasta_path, accession)

    if not lengths:
        lengths = fasta_lengths(fasta_path)

    assembly_length = sum(lengths)
    num_contigs = len(lengths)
    largest_contig = max(lengths) if lengths else 0
    n50 = compute_n50(lengths)

    return {
        "computed_assembly_length": assembly_length,
        "computed_num_contigs": num_contigs,
        "computed_largest_contig": largest_contig,
        "computed_n50": n50,
    }


def load_taxonomy(taxonomy_path):
    """
    Expected:
    accession,taxid,domain,phylum,class,order,family,genus,species
    """
    taxonomy = {}

    if not taxonomy_path or not taxonomy_path.exists():
        return taxonomy

    with open(taxonomy_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split(",")
            if len(parts) < 9:
                continue

            accession = parts[0].strip()
            taxonomy[accession] = {
                "taxid": parts[1].strip(),
                "domain": parts[2].strip(),
                "phylum": parts[3].strip(),
                "class": parts[4].strip(),
                "order": parts[5].strip(),
                "family": parts[6].strip(),
                "genus": parts[7].strip(),
                "species": parts[8].strip(),
            }

    return taxonomy


def load_prokaryotes(prokaryotes_path):
    """
    Parses NCBI-style prokaryotes.txt.

    Important fields:
      Assembly Accession
      Status
      Scaffolds
      Size (Mb)
      GC%
      Reference
    """
    metadata = {}

    if not prokaryotes_path or not prokaryotes_path.exists():
        return metadata

    with open(prokaryotes_path) as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        header = [h.strip() for h in header]
        col = {name: idx for idx, name in enumerate(header)}

        if "Assembly Accession" not in col:
            return metadata

        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            acc_idx = col["Assembly Accession"]

            if acc_idx >= len(parts):
                continue

            accession = parts[acc_idx].strip()
            if not accession:
                continue

            def get_field(name, default="NA"):
                idx = col.get(name)
                if idx is None or idx >= len(parts):
                    return default
                value = parts[idx].strip()
                return value if value else default

            metadata[accession] = {
                "organism": get_field("Organism/Name"),
                "size_mb": get_field("Size (Mb)"),
                "gc_percent": get_field("GC%"),
                "scaffolds": get_field("Scaffolds"),
                "genes": get_field("Genes"),
                "proteins": get_field("Proteins"),
                "status": get_field("Status"),
                "reference": get_field("Reference"),
                "ftp_path": get_field("FTP Path"),
                "strain": get_field("Strain"),
            }

    return metadata


def load_assembly_report(report_path):
    """
    Loads local NCBI assembly_data_report.jsonl.

    Supports older camelCase schema:
      assemblyInfo.assemblyAccession
      assemblyInfo.pairedAssembly.accession
      assemblyStats.contigN50

    Stores records by both GCA and GCF when available.
    """
    records = {}

    if not report_path or not report_path.exists():
        return records

    with open(report_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            try:
                record = json.loads(line)
            except Exception:
                continue

            assembly_info = record.get("assemblyInfo", {})
            accession = assembly_info.get("assemblyAccession")
            paired = assembly_info.get("pairedAssembly", {}).get("accession")
            refseq_acc = assembly_info.get("refseqAssmAccession")
            genbank_acc = assembly_info.get("genbankAssmAccession")
            current_acc = assembly_info.get("currentAssemblyAccession")

            keys = [accession, paired, refseq_acc, genbank_acc, current_acc]
            for key in keys:
                if key:
                    records[key] = record

    return records


def parse_summary_json_file(path):
    """
    Parses current datasets summary genome accession output.

    Supports:
      { "reports": [ ... ] }
    and JSONL records.
    """
    records = {}

    if not path or not path.exists():
        return records

    text = path.read_text().strip()
    if not text:
        return records

    parsed_any = False

    try:
        obj = json.loads(text)
        parsed_any = True
        reports = obj.get("reports", [])
        for record in reports:
            add_summary_record(records, record)
    except Exception:
        pass

    if parsed_any:
        return records

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except Exception:
                continue

            if "reports" in record:
                for rep in record["reports"]:
                    add_summary_record(records, rep)
            else:
                add_summary_record(records, record)

    return records


def add_summary_record(records, record):
    accession = record.get("accession")
    current = record.get("current_accession")
    paired = record.get("paired_accession")

    assembly_info = record.get("assembly_info", {})
    paired_assembly = assembly_info.get("paired_assembly", {})
    paired_from_info = paired_assembly.get("accession")

    keys = [accession, current, paired, paired_from_info]
    for key in keys:
        if key:
            records[key] = record


def run_datasets_summary(accessions, datasets_bin, output_path):
    """
    Runs one batch datasets summary call for all accessions.
    """
    accessions = sorted(set(a for a in accessions if a and a != "NA"))

    if not accessions:
        return False

    input_file = output_path.with_suffix(".accessions.txt")
    with open(input_file, "w") as out:
        for acc in accessions:
            out.write(acc + "\n")

    cmd = [
        datasets_bin,
        "summary",
        "genome",
        "accession",
        "--inputfile",
        str(input_file),
    ]

    sys.stderr.write("[INFO] Fetching datasets summary metadata for "
                     f"{len(accessions)} accessions\n")
    sys.stderr.write("[INFO] Command: " + " ".join(cmd) + "\n")

    try:
        with open(output_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)
    except FileNotFoundError:
        sys.stderr.write("[WARN] datasets executable not found. Skipping summary fetch.\n")
        return False
    except subprocess.CalledProcessError as exc:
        sys.stderr.write("[WARN] datasets summary failed. Skipping summary fetch.\n")
        sys.stderr.write(exc.stderr + "\n")
        return False

    return True


def get_nested(record, path, default="NA"):
    current = record
    for key in path:
        if not isinstance(current, dict):
            return default
        if key not in current:
            return default
        current = current[key]
    return normalize_value(current, default)


def extract_local_assembly_metadata(record):
    if not record:
        return {}

    assembly_info = record.get("assemblyInfo", {})
    assembly_stats = record.get("assemblyStats", {})
    annotation_info = record.get("annotationInfo", {})
    gene_counts = annotation_info.get("stats", {}).get("geneCounts", {})

    return {
        "assembly_accession_report": assembly_info.get("assemblyAccession", "NA"),
        "paired_accession_report": assembly_info.get("pairedAssembly", {}).get("accession", "NA"),
        "refseq_accession_report": assembly_info.get("refseqAssmAccession", "NA"),
        "genbank_accession_report": assembly_info.get("genbankAssmAccession", "NA"),
        "assembly_level": assembly_info.get("assemblyLevel", "NA"),
        "assembly_status": assembly_info.get("assemblyStatus", "NA"),
        "assembly_name": assembly_info.get("assemblyName", "NA"),
        "sequencing_tech": assembly_info.get("sequencingTech", "NA"),
        "submitter": assembly_info.get("submitter", "NA"),
        "organism_report": record.get("organismName", "NA"),
        "taxid_report": record.get("taxId", "NA"),
        "contig_l50": safe_int(assembly_stats.get("contigL50")),
        "contig_n50": safe_int(assembly_stats.get("contigN50")),
        "scaffold_l50": safe_int(assembly_stats.get("scaffoldL50")),
        "scaffold_n50": safe_int(assembly_stats.get("scaffoldN50")),
        "number_of_contigs": safe_int(assembly_stats.get("numberOfContigs")),
        "number_of_scaffolds": safe_int(assembly_stats.get("numberOfScaffolds")),
        "total_sequence_length": safe_int(assembly_stats.get("totalSequenceLength")),
        "total_ungapped_length": safe_int(assembly_stats.get("totalUngappedLength")),
        "gc_count": safe_int(assembly_stats.get("gcCount")),
        "gene_total": safe_int(gene_counts.get("total")),
        "gene_protein_coding": safe_int(gene_counts.get("proteinCoding")),
        "gene_pseudogene": safe_int(gene_counts.get("pseudogene")),
    }


def extract_summary_metadata(record):
    if not record:
        return {}

    assembly_info = record.get("assembly_info", {})
    assembly_stats = record.get("assembly_stats", {})
    annotation_info = record.get("annotation_info", {})
    gene_counts = annotation_info.get("stats", {}).get("gene_counts", {})
    checkm = record.get("checkm_info", {})
    ani = record.get("average_nucleotide_identity", {})

    return {
        "summary_accession": record.get("accession", "NA"),
        "summary_current_accession": record.get("current_accession", "NA"),
        "summary_paired_accession": record.get("paired_accession", "NA"),
        "summary_source_database": record.get("source_database", "NA"),
        "summary_assembly_level": assembly_info.get("assembly_level", "NA"),
        "summary_assembly_status": assembly_info.get("assembly_status", "NA"),
        "summary_assembly_name": assembly_info.get("assembly_name", "NA"),
        "summary_sequencing_tech": assembly_info.get("sequencing_tech", "NA"),
        "summary_organism": get_nested(record, ["organism", "organism_name"]),
        "summary_taxid": get_nested(record, ["organism", "tax_id"]),
        "summary_gc_percent": safe_float(assembly_stats.get("gc_percent")),
        "summary_genome_coverage": assembly_stats.get("genome_coverage", "NA"),
        "summary_number_of_contigs": safe_int(assembly_stats.get("number_of_contigs")),
        "summary_number_of_scaffolds": safe_int(assembly_stats.get("number_of_scaffolds")),
        "summary_contig_n50": safe_int(assembly_stats.get("contig_n50")),
        "summary_scaffold_n50": safe_int(assembly_stats.get("scaffold_n50")),
        "summary_contig_l50": safe_int(assembly_stats.get("contig_l50")),
        "summary_scaffold_l50": safe_int(assembly_stats.get("scaffold_l50")),
        "summary_total_sequence_length": safe_int(assembly_stats.get("total_sequence_length")),
        "summary_total_ungapped_length": safe_int(assembly_stats.get("total_ungapped_length")),
        "summary_gene_total": safe_int(gene_counts.get("total")),
        "summary_gene_protein_coding": safe_int(gene_counts.get("protein_coding")),
        "summary_gene_pseudogene": safe_int(gene_counts.get("pseudogene")),
        "checkm_completeness": safe_float(checkm.get("completeness")),
        "checkm_contamination": safe_float(checkm.get("contamination")),
        "checkm_marker_set": checkm.get("checkm_marker_set", "NA"),
        "checkm_marker_set_rank": checkm.get("checkm_marker_set_rank", "NA"),
        "ani_match_status": ani.get("match_status", "NA"),
        "ani_taxonomy_check_status": ani.get("taxonomy_check_status", "NA"),
        "ani_category": ani.get("category", "NA"),
        "best_ani": safe_float(get_nested(ani, ["best_ani_match", "ani"])),
        "best_ani_assembly": get_nested(ani, ["best_ani_match", "assembly"]),
        "best_ani_coverage": safe_float(get_nested(ani, ["best_ani_match", "assembly_coverage"])),
        "submitted_species": ani.get("submitted_species", "NA"),
    }


def pick_first(*values, default="NA"):
    for value in values:
        if value is None:
            continue
        if value == "NA":
            continue
        if value == "":
            continue
        return value
    return default


def build_query_accessions(accession, local_meta):
    """
    Prefer GCF/RefSeq accession for datasets summary when available.
    """
    candidates = []

    for key in [
        "refseq_accession_report",
        "paired_accession_report",
        "assembly_accession_report",
        "genbank_accession_report",
    ]:
        value = local_meta.get(key)
        if value and value != "NA":
            candidates.append(value)

    if accession:
        candidates.append(accession)

    return candidates


def group_key(accession, taxonomy, organism):
    tax = taxonomy.get(accession, {})

    species = tax.get("species", "NA")
    genus = tax.get("genus", "NA")

    if species and species != "NA":
        return "species:" + species

    if genus and genus != "NA":
        return "genus:" + genus

    if organism and organism != "NA":
        return "organism:" + organism

    return "accession:" + accession


def score_genome(metrics, prefer_refseq_representatives=False):
    score = 0

    ref_label = metrics["reference_label"]
    assembly_level = metrics["assembly_level"]

    if prefer_refseq_representatives:
        score += REFSEQ_CATEGORY_SCORE.get(ref_label, 0)

    score += STATUS_SCORE.get(assembly_level, 0)

    normalized_n50 = metrics["normalized_n50"]
    largest_fraction = metrics["largest_contig_fraction"]
    num_contigs = metrics["number_of_contigs"]
    num_scaffolds = metrics["number_of_scaffolds"]

    if normalized_n50 >= 0.50:
        score += 4
    elif normalized_n50 >= 0.20:
        score += 3
    elif normalized_n50 >= 0.05:
        score += 2

    if largest_fraction >= 0.50:
        score += 3
    elif largest_fraction >= 0.20:
        score += 2
    elif largest_fraction >= 0.05:
        score += 1

    if num_contigs <= 10:
        score += 3
    elif num_contigs <= 100:
        score += 2
    elif num_contigs <= 500:
        score += 1
    else:
        score -= 2

    if num_scaffolds <= 10:
        score += 2
    elif num_scaffolds <= 100:
        score += 1
    elif num_scaffolds > 500:
        score -= 1

    completeness = metrics["checkm_completeness"]
    contamination = metrics["checkm_contamination"]

    if completeness != "NA":
        if completeness >= 95:
            score += 3
        elif completeness >= 90:
            score += 2
        elif completeness >= 80:
            score += 1
        else:
            score -= 2

    if contamination != "NA":
        if contamination <= 1:
            score += 2
        elif contamination <= 5:
            score += 1
        elif contamination > 10:
            score -= 3
        elif contamination > 5:
            score -= 1

    ani_status = metrics["ani_taxonomy_check_status"]
    match_status = metrics["ani_match_status"]

    if ani_status == "OK":
        score += 1
    elif ani_status not in ["NA", "OK"]:
        score -= 2

    if match_status == "species_match":
        score += 1
    elif match_status not in ["NA", "species_match"]:
        score -= 1

    if ref_label in ["REFR", "REPR", "reference genome", "representative genome"]:
        score += 2

    return score


def quality_label(metrics, args):
    completeness = metrics["checkm_completeness"]
    contamination = metrics["checkm_contamination"]

    has_checkm = completeness != "NA" and contamination != "NA"

    assembly_level = metrics["assembly_level"]
    num_contigs = metrics["number_of_contigs"]
    num_scaffolds = metrics["number_of_scaffolds"]
    normalized_n50 = metrics["normalized_n50"]
    contig_n50 = metrics["contig_n50"]
    reference_label = metrics["reference_label"]

    is_refrep = reference_label in {
        "REFR",
        "REPR",
        "reference genome",
        "representative genome",
    }

    good_assembly_level = assembly_level in {"Complete Genome", "Chromosome"}
    acceptable_assembly_level = assembly_level in {"Complete Genome", "Chromosome", "Scaffold"}

    good_contiguity = (
        normalized_n50 >= 0.20
        and num_contigs <= 100
    )

    acceptable_contiguity = (
        normalized_n50 >= args.min_normalized_n50
        and num_contigs <= args.max_contigs_good
    )

    very_fragmented = (
        num_contigs > args.max_contigs_good
        or contig_n50 < 5000
    )

    if has_checkm:
        good_checkm = (
            completeness >= args.high_completeness
            and contamination <= args.max_high_contamination
        )

        acceptable_checkm = (
            completeness >= args.medium_completeness
            and contamination <= args.max_medium_contamination
        )

        poor_checkm = (
            completeness < args.medium_completeness
            or contamination > args.max_medium_contamination
        )

        if good_checkm and (good_contiguity or good_assembly_level or is_refrep):
            return "high"

        if acceptable_checkm and acceptable_contiguity:
            return "medium"

        if poor_checkm or very_fragmented:
            return "low"

        return "medium"

    # Fallback when CheckM is unavailable
    if (good_contiguity and acceptable_assembly_level) or (is_refrep and acceptable_contiguity):
        return "high"

    if acceptable_contiguity:
        return "medium"

    return "low"


def is_extreme_low_quality(metrics):
    if metrics["number_of_contigs"] > 2000:
        return True

    if metrics["contig_n50"] != "NA" and metrics["contig_n50"] < 5000:
        return True

    if metrics["checkm_completeness"] != "NA" and metrics["checkm_completeness"] < 70:
        return True

    if metrics["checkm_contamination"] != "NA" and metrics["checkm_contamination"] > 20:
        return True

    return False


def main():
    args = parse_args()

    taxonomy = load_taxonomy(args.taxonomy_info)
    prokaryotes = load_prokaryotes(args.prokaryotes)
    assembly_report = load_assembly_report(args.assembly_data_report)
    summary_records = parse_summary_json_file(args.summary_metadata)

    fasta_paths = read_matching_files(args.matching_files)

    preliminary_local_meta = {}
    accessions_to_fetch = set()

    for fasta_path in fasta_paths:
        accession = accession_from_path(fasta_path)
        if not accession:
            continue

        local_record = assembly_report.get(accession)
        local_meta = extract_local_assembly_metadata(local_record)
        preliminary_local_meta[accession] = local_meta

        summary_record = summary_records.get(accession)
        summary_meta = extract_summary_metadata(summary_record)

        has_checkm = (
            summary_meta.get("checkm_completeness") not in [None, "NA"]
            and summary_meta.get("checkm_contamination") not in [None, "NA"]
        )

        if not has_checkm:
            for query_acc in build_query_accessions(accession, local_meta):
                accessions_to_fetch.add(query_acc)
                break

    if args.fetch_missing_summary:
        if not args.summary_cache.exists():
            run_datasets_summary(
                accessions=accessions_to_fetch,
                datasets_bin=args.datasets_bin,
                output_path=args.summary_cache,
            )

        fetched_records = parse_summary_json_file(args.summary_cache)
        summary_records.update(fetched_records)

    metrics_by_acc = {}
    groups = defaultdict(list)

    for fasta_path in fasta_paths:
        accession = accession_from_path(fasta_path)
        if accession is None:
            sys.stderr.write(f"[WARN] Could not parse accession from {fasta_path}\n")
            continue

        fasta_stats = compute_fasta_stats(fasta_path, accession)

        local_record = assembly_report.get(accession)
        local_meta = extract_local_assembly_metadata(local_record)

        candidate_summary_keys = [accession]
        candidate_summary_keys.extend(build_query_accessions(accession, local_meta))

        summary_record = None
        for key in candidate_summary_keys:
            if key in summary_records:
                summary_record = summary_records[key]
                break

        summary_meta = extract_summary_metadata(summary_record)

        prok_meta = prokaryotes.get(accession, {})
        tax = taxonomy.get(accession, {})

        assembly_level = pick_first(
            summary_meta.get("summary_assembly_level"),
            local_meta.get("assembly_level"),
            prok_meta.get("status"),
        )

        reference_label = pick_first(
            prok_meta.get("reference"),
            summary_meta.get("summary_source_database"),
            default="NA",
        )

        total_length = pick_first(
            summary_meta.get("summary_total_sequence_length"),
            local_meta.get("total_sequence_length"),
            fasta_stats["computed_assembly_length"],
            default=0,
        )
        total_length = safe_int(total_length) or 0

        number_of_contigs = pick_first(
            summary_meta.get("summary_number_of_contigs"),
            local_meta.get("number_of_contigs"),
            fasta_stats["computed_num_contigs"],
            default=0,
        )
        number_of_contigs = safe_int(number_of_contigs) or 0

        number_of_scaffolds = pick_first(
            summary_meta.get("summary_number_of_scaffolds"),
            local_meta.get("number_of_scaffolds"),
            prok_meta.get("scaffolds"),
            number_of_contigs,
            default=0,
        )
        number_of_scaffolds = safe_int(number_of_scaffolds) or 0

        contig_n50 = pick_first(
            summary_meta.get("summary_contig_n50"),
            local_meta.get("contig_n50"),
            fasta_stats["computed_n50"],
            default=0,
        )
        contig_n50 = safe_int(contig_n50) or 0

        scaffold_n50 = pick_first(
            summary_meta.get("summary_scaffold_n50"),
            local_meta.get("scaffold_n50"),
            contig_n50,
            default=0,
        )
        scaffold_n50 = safe_int(scaffold_n50) or 0

        largest_contig = fasta_stats["computed_largest_contig"]
        normalized_n50 = contig_n50 / total_length if total_length > 0 else 0.0
        largest_contig_fraction = largest_contig / total_length if total_length > 0 else 0.0

        completeness = summary_meta.get("checkm_completeness")
        contamination = summary_meta.get("checkm_contamination")

        completeness = completeness if completeness not in [None, "NA"] else "NA"
        contamination = contamination if contamination not in [None, "NA"] else "NA"

        organism = pick_first(
            summary_meta.get("summary_organism"),
            local_meta.get("organism_report"),
            prok_meta.get("organism"),
            tax.get("species"),
        )

        metrics = {
            "accession": accession,
            "query_accession_used_for_summary": summary_meta.get("summary_accession", "NA"),
            "paired_accession": pick_first(
                summary_meta.get("summary_paired_accession"),
                local_meta.get("paired_accession_report"),
                default="NA",
            ),
            "fasta_path": fasta_path,
            "organism": organism,
            "taxid": pick_first(
                summary_meta.get("summary_taxid"),
                local_meta.get("taxid_report"),
                tax.get("taxid"),
            ),
            "domain": tax.get("domain", "NA"),
            "phylum": tax.get("phylum", "NA"),
            "class": tax.get("class", "NA"),
            "order": tax.get("order", "NA"),
            "family": tax.get("family", "NA"),
            "genus": tax.get("genus", "NA"),
            "species": tax.get("species", "NA"),
            "group": group_key(accession, taxonomy, organism),
            "reference_label": reference_label,
            "assembly_level": assembly_level,
            "assembly_status": pick_first(
                summary_meta.get("summary_assembly_status"),
                local_meta.get("assembly_status"),
                default="NA",
            ),
            "assembly_name": pick_first(
                summary_meta.get("summary_assembly_name"),
                local_meta.get("assembly_name"),
                default="NA",
            ),
            "sequencing_tech": pick_first(
                summary_meta.get("summary_sequencing_tech"),
                local_meta.get("sequencing_tech"),
                default="NA",
            ),
            "total_sequence_length": total_length,
            "number_of_contigs": number_of_contigs,
            "number_of_scaffolds": number_of_scaffolds,
            "contig_n50": contig_n50,
            "scaffold_n50": scaffold_n50,
            "largest_contig": largest_contig,
            "normalized_n50": normalized_n50,
            "largest_contig_fraction": largest_contig_fraction,
            "metadata_size_mb": prok_meta.get("size_mb", "NA"),
            "metadata_gc_percent": pick_first(
                summary_meta.get("summary_gc_percent"),
                prok_meta.get("gc_percent"),
                default="NA",
            ),
            "metadata_genome_coverage": summary_meta.get("summary_genome_coverage", "NA"),
            "metadata_scaffolds_prokaryotes": prok_meta.get("scaffolds", "NA"),
            "gene_total": pick_first(
                summary_meta.get("summary_gene_total"),
                local_meta.get("gene_total"),
                prok_meta.get("genes"),
                default="NA",
            ),
            "protein_total": pick_first(
                summary_meta.get("summary_gene_protein_coding"),
                local_meta.get("gene_protein_coding"),
                prok_meta.get("proteins"),
                default="NA",
            ),
            "checkm_completeness": completeness,
            "checkm_contamination": contamination,
            "checkm_marker_set": summary_meta.get("checkm_marker_set", "NA"),
            "checkm_marker_set_rank": summary_meta.get("checkm_marker_set_rank", "NA"),
            "ani_match_status": summary_meta.get("ani_match_status", "NA"),
            "ani_taxonomy_check_status": summary_meta.get("ani_taxonomy_check_status", "NA"),
            "ani_category": summary_meta.get("ani_category", "NA"),
            "best_ani": summary_meta.get("best_ani", "NA"),
            "best_ani_assembly": summary_meta.get("best_ani_assembly", "NA"),
            "best_ani_coverage": summary_meta.get("best_ani_coverage", "NA"),
            "submitted_species": summary_meta.get("submitted_species", "NA"),
            "metadata_source": (
                "datasets_summary"
                if summary_record
                else "assembly_data_report"
                if local_record
                else "prokaryotes_or_computed"
            ),
        }

        metrics["quality_score"] = score_genome(
            metrics,
            prefer_refseq_representatives=args.prefer_refseq_representatives,
        )
        metrics["quality_label"] = quality_label(metrics, args)
        metrics["extreme_low_quality"] = is_extreme_low_quality(metrics)

        metrics_by_acc[accession] = metrics
        groups[metrics["group"]].append(accession)

    keep = set()
    decisions = {}

    for group, accessions in groups.items():
        sorted_accessions = sorted(
            accessions,
            key=lambda acc: (
                metrics_by_acc[acc]["quality_score"],
                metrics_by_acc[acc]["checkm_completeness"]
                if metrics_by_acc[acc]["checkm_completeness"] != "NA"
                else -1,
                -metrics_by_acc[acc]["checkm_contamination"]
                if metrics_by_acc[acc]["checkm_contamination"] != "NA"
                else 0,
                metrics_by_acc[acc]["normalized_n50"],
                -metrics_by_acc[acc]["number_of_contigs"],
                metrics_by_acc[acc]["largest_contig"],
            ),
            reverse=True,
        )

        non_low = [
            acc for acc in sorted_accessions
            if metrics_by_acc[acc]["quality_label"] in {"high", "medium"}
        ]

        if non_low:
            selected = set(sorted_accessions[:args.top_n_per_group])
            for acc in accessions:
                m = metrics_by_acc[acc]

                if args.drop_extreme_low_quality and m["extreme_low_quality"]:
                    decisions[acc] = ("DROP", "extreme_low_quality")
                    continue

                if acc in selected:
                    keep.add(acc)
                    decisions[acc] = (
                        "KEEP",
                        "top_quality_representative_in_group_with_high_or_medium_quality",
                    )
                else:
                    decisions[acc] = (
                        "DROP",
                        "redundant_lower_quality_genome_in_group",
                    )
        else:
            selected = set(
                sorted_accessions[: max(1, min(args.top_n_per_group, len(sorted_accessions)))]
            )

            for acc in accessions:
                m = metrics_by_acc[acc]

                if args.drop_extreme_low_quality and m["extreme_low_quality"]:
                    decisions[acc] = ("DROP", "extreme_low_quality")
                    continue

                if acc in selected:
                    keep.add(acc)
                    decisions[acc] = (
                        "KEEP",
                        "best_available_low_quality_representative_for_group",
                    )
                else:
                    decisions[acc] = (
                        "DROP",
                        "redundant_low_quality_genome_in_group",
                    )

    with open(args.output_matching_files, "w") as out:
        for fasta_path in fasta_paths:
            acc = accession_from_path(fasta_path)
            if acc in keep:
                out.write(fasta_path + "\n")

    metric_fields = [
        "accession",
        "query_accession_used_for_summary",
        "paired_accession",
        "fasta_path",
        "organism",
        "taxid",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "group",
        "reference_label",
        "assembly_level",
        "assembly_status",
        "assembly_name",
        "sequencing_tech",
        "total_sequence_length",
        "number_of_contigs",
        "number_of_scaffolds",
        "contig_n50",
        "scaffold_n50",
        "largest_contig",
        "normalized_n50",
        "largest_contig_fraction",
        "metadata_size_mb",
        "metadata_gc_percent",
        "metadata_genome_coverage",
        "metadata_scaffolds_prokaryotes",
        "gene_total",
        "protein_total",
        "checkm_completeness",
        "checkm_contamination",
        "checkm_marker_set",
        "checkm_marker_set_rank",
        "ani_match_status",
        "ani_taxonomy_check_status",
        "ani_category",
        "best_ani",
        "best_ani_assembly",
        "best_ani_coverage",
        "submitted_species",
        "metadata_source",
        "quality_score",
        "quality_label",
        "extreme_low_quality",
    ]

    with open(args.output_metrics, "w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=metric_fields)
        writer.writeheader()
        for acc in sorted(metrics_by_acc):
            row = metrics_by_acc[acc].copy()
            row["normalized_n50"] = f"{row['normalized_n50']:.6f}"
            row["largest_contig_fraction"] = f"{row['largest_contig_fraction']:.6f}"
            writer.writerow({field: row.get(field, "NA") for field in metric_fields})

    decision_fields = [
        "accession",
        "decision",
        "reason",
        "group",
        "quality_score",
        "quality_label",
        "reference_label",
        "assembly_level",
        "checkm_completeness",
        "checkm_contamination",
        "normalized_n50",
        "number_of_contigs",
        "number_of_scaffolds",
        "contig_n50",
        "scaffold_n50",
        "ani_taxonomy_check_status",
        "metadata_source",
        "fasta_path",
    ]

    with open(args.output_decisions, "w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=decision_fields)
        writer.writeheader()

        for acc in sorted(metrics_by_acc):
            m = metrics_by_acc[acc]
            decision, reason = decisions.get(acc, ("DROP", "not_selected"))
            row = {
                "accession": acc,
                "decision": decision,
                "reason": reason,
                "group": m["group"],
                "quality_score": m["quality_score"],
                "quality_label": m["quality_label"],
                "reference_label": m["reference_label"],
                "assembly_level": m["assembly_level"],
                "checkm_completeness": m["checkm_completeness"],
                "checkm_contamination": m["checkm_contamination"],
                "normalized_n50": f"{m['normalized_n50']:.6f}",
                "number_of_contigs": m["number_of_contigs"],
                "number_of_scaffolds": m["number_of_scaffolds"],
                "contig_n50": m["contig_n50"],
                "scaffold_n50": m["scaffold_n50"],
                "ani_taxonomy_check_status": m["ani_taxonomy_check_status"],
                "metadata_source": m["metadata_source"],
                "fasta_path": m["fasta_path"],
            }
            writer.writerow(row)

    sys.stderr.write(f"[INFO] Total input references: {len(fasta_paths)}\n")
    sys.stderr.write(f"[INFO] References with metrics: {len(metrics_by_acc)}\n")
    sys.stderr.write(f"[INFO] Kept references: {len(keep)}\n")
    sys.stderr.write(f"[INFO] Dropped references: {len(metrics_by_acc) - len(keep)}\n")


if __name__ == "__main__":
    main()