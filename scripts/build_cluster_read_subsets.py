#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Set, TextIO


CLUSTER_RE = re.compile(r"cluster_(\d+)")
PAIR_SUFFIX_RE = re.compile(r"/[12]$")


def open_text(path: Path, mode: str = "rt") -> TextIO:
    """
    Open plain text or gzipped text transparently.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_read_id(raw: str) -> str:
    """
    Normalize read IDs so the same logical read can be matched across:
      - minimap2 output
      - FASTQ headers

    Rules:
      - remove leading '@' if present
      - keep only the first whitespace-delimited token
      - strip trailing /1 or /2
    """
    raw = raw.strip()
    if raw.startswith("@"):
        raw = raw[1:]
    raw = raw.split()[0]
    raw = PAIR_SUFFIX_RE.sub("", raw)
    return raw


def parse_cluster_id(value: str) -> int:
    """
    Parse cluster ID from strings like:
      cluster_7
      cluster_7|GCA_000247755.2|AFVV01000034.1
    """
    m = CLUSTER_RE.search(value)
    if not m:
        raise ValueError(f"Could not parse cluster id from value: {value}")
    return int(m.group(1))


def count_clusters(cluster_list_path: Path) -> int:
    """
    Count the number of non-empty rows in the original clusters.csv.
    This determines how many cluster_<id>.fq files we should create.
    """
    n = 0
    with open(cluster_list_path, "r", newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if row:
                n += 1
    return n


def load_read_to_clusters(ref_map_path: Path, num_clusters: int) -> DefaultDict[str, Set[int]]:
    """
    Load read -> cluster assignments.

    Expected main format:
        read_id<TAB>cluster_<id>|GCA_xxx|contig_id

    Also supports:
        read_id<TAB>ref_name<TAB>cluster_<id>
        read_id<TAB>ref_name<TAB>7

    Returns:
        dict: normalized_read_id -> set(cluster_ids)
    """
    read_to_clusters: DefaultDict[str, Set[int]] = defaultdict(set)
    bad_lines = 0
    total_lines = 0

    with open_text(ref_map_path, "rt") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue

            total_lines += 1
            parts = line.split("\t")

            if len(parts) < 2:
                bad_lines += 1
                print(
                    f"[WARN] Skipping malformed line {line_no} in {ref_map_path}: {line}",
                    file=sys.stderr,
                )
                continue

            read_id = normalize_read_id(parts[0])

            cluster_id = None

            # Preferred: parse cluster from 2nd column ref header
            try:
                cluster_id = parse_cluster_id(parts[1])
            except ValueError:
                pass

            # Fallback: parse from 3rd column if present
            if cluster_id is None and len(parts) >= 3:
                third = parts[2]
                try:
                    if third.isdigit():
                        cluster_id = int(third)
                    else:
                        cluster_id = parse_cluster_id(third)
                except ValueError:
                    cluster_id = None

            if cluster_id is None:
                bad_lines += 1
                print(
                    f"[WARN] Could not determine cluster for line {line_no} in {ref_map_path}: {line}",
                    file=sys.stderr,
                )
                continue

            if cluster_id < 1 or cluster_id > num_clusters:
                bad_lines += 1
                print(
                    f"[WARN] Cluster id {cluster_id} out of range 1..{num_clusters} "
                    f"at line {line_no} in {ref_map_path}",
                    file=sys.stderr,
                )
                continue

            read_to_clusters[read_id].add(cluster_id)

    print(
        f"[INFO] Loaded {len(read_to_clusters)} unique reads from {ref_map_path} "
        f"across {total_lines} mapping rows ({bad_lines} skipped).",
        file=sys.stderr,
    )
    return read_to_clusters


def initialize_output_files(outdir: Path, num_clusters: int) -> None:
    """
    Create empty cluster_<id>.fq files up front so downstream join logic
    can rely on their presence even when a cluster gets zero reads.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    for cluster_id in range(1, num_clusters + 1):
        (outdir / f"cluster_{cluster_id}.fq").touch()


def write_cluster_fastqs(
    fq_path: Path,
    read_to_clusters: Dict[str, Set[int]],
    outdir: Path,
) -> None:
    """
    Single-pass FASTQ splitter.

    For each FASTQ record:
      - normalize read id
      - look up which cluster(s) it belongs to
      - write the full record to each corresponding cluster_<id>.fq

    Because cluster membership is stored as a set, a read will not be written
    twice to the same cluster even if it mapped to multiple refs in that cluster.
    """
    handles: Dict[int, TextIO] = {}
    total_records = 0
    matched_records = 0
    cluster_record_counts: DefaultDict[int, int] = defaultdict(int)

    try:
        with open_text(fq_path, "rt") as fh:
            while True:
                header = fh.readline()
                if not header:
                    break

                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()

                if not seq or not plus or not qual:
                    raise ValueError(f"Malformed FASTQ: incomplete record near record {total_records + 1} in {fq_path}")

                total_records += 1
                read_id = normalize_read_id(header)
                clusters = read_to_clusters.get(read_id)

                if not clusters:
                    continue

                record = f"{header}{seq}{plus}{qual}"

                for cluster_id in clusters:
                    if cluster_id not in handles:
                        handles[cluster_id] = open(outdir / f"cluster_{cluster_id}.fq", "a")
                    handles[cluster_id].write(record)
                    cluster_record_counts[cluster_id] += 1
                    matched_records += 1

    finally:
        for handle in handles.values():
            handle.close()

    print(f"[INFO] FASTQ records scanned: {total_records}", file=sys.stderr)
    print(f"[INFO] FASTQ record writes performed: {matched_records}", file=sys.stderr)

    nonempty = sum(1 for _, count in cluster_record_counts.items() if count > 0)
    print(f"[INFO] Non-empty cluster FASTQs: {nonempty}", file=sys.stderr)

    for cluster_id in sorted(cluster_record_counts):
        print(
            f"[INFO] cluster_{cluster_id}.fq records: {cluster_record_counts[cluster_id]}",
            file=sys.stderr,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build one FASTQ subset per cluster using read-to-reference mappings "
            "from reduceClusters."
        )
    )
    parser.add_argument(
        "--cluster-list",
        required=True,
        type=Path,
        help="Original clusters.csv file with one row per cluster",
    )
    parser.add_argument(
        "--bam",
        required=False,
        type=Path,
        default=None,
        help="Mapped BAM from global alignment (currently unused, kept for future use)",
    )
    parser.add_argument(
        "--fq",
        required=True,
        type=Path,
        help="Global mapped FASTQ produced by reduceClusters (e.g. concat_refs_mapped.fq)",
    )
    parser.add_argument(
        "--ref-map",
        required=True,
        type=Path,
        help=(
            "Read-to-reference mapping file. Expected primary format: "
            "read_id<TAB>cluster_<id>|GCA_xxx|contig_id"
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for cluster FASTQs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.cluster_list.exists():
        raise FileNotFoundError(f"Cluster list not found: {args.cluster_list}")
    if not args.fq.exists():
        raise FileNotFoundError(f"Mapped FASTQ not found: {args.fq}")
    if not args.ref_map.exists():
        raise FileNotFoundError(f"Reference map not found: {args.ref_map}")

    num_clusters = count_clusters(args.cluster_list)
    if num_clusters == 0:
        raise ValueError(f"No clusters found in: {args.cluster_list}")

    print(f"[INFO] Number of clusters: {num_clusters}", file=sys.stderr)

    initialize_output_files(args.outdir, num_clusters)

    read_to_clusters = load_read_to_clusters(args.ref_map, num_clusters)
    write_cluster_fastqs(args.fq, read_to_clusters, args.outdir)

    print("[INFO] build_cluster_read_subsets.py completed successfully", file=sys.stderr)


if __name__ == "__main__":
    main()