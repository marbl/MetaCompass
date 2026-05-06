#!/usr/bin/env python3

import os
import re
import argparse
from urllib.parse import unquote, quote


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate MetaCompass reference-guided contigs using RefSeq GFF files."
    )

    parser.add_argument(
        "--gff-base",
        required=True,
        help=(
            "Base directory containing per-genome RefSeq GFF files. "
            "Example: /path/to/RefSeq_V2_db/collect_refs/ncbi_dataset/data"
        ),
    )

    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Output TSV file with detailed contig-gene annotation.",
    )

    parser.add_argument(
        "--output-gff",
        required=True,
        help="Output GFF3 file with contig-relative feature coordinates.",
    )

    parser.add_argument(
        "--refctgs",
        required=True,
        nargs="+",
        help="Per-genome MetaCompass *.refctgs.fna files.",
    )

    return parser.parse_args()


def parse_attrs(attr_string):
    attrs = {}

    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = unquote(value)

    return attrs


def gff_escape(value):
    if value is None:
        value = "NA"

    value = str(value)

    if value == "":
        value = "NA"

    return quote(value, safe=":/._|-")


def normalize_refseq_id(raw_ref):
    """
    MetaCompass contig headers usually encode RefSeq accessions like:
      CP021391_1

    RefSeq GFF files usually use:
      CP021391.1

    This converts only the final underscore-number suffix to dot-number.
    """
    m = re.match(r"^(.+)_([0-9]+)$", raw_ref)

    if m:
        return m.group(1) + "." + m.group(2)

    return raw_ref


def parse_contig_header(header):
    """
    Example MetaCompass header:
      >CP021391_1_382803_418215_pilon

    Interpreted as:
      reference seqid: CP021391.1
      reference start: 382803
      reference end:   418215
    """
    contig_id = header.lstrip(">").strip().split()[0]

    m = re.match(r"^(.+)_([0-9]+)_([0-9]+)_pilon$", contig_id)

    if not m:
        return None

    raw_ref = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))

    if start > end:
        start, end = end, start

    return {
        "contig_id": contig_id,
        "ref_seqid": normalize_refseq_id(raw_ref),
        "contig_ref_start": start,
        "contig_ref_end": end,
    }


def read_fasta_records(fasta_path):
    records = []
    header = None
    seq_parts = []

    with open(fasta_path) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))

                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        if header is not None:
            records.append((header, "".join(seq_parts)))

    return records


def load_gff_features(gff_path):
    """
    Load functional gene features from a RefSeq GFF file.

    We intentionally skip generic 'gene' rows because they usually duplicate
    CDS/RNA rows and often do not contain product or protein_id information.
    """
    keep_types = {"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"}
    features_by_seqid = {}

    with open(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts

            if feature_type not in keep_types:
                continue

            attr_dict = parse_attrs(attrs)

            feature = {
                "seqid": seqid,
                "source": source,
                "feature_type": feature_type,
                "gene_start": int(start),
                "gene_end": int(end),
                "strand": strand,
                "phase": phase,
                "gene": attr_dict.get("gene", "NA"),
                "locus_tag": attr_dict.get("locus_tag", "NA"),
                "protein_id": attr_dict.get("protein_id", "NA"),
                "product": attr_dict.get("product", "NA"),
                "feature_id": attr_dict.get("ID", "NA"),
                "parent": attr_dict.get("Parent", "NA"),
            }

            features_by_seqid.setdefault(seqid, []).append(feature)

    return features_by_seqid


def overlap_interval(a_start, a_end, b_start, b_end):
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)

    if overlap_start > overlap_end:
        return 0, None, None

    overlap_bp = overlap_end - overlap_start + 1
    return overlap_bp, overlap_start, overlap_end


def genome_from_refctg_filename(path):
    """
    Example:
      GCA_002838365.1_ASM283836v1_genomic.refctgs.fna

    Returns:
      GCA_002838365.1
    """
    name = os.path.basename(path)
    m = re.match(r"^(GCA_\d+\.\d+)", name)

    if not m:
        return None

    return m.group(1)


def make_unique_gff_id(genome, contig_id, feat, ref_overlap_start, ref_overlap_end):
    base_id = feat.get("feature_id", "NA")

    if base_id == "NA":
        base_id = (
            feat["feature_type"]
            + "-"
            + feat.get("locus_tag", "unknown")
            + "-"
            + str(ref_overlap_start)
            + "-"
            + str(ref_overlap_end)
        )

    raw_id = (
        genome
        + "|"
        + contig_id
        + "|"
        + base_id
        + "|"
        + str(ref_overlap_start)
        + "-"
        + str(ref_overlap_end)
    )

    return gff_escape(raw_id)


def make_gff_attributes(genome, contig_id, ref_seqid, feat, ref_overlap_start, ref_overlap_end):
    attrs = []

    unique_id = make_unique_gff_id(
        genome,
        contig_id,
        feat,
        ref_overlap_start,
        ref_overlap_end,
    )

    attrs.append("ID=" + unique_id)

    if feat["parent"] != "NA":
        attrs.append("Parent=" + gff_escape(feat["parent"]))

    if feat["locus_tag"] != "NA":
        attrs.append("locus_tag=" + gff_escape(feat["locus_tag"]))

    if feat["gene"] != "NA":
        attrs.append("gene=" + gff_escape(feat["gene"]))

    if feat["protein_id"] != "NA":
        attrs.append("protein_id=" + gff_escape(feat["protein_id"]))

    if feat["product"] != "NA":
        attrs.append("product=" + gff_escape(feat["product"]))

    attrs.append("ref_genome=" + gff_escape(genome))
    attrs.append("ref_seqid=" + gff_escape(ref_seqid))
    attrs.append("ref_start=" + str(ref_overlap_start))
    attrs.append("ref_end=" + str(ref_overlap_end))

    return ";".join(attrs)


def write_outputs(gff_base, out_tsv, out_gff, refctg_files):
    gff_cache = {}

    with open("missing_gff_files.txt", "w") as missing_gff, \
         open("unparsed_contig_headers.txt", "w") as unparsed, \
         open(out_tsv, "w") as tsv, \
         open(out_gff, "w") as gffout:

        tsv.write(
            "genome_accession\t"
            "refctg_file\t"
            "contig_id\t"
            "ref_seqid\t"
            "contig_ref_start\t"
            "contig_ref_end\t"
            "contig_length\t"
            "feature_type\t"
            "gene_start_on_reference\t"
            "gene_end_on_reference\t"
            "overlap_bp\t"
            "gene_fraction_covered\t"
            "contig_feature_start\t"
            "contig_feature_end\t"
            "feature_strand\t"
            "gene\t"
            "locus_tag\t"
            "protein_id\t"
            "product\t"
            "feature_id\t"
            "parent\n"
        )

        gffout.write("##gff-version 3\n")

        for refctg in sorted(refctg_files):
            genome = genome_from_refctg_filename(refctg)

            if genome is None:
                unparsed.write("Could not parse genome accession from file: " + refctg + "\n")
                continue

            gff_path = os.path.join(gff_base, genome, "genomic.gff")

            if not os.path.exists(gff_path):
                missing_gff.write(genome + "\t" + gff_path + "\n")
                continue

            if genome not in gff_cache:
                gff_cache[genome] = load_gff_features(gff_path)

            features_by_seqid = gff_cache[genome]
            records = read_fasta_records(refctg)

            for header, sequence in records:
                contig = parse_contig_header(header)

                if contig is None:
                    unparsed.write(
                        genome
                        + "\t"
                        + os.path.basename(refctg)
                        + "\t"
                        + header.strip()
                        + "\n"
                    )
                    continue

                contig_id = contig["contig_id"]
                ref_seqid = contig["ref_seqid"]
                contig_ref_start = contig["contig_ref_start"]
                contig_ref_end = contig["contig_ref_end"]
                contig_len = len(sequence)

                gffout.write(
                    "##sequence-region "
                    + contig_id
                    + " 1 "
                    + str(contig_len)
                    + "\n"
                )

                features = features_by_seqid.get(ref_seqid, [])

                for feat in features:
                    overlap_bp, ref_overlap_start, ref_overlap_end = overlap_interval(
                        contig_ref_start,
                        contig_ref_end,
                        feat["gene_start"],
                        feat["gene_end"],
                    )

                    if overlap_bp == 0:
                        continue

                    gene_len = feat["gene_end"] - feat["gene_start"] + 1
                    gene_fraction_covered = overlap_bp / gene_len if gene_len > 0 else 0.0

                    contig_feature_start = ref_overlap_start - contig_ref_start + 1
                    contig_feature_end = ref_overlap_end - contig_ref_start + 1

                    if contig_feature_start < 1:
                        contig_feature_start = 1

                    if contig_feature_end > contig_len:
                        contig_feature_end = contig_len

                    tsv.write(
                        genome + "\t"
                        + os.path.basename(refctg) + "\t"
                        + contig_id + "\t"
                        + ref_seqid + "\t"
                        + str(contig_ref_start) + "\t"
                        + str(contig_ref_end) + "\t"
                        + str(contig_len) + "\t"
                        + feat["feature_type"] + "\t"
                        + str(feat["gene_start"]) + "\t"
                        + str(feat["gene_end"]) + "\t"
                        + str(overlap_bp) + "\t"
                        + format(gene_fraction_covered, ".6f") + "\t"
                        + str(contig_feature_start) + "\t"
                        + str(contig_feature_end) + "\t"
                        + feat["strand"] + "\t"
                        + feat["gene"] + "\t"
                        + feat["locus_tag"] + "\t"
                        + feat["protein_id"] + "\t"
                        + feat["product"] + "\t"
                        + feat["feature_id"] + "\t"
                        + feat["parent"] + "\n"
                    )

                    gff_attrs = make_gff_attributes(
                        genome,
                        contig_id,
                        ref_seqid,
                        feat,
                        ref_overlap_start,
                        ref_overlap_end,
                    )

                    gffout.write(
                        contig_id + "\t"
                        + "MetaCompass" + "\t"
                        + feat["feature_type"] + "\t"
                        + str(contig_feature_start) + "\t"
                        + str(contig_feature_end) + "\t"
                        + "." + "\t"
                        + feat["strand"] + "\t"
                        + feat["phase"] + "\t"
                        + gff_attrs + "\n"
                    )


def main():
    args = parse_args()

    write_outputs(
        gff_base=args.gff_base,
        out_tsv=args.output_tsv,
        out_gff=args.output_gff,
        refctg_files=args.refctgs,
    )


if __name__ == "__main__":
    main()