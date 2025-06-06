def test_align_reads(read_aligner_integration):
    read_aligner_integration.map_clusters()

    actual_ref_assembly_dir = read_aligner_integration.inputs["out"]

    assert (actual_ref_assembly_dir / "GCA_009556455.1_ASM955645v1_genomic.fna").exists()
    assert (actual_ref_assembly_dir / "cluster_mapped.fq").exists()
    assert (actual_ref_assembly_dir / "unmapped.fq").exists()
