// Main script to parse through data file and update html report.
// DataTables: https://datatables.net/examples/data_sources/js_array.html

$(document).ready(function () {
    'use strict';

    // Eventually useful for multi-file upload
    var fileread_flag = 0;
    var numFiles = 0;


    $(document).ready(function () {

        // Upon File upload this function is called
        //
        // TODO: Implement multiple file upload handling.
        // Right now, can read-in multiple files but display 
        // functions can't parse through them.
        //
        let fileUpload = document.getElementById("inputfile");
        fileUpload.onchange = function (event) {
            numFiles = event.target.files.length;
            var i;
            for (i = 0; i < numFiles; i++) {
                var reader = new FileReader();
                reader.onload = onReaderLoad;
                reader.readAsText(event.target.files[i]);
            }
        }

    });

    //Function to read the files uploaded
    function onReaderLoad(event) {

        var reportInfo = JSON.parse(event.target.result);
        display_header(reportInfo.header);
        display_ref_selection(reportInfo.reference_selection);
        display_ref_culling(reportInfo.reference_culling);
        display_ref_assembly(reportInfo.ref_assembly);
        display_denovo_assembly(reportInfo.denovo_assembly);
        display_merged_assembly(reportInfo.merged_assembly);

        fileread_flag += 1;

        if (fileread_flag == numFiles) {
            document.getElementById('upload_msg').innerHTML = "File reading completed!";
        }
    }

    // Alter header elements all by id.
    // To add to header, must add to data file, and added to header with appropriate ids
    function display_header(header) {
        document.getElementById('forward').innerHTML = header.reads.forward;
        document.getElementById('reverse').innerHTML = header.reads.reverse;
        document.getElementById('unpaired').innerHTML = header.reads.unpaired;

        document.getElementById('num_forward').innerHTML = header.summary.num_forward;
        document.getElementById('num_reverse').innerHTML = header.summary.num_reverse;
        document.getElementById('num_unpaired').innerHTML = header.summary.num_unpaired;

        document.getElementById('ref_db').innerHTML = header.parameters.ref_db;
        document.getElementById('ref_filtering').innerHTML = header.parameters.ref_filtering;
        document.getElementById('ref_filtering_kmer_size').innerHTML = header.parameters.ref_filtering_kmer_size;
        document.getElementById('depth_of_cov').innerHTML = header.parameters.depth_of_cov;
        document.getElementById('breadth_of_cov').innerHTML = header.parameters.breadth_of_cov;
        document.getElementById('ref_culling_kmer_size').innerHTML = header.parameters.ref_culling_kmer_size;
    }

    function display_ref_selection(info) {
        var res = '';
        var num_refs = info.num_selected;
        var refs = info.references;
        var ref_info = [];

        document.getElementById('reference_selection-info').innerHTML = "Number of Selected References: " + num_refs;

        // parsing data to fit format for datatables
        refs.forEach(function (ref) {
            var curr_ref = [];
            curr_ref.push(ref.id);
            curr_ref.push(ref.name);
            curr_ref.push(ref.length);
            curr_ref.push(ref.num_seqs_covered);
            curr_ref.push(ref.depth);
            curr_ref.push(ref.breadth);
            ref_info.push(curr_ref);
        });

        // populating datatable
        $('#referenceSelection').DataTable({
            data: ref_info,
            columns: [
                { title: "Reference ID" },
                { title: "Scientific Name" },
                { title: "Genome Length" },
                { title: "# Sequences Mapped" },
                { title: "Depth of Cov." },
                { title: "Breadth of Cov." }
            ]
        });
    }

    function display_ref_culling(info) {
        var res = '';
        var num_refs = info.num_culled;
        var refs = info.references;
        var ref_info = [];

        document.getElementById('reference_culling-info').innerHTML = "Number of References Culled From Selection Step: " + num_refs;

        // parsing data to fit format for datatables
        refs.forEach(function (ref) {
            var curr_ref = [];
            curr_ref.push(ref.id);
            curr_ref.push(ref.name);
            curr_ref.push(ref.length);
            curr_ref.push(ref.num_seqs_covered);
            curr_ref.push(ref.depth);
            curr_ref.push(ref.breadth);
            ref_info.push(curr_ref);
        });

        // populating datatable
        $('#referenceCulling').DataTable({
            data: ref_info,
            columns: [
                { title: "Reference ID" },
                { title: "Scientific Name" },
                { title: "Genome Length" },
                { title: "# Sequences Mapped" },
                { title: "Depth of Cov." },
                { title: "Breadth of Cov." }
            ]
        });
    }

    function display_ref_assembly(info) {
        var res = '';
        var refs = info.references;
        var ref_info = [];

        document.getElementById('num_refs').innerHTML = info.summary.num_refs;
        document.getElementById('total_num_contigs').innerHTML = info.summary.total_num_contigs;
        document.getElementById('total_length_of_assembly').innerHTML = info.summary.total_length_of_assembly;
        document.getElementById('max_contig_length').innerHTML = info.summary.max_contig_length;
        document.getElementById('min_contig_length').innerHTML = info.summary.min_contig_length;

        // parsing data to fit format for datatables
        refs.forEach(function (ref) {
            var curr_ref = [];
            curr_ref.push(ref.id);
            curr_ref.push(ref.name);
            curr_ref.push(ref.length);
            curr_ref.push(ref.num_contigs);
            curr_ref.push(ref.assembly_length);
            curr_ref.push(ref.max_contig_length);
            curr_ref.push(ref.min_contig_length);
            curr_ref.push(ref.ng50);
            ref_info.push(curr_ref);
        });

        // populating datatable
        $('#referenceAssembly').DataTable({
            data: ref_info,
            columns: [
                { title: "Reference ID" },
                { title: "Scientific Name" },
                { title: "Genome Length" },
                { title: "Number of Contigs" },
                { title: "Length of Assembly" },
                { title: "Max Contig Length" },
                { title: "Min Contig Length" },
                { title: "NG50" }
            ]
        });
    }

    function display_denovo_assembly(info) {
        var res = '';

        document.getElementById('denovo_total_num_contigs').innerHTML = info.total_num_contigs;
        document.getElementById('denovo_total_length_of_assembly').innerHTML = info.total_length_of_assembly;
        document.getElementById('denovo_max_contig_length').innerHTML = info.max_contig_length;
        document.getElementById('denovo_min_contig_length').innerHTML = info.min_contig_length;

    }

    function display_merged_assembly(info) {
        var res = '';

        document.getElementById('merged_total_num_contigs').innerHTML = info.total_num_contigs;
        document.getElementById('merged_total_length_of_assembly').innerHTML = info.total_length_of_assembly;
        document.getElementById('merged_max_contig_length').innerHTML = info.max_contig_length;
        document.getElementById('merged_min_contig_length').innerHTML = info.min_contig_length;

    }
});