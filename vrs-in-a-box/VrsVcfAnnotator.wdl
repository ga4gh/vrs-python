version 1.0

workflow VrsVcfAnnotator {
    input {
        File vcf_file
        String reference_assembly = "GRCh38"
    }

    call VrsVcfAnnotatorTask {
        input:
            vcf_file = vcf_file,
            reference_assembly = reference_assembly
    }

    output {
        File output_vcf_file = VrsVcfAnnotatorTask.output_vcf_file
        File output_vrs_objects = VrsVcfAnnotatorTask.output_vrs_objects
    }

}

task VrsVcfAnnotatorTask {
    input {
        File vcf_file
        String reference_assembly = "GRCh38"
    }

    String ref_asm_lc = reference_assembly.toLower()

    command <<<
        vrs-annotate vcf --assembly ~{reference_assembly} "~{vcf_file}" --vcf-out "with_vrs_ids.vcf" --ndjson-out "vrs_objects.json"
    >>>

    runtime {
        docker: "ga4gh/vrs-vcf-annotator-~{ref_asm_lc}:2.1.2"
        memory: "4GB"
    }

    output {
        File output_vcf_file = "with_vrs_ids.vcf"
        File output_vrs_objects = "vrs_objects.json"
    }
}
