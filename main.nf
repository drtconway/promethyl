#!/usr/bin/env nextflow
/*
 * promethyl cohort pipeline
 *
 * Phase 1 (MODKIT, parallel per sample): modkit pileup + island aggregation.
 * Phase 2 (COHORT, single process):      cross-sample outlier analysis.
 */

nextflow.enable.dsl = 2

params.samples     = null                 // YAML: samples: [{id, bam}, ...]
params.annotation  = null                 // CpGs_with_promoters.bed
params.ref         = null                 // reference FASTA for modkit pileup
params.include_bed = null                 // optional BED to restrict modkit pileup
params.region      = null                 // optional region to restrict modkit pileup
params.outdir      = "results"

params.min_coverage = 10
params.mod_code     = "m"
params.min_delta    = 0.3
params.z_threshold  = 2.0
params.threads      = 10

if (!params.samples)    { error "Provide --samples cohort.yaml" }
if (!params.annotation) { error "Provide --annotation CpGs_with_promoters.bed" }

process MODKIT {
    tag "$sample_id"
    publishDir "${params.outdir}/modkit", mode: 'copy', pattern: "*_modkit.bed"
    cpus params.threads

    input:
    tuple val(sample_id), path(bam)
    path annotation
    path ref
    path include_bed

    output:
    tuple val(sample_id), path("${sample_id}.islands.tsv"), emit: islands
    path "*_modkit.bed"

    script:
    def ref_opt         = ref.name != 'NO_REF'         ? "--ref ${ref}"                 : ""
    def include_bed_opt = include_bed.name != 'NO_BED' ? "--include-bed ${include_bed}" : ""
    def region_opt      = params.region                ? "--region ${params.region}"    : ""
    """
    run_sample.py \
        --id ${sample_id} \
        --bam ${bam} \
        --modkit-dir . \
        --annotation ${annotation} \
        --output ${sample_id}.islands.tsv \
        --threads ${task.cpus} \
        --min-coverage ${params.min_coverage} \
        --mod-code ${params.mod_code} \
        ${ref_opt} ${include_bed_opt} ${region_opt}
    """
}

process COHORT {
    publishDir params.outdir, mode: 'copy'

    input:
    val sample_args
    path annotation

    output:
    path "cohort_methylation.tsv"

    script:
    """
    run_cohort.py \
        --samples ${sample_args.join(' ')} \
        --annotation ${annotation} \
        --output cohort_methylation.tsv \
        --min-delta ${params.min_delta} \
        --z-threshold ${params.z_threshold}
    """
}

workflow {
    annotation  = file(params.annotation)
    ref         = params.ref         ? file(params.ref)         : file('NO_REF')
    include_bed = params.include_bed ? file(params.include_bed) : file('NO_BED')

    cfg = new org.yaml.snakeyaml.Yaml().load(file(params.samples).text)
    samples_ch = Channel
        .fromList(cfg.samples)
        .map { s -> tuple(s.id, file(s.bam)) }

    MODKIT(samples_ch, annotation, ref, include_bed)

    // Collect all "id path" pairs into one list, then run the cohort process once.
    sample_args = MODKIT.out.islands
        .map { id, tsv -> "${id}=${tsv}" }
        .collect()

    COHORT(sample_args, annotation)
}
