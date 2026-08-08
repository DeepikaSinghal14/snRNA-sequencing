#!/usr/bin/env nextflow
/*
 * snRNA-sequencing pipeline
 * Sex differences in differentially expressed genes in early-life stress
 * induced chronic primary low back pain in rats.
 *
 * Stages:
 *   1. SOUPX_QC            - ambient RNA removal + Seurat object per sample
 *   2. SEURAT_INTEGRATION  - CellBender-denoised counts, integrated across samples
 *   3. DESEQ2_PSEUDOBULK   - pseudobulk differential expression (stress vs control)
 *   4. ENSEMBL_ANNOTATION  - annotate DE results with Ensembl gene metadata
 */

nextflow.enable.dsl = 2

include { SOUPX_QC }            from './modules/soupx_qc.nf'
include { SEURAT_INTEGRATION }  from './modules/seurat_integration.nf'
include { DESEQ2_PSEUDOBULK }   from './modules/deseq2_pseudobulk.nf'
include { ENSEMBL_ANNOTATION }  from './modules/ensembl_annotation.nf'

workflow {

    // --- Step 1: per-sample ambient RNA removal --------------------------
    // samplesheet columns: sample_id,cellranger_dir,cellbender_h5,sex,group
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.cellranger_dir)) }
        .set { cellranger_ch }

    SOUPX_QC(cellranger_ch)

    // --- Step 2: integration across all samples ---------------------------
    // Integration reads CellBender h5 paths directly from the samplesheet,
    // independent of the SoupX branch above (see README for why these two
    // QC paths are kept separate, mirroring the original scripts).
    cellbender_files_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> file(row.cellbender_h5) }
        .collect()

    SEURAT_INTEGRATION(
        file(params.samplesheet),
        cellbender_files_ch
    )

    // --- Step 3: pseudobulk DE ---------------------------------------------
    DESEQ2_PSEUDOBULK(SEURAT_INTEGRATION.out.integrated_rds)

    // --- Step 4: gene annotation --------------------------------------------
    ENSEMBL_ANNOTATION(DESEQ2_PSEUDOBULK.out.de_results)
}

workflow.onComplete {
    log.info "Pipeline complete -- success: ${workflow.success}, duration: ${workflow.duration}, results: ${params.outdir}"
}
