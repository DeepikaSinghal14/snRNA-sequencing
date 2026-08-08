process SOUPX_QC {
    tag "$sample_id"
    label 'process_medium'
    container params.container_r

    publishDir "${params.outdir}/01_soupx_qc", mode: 'copy'

    input:
    tuple val(sample_id), path(cellranger_dir)

    output:
    tuple val(sample_id), path("${sample_id}.soupx_clean.rds"), emit: cleaned_rds

    script:
    """
    Rscript ${projectDir}/bin/01_soupx_qc.R \\
        --cellranger_dir ${cellranger_dir} \\
        --sample_id ${sample_id} \\
        --outdir .
    """

    stub:
    """
    touch ${sample_id}.soupx_clean.rds
    """
}
