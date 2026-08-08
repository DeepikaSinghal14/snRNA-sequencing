process SEURAT_INTEGRATION {
    label 'process_high'
    container params.container_r

    publishDir "${params.outdir}/02_seurat_integration", mode: 'copy'

    input:
    path samplesheet
    path cellbender_h5_files

    output:
    path "sex.combined.rds", emit: integrated_rds

    script:
    """
    Rscript ${projectDir}/bin/02_seurat_integration.R \\
        --samplesheet ${samplesheet} \\
        --min_features ${params.min_features} \\
        --nfeatures ${params.nfeatures} \\
        --outdir .
    """

    stub:
    """
    touch sex.combined.rds
    """
}
