process DESEQ2_PSEUDOBULK {
    label 'process_medium'
    container params.container_r

    publishDir "${params.outdir}/03_deseq2_pseudobulk", mode: 'copy'

    input:
    path integrated_rds

    output:
    path "deseq2_pseudobulk_results.csv", emit: de_results

    script:
    """
    Rscript ${projectDir}/bin/03_deseq2_pseudobulk.R \\
        --integrated_rds ${integrated_rds} \\
        --group_var ${params.group_var} \\
        --sample_var ${params.sample_var} \\
        --outdir .
    """

    stub:
    """
    touch deseq2_pseudobulk_results.csv
    """
}
