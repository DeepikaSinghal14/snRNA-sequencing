process ENSEMBL_ANNOTATION {
    label 'process_low'
    container params.container_r

    publishDir "${params.outdir}/04_ensembl_annotation", mode: 'copy'

    input:
    path de_results_csv

    output:
    path "de_results_annotated.csv", emit: annotated_results

    script:
    """
    Rscript ${projectDir}/bin/04_ensembl_annotation.R \\
        --de_results_csv ${de_results_csv} \\
        --species_dataset ${params.species_dataset} \\
        --outdir .
    """

    stub:
    """
    touch de_results_annotated.csv
    """
}
