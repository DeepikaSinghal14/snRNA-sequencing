process DIFFERENTIAL_EXPRESSION {

    tag "Differential expression"

    input:
    path seurat_object

    output:
    path "de_results"

    conda "${projectDir}/environment.yml"

    script:
    """
    mkdir -p de_results

    Rscript ${projectDir}/scripts/differential_expression.R \
        --input ${seurat_object} \
        --outdir de_results
    """
}
