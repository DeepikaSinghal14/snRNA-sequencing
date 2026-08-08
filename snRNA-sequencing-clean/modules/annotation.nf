process ANNOTATE_GENES {

    tag "Gene annotation"

    input:
    path de_results

    output:
    path "annotation_results"

    conda "${projectDir}/environment.yml"

    script:
    """
    mkdir -p annotation_results

    Rscript ${projectDir}/scripts/annotate_genes.R \
        --input ${de_results} \
        --outdir annotation_results
    """
}
