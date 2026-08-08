process INTEGRATE_SEURAT {

    tag "Seurat integration"

    input:
    path ctrl_f
    path stress_f
    path stress_m
    path ctrl_m

    output:
    path "sex.combined.rds"
    path "integration_qc"

    conda "${projectDir}/environment.yml"

    script:
    """
    mkdir -p integration_qc

    Rscript ${projectDir}/scripts/integrate_seurat.R \
        --ctrl-f ${ctrl_f} \
        --stress-f ${stress_f} \
        --stress-m ${stress_m} \
        --ctrl-m ${ctrl_m} \
        --out-rds sex.combined.rds \
        --outdir integration_qc
    """
}
