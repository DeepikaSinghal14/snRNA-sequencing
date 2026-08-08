nextflow.enable.dsl=2

params.ctrl_f   = null
params.stress_f = null
params.stress_m = null
params.ctrl_m   = null

include { INTEGRATE_SEURAT } from './modules/integrate'
include { DIFFERENTIAL_EXPRESSION } from './modules/differential_expression'
include { ANNOTATE_GENES } from './modules/annotation'

workflow {
    if (!params.ctrl_f || !params.stress_f || !params.stress_m || !params.ctrl_m) {
        error 'Provide --ctrl_f, --stress_f, --stress_m and --ctrl_m'
    }

    ctrl_f   = channel.fromPath(params.ctrl_f, checkIfExists: true)
    stress_f = channel.fromPath(params.stress_f, checkIfExists: true)
    stress_m = channel.fromPath(params.stress_m, checkIfExists: true)
    ctrl_m   = channel.fromPath(params.ctrl_m, checkIfExists: true)

    integrated = INTEGRATE_SEURAT(ctrl_f, stress_f, stress_m, ctrl_m)
    de_results = DIFFERENTIAL_EXPRESSION(integrated)
    ANNOTATE_GENES(de_results)
}
