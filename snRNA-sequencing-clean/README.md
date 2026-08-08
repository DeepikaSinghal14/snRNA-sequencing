# snRNA-seq analysis of sex differences in early-life stress

A reproducible R/Nextflow workflow for single-nucleus RNA sequencing (snRNA-seq) analysis of rat spinal dorsal horn tissue from control and early-life-stress conditions.

## Biological question

This project investigates transcriptional differences associated with sex and early-life stress in the rat spinal dorsal horn, with the broader aim of understanding molecular mechanisms relevant to chronic primary low back pain.

## Workflow

```text
Cell Ranger → SoupX → CellBender
                  ↓
        Seurat object + QC
                  ↓
     Normalization + integration
                  ↓
       PCA → UMAP → clustering
                  ↓
        Cell-type annotation
                  ↓
       ┌──────────┴──────────┐
       ↓                     ↓
 Cluster markers       Sex / stress DE
                             ↓
                      Gene annotation
                             ↓
                       GO enrichment
```

## Input data

The original analysis used four CellBender-filtered matrices:

| Sample | Sex | Condition |
|---|---|---|
| `1` | Female | Control |
| `4` | Female | Stress |
| `A` | Male | Stress |
| `B` | Male | Control |

Raw sequencing data and large intermediate objects are intentionally excluded from the repository.

## Repository structure

```text
.
├── README.md
├── main.nf
├── nextflow.config
├── modules/
├── scripts/
├── conf/
├── assets/
├── data/
├── results/
└── environment.yml
```

## Run

```bash
nextflow run main.nf   --ctrl_f data/1/output_1_filtered.h5   --stress_f data/4/output_4_filtered.h5   --stress_m data/A/output_A_filtered.h5   --ctrl_m data/B/output_B_filtered.h5
```

## Outputs

The refactored workflow is designed to produce:

- integrated Seurat object
- QC and UMAP plots
- cluster marker tables
- sex-specific differential-expression tables
- stress/control differential-expression tables
- Ensembl gene annotations
- downstream enrichment results

## Reproducibility

- **Workflow orchestration:** Nextflow
- **Statistical analysis:** R / Seurat
- **Environment:** Conda
- **Version control:** Git / GitHub

## Important methodological note

The original repository contains exploratory scripts and notebooks. During refactoring, analysis decisions should be preserved deliberately rather than silently changed.

The original `DESeq2.Rmd` is labelled as pseudobulk DESeq2 but contains Seurat/SCT and `FindMarkers()` analyses. A true pseudobulk DESeq2 workflow should aggregate counts by biological replicate and specify the experimental design explicitly before being described as pseudobulk DESeq2.

## Status

Refactoring from exploratory analysis scripts into a reproducible Nextflow/R workflow.

## Author

Deepika Singhal
