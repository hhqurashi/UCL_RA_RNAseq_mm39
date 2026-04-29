# UCL RA RNA-seq mm39 Analysis

This repository contains scripts used for downstream bulk RNA-seq analysis of mouse RNA-seq data aligned to the `mm39` reference genome.

The workflow focuses on differential expression analysis, PCA visualisation, volcano plots, and pairwise comparisons between sample groups using processed RNA-seq count outputs.

## Repository contents

```text
deseq2.R
deseq2_all_pairs.R
deseq2_mm39_pairs.R
deseq2_mm39_pairs_tximport.R
deseq2_mm39_pairs_tximport_updated.R
deseq2_mm39_pairs_tximport_updated_remov...
ensmbl_vol.R
pca_HAU_controls.R
run.sh
vol.R
```

## Workflow summary

The scripts are used to:

1. Load processed RNA-seq count or Salmon/tximport outputs
2. Prepare sample metadata and group comparisons
3. Run DESeq2 differential expression analysis
4. Perform pairwise comparisons between experimental groups
5. Generate PCA plots for sample-level quality control
6. Create volcano plots for differential expression results
7. Export result tables and figures for downstream interpretation

## Main script groups

### DESeq2 analysis

```text
deseq2*.R
```

These scripts run DESeq2-based differential expression analysis, including pairwise comparisons and versions using `tximport`-processed inputs.

### PCA analysis

```text
pca_HAU_controls.R
```

This script generates PCA plots to assess sample clustering and check relationships between control/sample groups.

### Volcano plotting

```text
vol.R
ensmbl_vol.R
```

These scripts generate volcano plots from differential expression results, including versions using Ensembl-style gene identifiers.

### Run script

```text
run.sh
```

Shell script used to launch or coordinate parts of the analysis workflow.

## Input

Expected inputs may include:

- gene-level count matrices
- Salmon/tximport output files
- sample metadata
- group/condition information
- mouse gene annotation linked to the `mm39` reference

Raw FASTQ files and large upstream RNA-seq pipeline outputs are not tracked in this repository.

## Output

Outputs include:

- DESeq2 result tables
- significant gene lists
- normalised count tables
- PCA plots
- volcano plots
- comparison-specific result folders

Large result files and intermediate outputs should remain outside GitHub unless they are small, non-sensitive, and useful as examples.
