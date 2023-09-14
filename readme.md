# Leukemia outlier

Aberrant transcriptome analysis of 3,760 hematologic malignancies.

Pipeline to identify leukemia driver gene.


## Setup
`projectPath` needs to be set in `wbuild.yaml` and `Scripts/experiment_design_mll.py`

conda env: `vale_202204` can be installed using `envs/vale_202204.yaml`


## Run the pipeline
1. Add experimental design using `Scripts/experiment_design.py`
2. Dry run: `snakemake -n`
3. Run the pipeline: `snakemake --cores 4 --use-conda`
