## Driver gene prediction
Pipeline to identify leukemia driver gene.

### Setup project
1. Specify project directory path `projectPath` in `wbuild.yaml` and `Scripts/experiment_design_mll.py`
2. Install conda environment `conda env create -f driver_gene_prediction/envs/vale_202204.yaml`

### Run driver gene prediction
1. Add experimental design by running `Scripts/experiment_design.py`
2. Dry run: `snakemake -n`
3. Run the pipeline: `snakemake --cores 4 --use-conda`
