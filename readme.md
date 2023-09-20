# Leukemia outlier
Aberrant genome and transcriptome analysis of 3,760 hematologic malignancies.



## OUTRIDER
OUTRIDER was run using branch `outrier2` (https://github.com/gagneurlab/OUTRIDER/tree/outrider2) and DROP v 1.1.0 (https://github.com/gagneurlab/drop/releases/tag/1.1.0)



## FRASER
FRASER was run using v 1.99.0 (https://github.com/c-mertes/FRASER/releases/tag/1.99.0) and DROP v 1.2.3 (https://github.com/gagneurlab/drop/releases/tag/1.2.3)



## AbSplice
AbSplice was run using v 1.0.0 (https://github.com/gagneurlab/absplice)



## intOGen
intOGen was run on version commit 437a047 (https://bitbucket.org/intogen/intogen-plus/src/437a0473b3ccb187415b6b2fc68c5b7fbe15ae3e/)

### Run intOGen
1. Move to intOGen folder 
`cd intogen-plus`
2. Download file 
`wget https://bitbucket.org/intogen/intogen-plus/src/437a0473b3ccb187415b6b2fc68c5b7fbe15ae3e/build/containers/dndscv/dndscv.tar.gz`
3. Move file 
`mv dndscv.tar.gz build/containers/dndscv/`
4. Run intOGen as `README.md`



## Driver gene prediction
Pipeline to identify leukemia driver gene.

### Setup project
1. Specify project directory path `projectPath` in `wbuild.yaml` and `Scripts/experiment_design_mll.py`
2. Install conda environment `conda env create -f driver_gene_prediction/envs/vale_202204.yaml`

### Run driver gene prediction
1. Move to driver gene prediction folder 
`cd driver_gene_prediction`
2. Add experimental design by running `Scripts/experiment_design.py`
3. Dry run: `snakemake -n`
4. Run the pipeline: `snakemake --cores 4 --use-conda`
