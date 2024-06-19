#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --partition=largemem
#SBATCH -o cellranger_organoid.log

source ~/.bashrc
PROJECT_DIR="/projects/compbio/users/wma36/collaborations/JieXu"
DATA_DIR=$PROJECT_DIR/PJI23088/PJI23088-99415

declare -a elems=(
    "p22262-s009_C12N-1"
    "p22262-s010_C12N-2"
    "p22262-s001_C1-2-1"
    "p22262-s002_C1-2-2"
    "p22262-s005_IIC1-1"
    "p22262-s006_IIC1-2"
    "p22262-s011_FXS-1"
    "p22262-s012_FXS-2"
    "p22262-s003_FXSB2-1"
    "p22262-s004_FXSB2-2"
    "p22262-s008_soloc2-2"
    "p22262-s007_soloc3-1"
)

for elem in ${elems[@]}; do
    ## run cellranger count
    echo "Running 10X pipeline for $elem"
    cellranger count --id=$elem \
        --fastqs=$DATA_DIR/$elem \
        --sample=$elem \
        --transcriptome=$PROJECT_DIR/refdata-gex-GRCh38-2020-A
done
