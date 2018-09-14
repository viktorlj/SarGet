#!/bin/bash
#SBATCH -A projID
#SBATCH -p core
#SBATCH -t 12:00:00
##SBATCH --qos=short
##SBATCH -t 00:14:00
set -x
# skeleton script to launch nextflow/singularity Sarek jobs with slurm on irma
#
# PARAMETERS:
#  add $1 as -resume to resume an old run
#
SARGETDIR=~/Sarget
DATE=`date +%Y-%b-%d-%H%M`
PREFIX=UMI_calling_${DATE}

# nextflow specific stuff to save everything on /scratch
module load bioinfo-tools Nextflow
export NXF_SINGULARITY_CACHEDIR=`pwd`
export PROJECT=projID
export NXF_LAUNCHBASE=/scratch
export NXF_WORK=`pwd`/work
#export NXF_WORK=/scratch
export NXF_HOME=/sw/apps/bioinfo/Nextflow/0.31.1/bianca/nxf_home

export CONTAINERPATH=~/Sarget/containers/
export GENOMEBASE=/proj/projID/nobackup/SarGet/reference
export NXF_LAUNCHER=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP

#mkdir ${1} 
#cd ${1}

nextflow run ${SARGETDIR}/main.nf \
               --project $PROJECT \
               --sample FullRun.tsv \
               --ampliconFile AgilentPanelID_Amplicons.bed \
               --regionsFile AgilentPanelID_Regions.bed \
               -with-report ${PREFIX}.report.html \
               --containerPath ${CONTAINERPATH} \
               --genome GRCh37 \
               --genome_base ${GENOMEBASE} \
               --outDir `pwd`/results \
               -profile slurm ${1}