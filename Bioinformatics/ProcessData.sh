#!/bin/bash
#!
#! Example SLURM job script for Gurdon Institute Cluster
#! Last updated: Sat Apr 18 13:05:53 BST 2015
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J PGCmeth
#! How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH -n 6
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! Job farming
#SBATCH --array=1-103%50

cd /mnt/scratch/gurdon/cap76/Thorsten/AllTics

m=( $(ls *.bam | sed -e 's/\.bam$//' ) )
FILES=( $(ls *.bam | sed -e 's/\.bam$//' ) )
a=($(echo $t | tr ',' "\n"))

T=8
wdir='/mnt/scratch/gurdon/cap76/Thorsten/UnknownTics/'
FASTA='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.fa'
GTF='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/newGTF_modified_extended5kb2.gtf'
INDICES='/mnt/scratch/surani/cap76/TEST2'
CHRSIZES='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.sizes'
GENEMODEL='/mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.bed'

mystring="cd ${wdir}fastqc_trimmed_results ; STAR --runThreadN $T --runMode alignReads --genomeDir $INDICES --readFilesIn ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1_val_1.fq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2_val_2.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix ${wdir}STAR_results/${FILES[$SLURM_ARRAY_TASK_ID-1]} --sjdbGTFfile $GTF --sjdbOverhang 149 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM"

echo $mystring

eval $mystring
