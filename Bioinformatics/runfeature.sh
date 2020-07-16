GTF='/mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf'
A=$(cat selectSamples_and_featurecounts.txt)
featureCounts -p -T 16 -t exon -g gene_id -a $GTF -o featurecountsMouse.txt $(echo $A)
