# Bioinformatics processing.

1) Data was trimmed and aligned by running ProcessData.sh script with appropriate modifications to file paths

2) Gene counts were calculated by running the runfeature.sh script.

3) Processed data is icluded in the Data folder. Unzip data files and data anotation keys for further analysis.

4) Quality control metrics from (1) were aggregated in a summary spreadsheet. Appropriate plots of QC metrics can be plotted using the QC.R script.

5) Dimensionality reduction on data could be run using the DimRed.R script.

6) Heatmaps were generated using the pltHeatmap.R script

7) Seurat alignment of marmoset, cynomolgus, and human dataset could be done using Cross_species_alignment.R script
