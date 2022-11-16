# mgCST-classifier
When used in concert with the vaginal non-redundant gene database (VIRGO), assigns metagenomic community state types (mgCSTs) to vaginal metagenomes. 
This script *requires* VIRGO output files "summary.Abundance.txt" and "summary.NR.abundance.txt" as well as VIRGO reference file 0.geneLength.txt.

To Use: 
  1. Clone this repository
  2. Download the metagenomic subspecies random forest classifier from https://figshare.com/ndownloader/files/36123449 and place in the github repository directory. 
  3. Run data through VIRGO (virgo.igs.umaryland.edu).
  4. Run classifier. The classifier will first assign metagenomic subspecies using the mgss_classifier.RDS file downloaded from figshare. It will then assign metagenomic CSTs to each sample using a nearest-neighbor, centroid-based classifier. 

Full paths to the following files/directories must be supplied in this order:
(1) summary.Abundance.txt
(2) summary.NR.abundance.txt
(3) VIRGO-master directory
(4) mgCST-classifier-master directory

Example: 

Rscript classify_mgCST_centroid.R   example/summary.Abundance.txt   example/summary.NR.abundance.txt   /full/path/to/VIRGO-master   /full/path/to/mgCST-classifier-master

Output is written to current directory. Each output file is dated.

Output files: 

*relabund_w_mgCSTs_1Jul2022.csv*: Relative abundances of metagenomic subspecies (columns) by samples (rows). Right-most column contains mgCST classification.

*mgCST_heatmap_1Jul2022.pdf*: A heatmap of relative abundances of 50 most abundant metagenomic subspecies with samples colored by mgCST. 

*mgCSTs_1Jul2022.csv*: Two-column table with sample ID and mgCST assignment

*norm_counts_mgSs_mgCST_1Jul2022.csv*: Gene abundances normalized by gene length of metagenomic subspecies (columns) by samples (rows). Right-most column contains mgCST classification. Normalization assumed 150 bp long reads: (Gene_Abundance*150)/Gene_Length
