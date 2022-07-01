# mgCST-classifier
When used in concert with the vaginal non-redundant gene database (VIRGO), assigns metagenomic community state types (mgCSTs) to vaginal metagenomes. 
This script *requires* VIRGO output files "summary.Abundance.txt" and "summary.NR.abundance.txt" as well as VIRGO reference file 0.geneLength.txt.

To Use: 
  1. Clone this repository
  2. Download the metagenomic subspecies random forest classifier from https://figshare.com/ndownloader/files/36123449 and place in the github repository directory. 
  3. Run data through VIRGO.
  4. Run classifier. The classifier will first assign metagenomic subspecies using the mgss_classifier.RDS file downloaded from figshare. It will then assign metagenomic CSTs to each sample using a centroid-based classifier. 
  
Rscript classify_mgCST_centroid.R

Full paths to the following files/directories must be supplied in this order:
(1) summary.Abundance.txt
(2) summary.NR.abundance.txt
(3) VIRGO-master directory
(4) mgCST-classifier-master directory

Example: Rscript classify_mgCST_centroid.R summary.Abundance.txt summary.NR.abundance.txt /full/path/to/VIRGO-master /full/path/to/mgCST-classifier-master
