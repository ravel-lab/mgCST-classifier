#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#args<-c("~/Dropbox (IGS)/mgss_mcst_devel/mgCST\ classifier/VIRGO-master/mgCST_classifier/LSVF/summary.Abundance.txt", "~/Dropbox (IGS)/mgss_mcst_devel/mgCST\ classifier/VIRGO-master/mgCST_classifier/LSVF/summary.NR.abundance.txt", "~/Dropbox (IGS)/mgss_mcst_devel/mgCST\ classifier/VIRGO-master/", "~/Dropbox (IGS)/MSL/VALENCIA-master/mgCST_centroids_22Jun2022.csv")
wd<-getwd()

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Full paths to the following files/directories must be supplied in this order:\n(1) summary.Abundance.txt\n(2) summary.NR.abundance.txt\n(3) VIRGO-master directory\n(4) reference centroids", call.=FALSE)
}

mgss.classifiers<-readRDS(paste(args[3], "/mgCST_classifier/mgss_classifier.RDS", sep=""))
if(!exists("mgss.classifiers")){
  print("Path to VIRGO-master not correct")
  stop()
}

# Capture date for output
today <- strsplit(date(), " ")
month <- today[[1]][2]
day <- today[[1]][3]
year <- today[[1]][5]
today2 <- paste(day,month,year, sep="")

## read in data
counts.mgss<-as.data.frame(read.delim(args[1], header = TRUE, row.names = 1, as.is = TRUE))
## Update bvab1 to Ca._Lachnocurva vaginae
bvab1.col<-which(names(counts.mgss) %in% "BVAB1")
names(counts.mgss)[bvab1.col]<-"Ca._Lachnocurva_vaginae"
counts.mgss.ngl<-counts.mgss
counts.genes<-as.data.frame(read.delim(args[2], header = TRUE, check.names = FALSE))
gene.length<-as.data.frame(read.delim(paste(args[3], "/1_VIRGO/", "0.geneLength.txt", sep=""), header = FALSE))
names(gene.length)<-c("PC","Length")
taxon.tbl<-as.data.frame(read.delim(paste(args[3], "/1_VIRGO/", "1.taxon.tbl_Jun_20_2022.txt", sep=""), header = FALSE))
names(taxon.tbl)<-c("cluster", "PC", "Taxa", "Length")

## format gene table
genes<-merge(gene.length, taxon.tbl[,c("PC", "Taxa")], all=TRUE, by="PC")
genes<-merge(genes, counts.genes, all.y=TRUE, by="PC")

a<-which(names(genes) %in% "Taxa")
a<-a+1
genes.ngl<-genes
genes.ngl[,a:ncol(genes.ngl)]<-apply(genes[,a:ncol(genes)], 2, function(x) x*150/genes[["Length"]])
write.csv(genes.ngl, paste(wd, "/norm_counts_", today2, ".csv", sep=""), quote=F)

## Make gene table presence/absence
genes.ngl.pa<-genes
genes.ngl.pa[,a:ncol(genes.ngl.pa)]<-apply(genes[,a:ncol(genes)], 2, function(x) ifelse(x*150/genes[["Length"]] >= 0.5, 1, 0))


if (!require(randomForestSRC, quietly = TRUE)) install.packages("randomForestSRC", quiet = TRUE)
require(randomForestSRC)

## classify mgss
for (taxon in names(counts.mgss[names(counts.mgss) %in% names(mgss.classifiers)])){
  ## reformat taxon-specific gene count table to match classifier
  table<-as.data.frame(t(genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, a:ncol(genes.ngl.pa)]))
  names(table)<-genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, 1]
  ngl.sum<-rowSums(as.data.frame(t(genes.ngl[genes.ngl[["Taxa"]] %in% taxon, a:ncol(genes.ngl)])))
  ## rename base taxon to be taxon_0
  names(counts.mgss.ngl)[names(counts.mgss.ngl) %in% taxon]<-paste(taxon, 0, sep="_")
  ## capture samples with too few species genes to qualify for mgss
  samples.for.0<-rownames(table)[rowSums(table) < 750]
  ## make gene count table match format of classifier
  if(length(samples.for.0) < nrow(table)){ ## if there are any samples with enough genes to classify mgss:
    for(gene in mgss.classifiers[[taxon]]$xvar.names){
      if(is.null(table[[gene]])){
        table[[gene]]<-0
      }
    }
    ## reorder gene count table columns  match format of classifier
    table<-table[,mgss.classifiers[[taxon]]$xvar.names]
    ## predict mgss and get the estimated class probability for each class (mgss).
    predicted.mgss<-predict.rfsrc(mgss.classifiers[[taxon]], table)$predicted
    ss<-as.data.frame(cbind(mgss=colnames(predicted.mgss)[apply(predicted.mgss,1,which.max)]), row.names = rownames(table))
    for (sample in rownames(table)){
      ## transfer counts to taxon_mgss (as long as it has enough species genes for mgss)
      new<-ifelse(sample %in% samples.for.0, paste(taxon, 0, sep="_"), paste(taxon, ss[sample, "mgss"], sep="_"))
      if(!new %in% paste(taxon, 0, sep="_")){
        counts.mgss.ngl[sample,new]<-as.numeric(ngl.sum[sample])
        counts.mgss.ngl[sample,paste(taxon, 0, sep="_")]<-0
      }
    }
  }
  for (sample in samples.for.0){
    counts.mgss.ngl[sample,paste(taxon, 0, sep="_")]<-as.numeric(ngl.sum[sample])
  }
}
## Transfer all none mgss taxon counts to ngl (normalized by gene length)
for (taxon in names(counts.mgss[!names(counts.mgss) %in% names(mgss.classifiers)])){
  ngl.sum<-rowSums(as.data.frame(t(genes.ngl[genes.ngl[["Taxa"]] %in% taxon, a:ncol(genes.ngl)])))
  table<-as.data.frame(t(genes.ngl.pa[genes.ngl.pa[["Taxa"]] %in% taxon, a:ncol(genes.ngl.pa)]))
  for (sample in rownames(table)){
    counts.mgss.ngl[sample,taxon]<-as.numeric(ngl.sum[sample])
  }
}
## make relative abundance table
counts.mgss.ngl[is.na(counts.mgss.ngl)]<-0
write.csv(counts.mgss.ngl, paste(wd, "/norm_counts_mgSs_mgCST_", today2, ".csv", sep=""), quote=F)

## ## Centroid classifier
#defining function to determine yue-clayton theta
yue_distance<-function(row, median){
  #creating a counting variable to index the median list
  taxon_count = 1
  #creating lists to iteratively store output
  median_times_obs <- vector()
  median_minus_obs_sq <- vector()
  #looping through the row and calculating product and difference squared between row data and median data
  for (taxon_abund in row) {
    #calculate p * q
    median_times_obs<-append(median_times_obs, as.numeric(median[taxon_count])*taxon_abund)
    #calculate p-q squared
    median_minus_obs_sq<-append(median_minus_obs_sq, as.numeric((median[taxon_count]-taxon_abund)**2))
    taxon_count <- taxon_count+1
  }
  #calculate sum p* q
  product = sum(median_times_obs)
  #calculate sum p-q squared
  diff_sq = sum(median_minus_obs_sq)
  #calculate yue_med_dist
  yue_med_dist = product / (diff_sq + product)
  #return the value of yue distance
  return(yue_med_dist)
}

## read in reference centroids
reference_centroids<-as.data.frame(read.csv(args[4], header = TRUE, row.names = 1))
## make relabund table fresh
relabund<-counts.mgss.ngl/rowSums(counts.mgss.ngl)
## reformat relabund to include all expected column names (xvar.names)
relabund<-relabund[, names(relabund) %in% names(reference_centroids)]
for(taxa in names(reference_centroids)){
  if(is.null(relabund[[taxa]])){
    relabund[[taxa]]<-0
  }
}
relabund<-relabund[,names(reference_centroids)]
n<-ncol(relabund)
## for each mgCST measure the similarity of each sample to each mgCST centroid using yue + clayton theta
for(i in 1:27){
  mgCST<-paste("mgCST", i, sep=" ")
  relabund[[mgCST]]<-apply(relabund[,1:n], 1, function(x) yue_distance(x, reference_centroids[mgCST,]))
}
m<-n+1
relabund[["mgCST"]]<-colnames(relabund[,m:which(colnames(relabund) %in% "mgCST 27")])[apply(relabund[,m:which(colnames(relabund) %in% "mgCST 27")],1,which.max)]
write.csv(relabund, paste(wd, "/relabund_w_mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)
write.csv(relabund["mgCST"], paste(wd, "/mgCSTs_", today2, ".csv", sep=""), row.names = TRUE, quote=F)

## plot heatmap
mgCST<-as.data.frame(rbind(c("1", "#FE0308"),c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),c("5", "#F0BCCC"),c("6", "#F6D3DA"),c("7", "#86C61A"), c("8", "#B4DB29"),c("9", "#DBEA77"), c("10", "#FF7200"), c("11", "#F68A11"),c("12", "#F8A40E"),c("13", "#F3BC11"),c("14", "#f7d15a"), c("15", "#FAE50D"),c("16", "#F3F46E"),c("17", "#448A73"),c("18", "#89BEAB"), c("19", "#BCD6CD"),c("20", "#221886"),c("21", "#3E3792"),c("22", "#5D579E"),c("23", "#7C76AC"),c("24", "#9A98BF"),c("25", "#C9C8D8"),c("26", "#98C999"), c("27", "#989898"), c("", "white"), c("NA", "white")))
names(mgCST)<-c("mgCST", "color")
colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
relabund.mgCST<-relabund[,1:n]
relabund.mgCST<-relabund.mgCST[,order(colSums(relabund.mgCST), decreasing = TRUE)]
relabund.mgCST$mgCST<-gsub("mgCST ", "", relabund$mgCST)
relabund.mgCST<-relabund.mgCST[order(as.numeric(relabund.mgCST[["mgCST"]])),]
relabund.mgCST[["color"]]<-mgCST[match(relabund.mgCST[["mgCST"]], mgCST$mgCST), "color"]
names(relabund.mgCST)<-gsub("_", " ", names(relabund.mgCST))

if (!require(gplots, quietly = TRUE, warn.conflicts = FALSE)) install.packages("gplots", quiet = TRUE)
pdf(paste(wd, "/mgCST_heatmap_", today2, ".pdf", sep=""), width=7, height=10)
gplots::heatmap.2(t(as.matrix(relabund.mgCST[,1:50])), Colv = FALSE, Rowv = FALSE, col=colfunc(100), keysize= 1.0, densadj=0, density.info='none', key = TRUE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.7, cexCol = 0.1, adjRow = c(1, NA),offsetRow = -38,  main = paste("mgCST Heatmap\nnSamples=", paste(nrow(relabund))), title(main = paste("mgCST Heatma\nnSamples=", paste(nrow(relabund))), line = -2), ColSideColors = as.vector(relabund.mgCST[["color"]]), lhei = c(1,7), dendrogram = "none")
dev.off()

