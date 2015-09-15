library(rtracklayer)
library(genomation)
library(GenomicFeatures)
library(parallel)

### PATHS
data_dir <- "/data/akalin/kasia/data/"
gene_annotation_path = "/data/akalin/kasia/data/gencode.v19.annotation.gtf.gz"
lincRNA_annotation_path = paste0(data_dir, "gencode.v19.long_noncoding_RNAs.gtf")
metadata_path = "/data/akalin/kasia/data/metadata.tsv"

### FUNCTIONS

readBigWig = function(target, windows=NULL, ...){
  ## read BigWig file
  if(!is.null(windows) & class(windows) != 'GRanges')
    stop('windows argument needs to be a GRanges object')
  if(is.null(windows)){
    bw = import(target, asRangedData = FALSE)
  }else{
    bw = import(target, asRangedData = FALSE, which=windows)
  }
  if(length(bw) == 0)
    stop('There are no ranges selected')
  
  covs = coverage(bw, weight=bw$score)
  return(covs)
}

binned.arithmfunc <- function(bins, numvar, mcolname, func)
{
  stopifnot(is(bins, "GRanges"))

  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         func(views)
                       })
  mcols(bins)[[mcolname]] <- unlist(means_list)[names(bins)]
  bins
}


### MAIN

# read lincRNAs
lincRNAs <- gffToGRanges(lincRNA_annotation_path)
lincRNAs <- unlist(lincRNAs)
lincRNAs = lincRNAs[lincRNAs$type=="transcript"]
names(lincRNAs) = lincRNAs$transcript_id	      

# read introns and exons
txdb <- makeTxDbFromGFF(gene_annotation_path, format="gtf", circ_seqs=character())
introns <- intronsByTranscript(txdb, use.names=TRUE)
exons <- exonsBy(txdb, use.names=TRUE)

a <- mclapply(introns[1:100],function(ranges){if(length(ranges)) ranges$idx=1:length(ranges);ranges}, mc.cores=30)
introns <- Reduce('c', a)

exons <- unlist(exons)
exons$transcript_id <- names(exons)


#read info about samples
metadata <- read.table(metadata_path, 
			sep="\t",
			header=TRUE)	
metadata$strand = "dummy"
metadata[which(metadata$Output.type=="plus strand signal of unique reads"),]$strand = "+"
metadata[which(metadata$Output.type=="minus strand signal of unique reads"),]$strand = "-"
metadata[which(metadata$Output.type=="plus strand signal of all reads"),]$strand = "+"
metadata[which(metadata$Output.type=="minus strand signal of all reads"),]$strand = "-"
metadata$amount = "dummy"
metadata[which(metadata$Output.type=="plus strand signal of unique reads"),]$amount = "unique"
metadata[which(metadata$Output.type=="minus strand signal of unique reads"),]$amount = "unique"
metadata[which(metadata$Output.type=="plus strand signal of all reads"),]$amount = "all"
metadata[which(metadata$Output.type=="minus strand signal of all reads"),]$amount = "all"

#only uniquely mapped reads
metadata = metadata[which(metadata$amount=="unique"),]
# only aortic smooth muscle cell
metadata = metadata[which(metadata$Biosample.term.name=="aortic smooth muscle cell"),]

metadata$id <- 1:nrow(metadata)


### calculate sum of coverage within introns, exons, lncRNAs

file_acc = as.character(metadata$File.accession)
exp_acc = as.character(metadata$Experiment.accession)
bw_files = paste0(data_dir, file_acc, ".bigWig")

readBigWig.parallel <- function(i, targets=targets){

    readBigWig(targets[i])
}    
covs <- mclapply(1:nrow(metadata), 
                  readBigWig.parallel,
                  targets=bw_files,
                  mc.cores=30
                  )
names(covs) <- file_acc

normalize = function(cov, nmbr){
  # for each chromosome separately
  for(i in 1:length(cov)){
    cov[[i]] <- cov[[i]] / nmbr  
  }
  cov
}
# normalize coverage of samples on + and - strand by sum of both of them
aa <- unique(metadata$Derived.from)
for(i in 1:length(aa)){
   print(i)
   samples_both_strands <- metadata [ metadata$Derived.from==aa[i] ,]
   sum_covs <- sum(sapply(covs[samples_both_strands$id], sum))
   covs[[samples_both_strands$id[1]]] <- normalize(covs[[samples_both_strands$id[1]]], sum_covs)
   covs[[samples_both_strands$id[2]]] <- normalize(covs[[samples_both_strands$id[2]]], sum_covs)
}

      
for(i in 1:length(metadata$File.accession)){

  introns <- binned.arithmfunc(introns, covs[[i]], names(covs)[i], viewSums)
  exons <- binned.arithmfunc(exons, covs[[i]], names(covs)[i], viewSums)
  lincRNAs <- binned.arithmfunc(lincRNAs, covs[[i]], names(covs)[i], viewSums)
  
  #normalization; normalize by sum within sample
  mcols(introns)[,names(covs)[i]] <- mcols(introns)[,names(covs)[i]]*1. / sum(mcols(introns)[,names(covs)[i]])
  mcols(exons)[,names(covs)[i]] <- mcols(exons)[,names(covs)[i]]*1. / sum(mcols(exons)[,names(covs)[i]])
  mcols(lincRNAs)[,names(covs)[i]] <- mcols(lincRNAs)[,names(covs)[i]]*1. / sum(mcols(lincRNAs)[,names(covs)[i]])

}

introns_ave <- data.frame(1:length(introns))
exons_ave <- data.frame(1:length(exons))
lincRNA_ave <- data.frame(1:length(lincRNAs))

aa <- as.character(unique(metadata$Derived.from))
for(i in 1:length(aa)){
   samples_both_strands <- as.character(metadata [ metadata$Derived.from==aa[i] ,]$File.accession)
   x1 <- rowSums(as.data.frame(mcols(introns)[samples_both_strands]), na.rm = TRUE)
   y1 <- rowSums(as.data.frame(mcols(exons)[samples_both_strands]), na.rm = TRUE)
   z1 <- rowSums(as.data.frame(mcols(lincRNAs)[samples_both_strands]), na.rm = TRUE)
   introns_ave <- cbind(introns_ave, x1)
   exons_ave <- cbind(exons_ave, y1)
   lincRNA_ave <- cbind(lincRNA_ave, z1)
}

introns$average = rowMeans(introns_ave)
exons$average = rowMeans(exons_ave)
lincRNAs$average = rowMeans(lincRNA_ave)

# take average of replicates
df.introns = data.frame(txid=introns$transcript_id, type="intron", idx=introns$intron_rank, expression=introns$average)
df.exons = data.frame(txid=exons$transcript_id, type="exon", idx=exons$exon_rank, expression=exons$average)
df.lincRNA = data.frame(txid=lincRNAs$transcript_id, type="lincRNA", idx=0, expression=lincRNAs$average)

df <- rbind(df.introns, df.exons, df.lincRNA)
write.table(df, paste0(data_dir, "expression_aortic_smooth_muscle_cell.txt"))
    
    

