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
  ## function copied from the genomation package
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
  ## count e.g. sum of coverage within bins (windows)
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
    
  stopifnot(identical(sort(seqlevels(bins)), sort(names(numvar))))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         func(views)
                       })
  new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}


### MAIN

# read lincRNAs
lincRNAs <- gffToGRanges(lincRNA_annotation_path)
lincRNAs <- unlist(lincRNAs)
lincRNAs = lincRNAs[lincRNAs$type=="transcript"]

# read introns and exons
#txdb <- makeTxDbFromGFF(gene_annotation_path, format="gtf", circ_seqs=character())
#trak2_txs <- transcriptsBy(txdb, by="gene")
#introns <- intronsByTranscript(txdb, use.names=TRUE)
#exons <- exonsBy(txdb, use.names=TRUE)
#saveRDS(trak2_txs, paste0(data_dir, "transcripts.rds"))
#saveRDS(introns, paste0(data_dir, "introns.rds"))
#saveRDS(exons, paste0(data_dir, "exons.rds"))
exons <- readRDS(paste0(data_dir, "exons.rds"))
#introns <- readRDS(paste0(data_dir, "introns.rds"))
#introns_transcript_id = names(introns)

#a <- mclapply(introns[1:100],function(ranges){if(length(ranges)) ranges$idx=1:length(ranges);ranges}, mc.cores=30)
#b <- Reduce('c', a)

#introns$intron_rank <- sapply(elementLengths(introns),seq,from=1)
#saveRDS(introns, paste0(data_dir, "introns_ranks1.rds"))
introns <- readRDS(paste0(data_dir, "introns_ranks1.rds"))

#introns <- unlist(introns)
#introns$transcript_id <-introns_name
exons <- unlist(exons)
exons$transcript_id <- names(exons)


#saveRDS(introns, paste0(data_dir, "introns_grl_unlist.rds"))
#saveRDS(exons, paste0(data_dir, "exons_grl_inlist.rds"))
#introns <- readRDS(paste0(data_dir, "introns_grl_unlist.rds"))
#introns <- introns[ which(!is.na(introns$transcript_id)), ] #TODO:
#exons <- readRDS(paste0(data_dir, "exons_grl_inlist.rds"))


#read data from Marcel
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
#derived_from = as.character(metadata$Derived.from)
bw_files = paste0(data_dir, file_acc, ".bigWig")

readBigWig.parallel <- function(i, targets=targets){
    #chr1.size = seqinfo(BigWigFile(fil))["chr1"]
    #chr1.gr <- GRanges(seqnames="chr1", ranges=IRanges(start=1, end=seqlengths(chr1.size)), strand="*")
    #readBigWig(targets[i], windows=chr1.gr)$chr1
    readBigWig(targets[i])
}    
#covs <- mclapply(1:nrow(metadata), 
#                  readBigWig.parallel,
#                  targets=bw_files,
#                  mc.cores=30
#                  )
#names(covs) <- file_acc
#saveRDS(covs, paste0(data_dir, "covs.rds"))
#covs <- readRDS(paste0(data_dir, "covs.rds"))

normalize = function(cov, nmbr){
  # for each chromosome separetly
  #TODO: too slow
  for(i in 1:length(cov)){
    cov[[i]] <- cov[[i]] / nmbr  
  }
  cov
}


## normalize 
aa <- unique(metadata$Derived.from)
for(i in 1:length(aa)){
   print(i)
   samples_both_strands <- metadata [ metadata$Derived.from==aa[i] ,]
   sum_covs <- sum(sapply(covs[samples_both_strands$id], sum))
   covs[[samples_both_strands$id[1]]] <- normalize(covs[[samples_both_strands$id[1]]], sum_covs)
   covs[[samples_both_strands$id[2]]] <- normalize(covs[[samples_both_strands$id[2]]], sum_covs)
}


#saveRDS(covs, paste0(data_dir, "covs_normalized.rds"))
covs <- readRDS(paste0(data_dir, "covs_normalized.rds"))
names(covs) <- file_acc



#old
#normalization. Samples from every strands are normalized by sum of coverage from both strands
#sum_covs <- sum(sapply(covs, sum)) #sum within derived from
#covs <- lapply(1:length(covs), function(x) covs[[x]] / sum_covs)
#names(covs) <- file_acc


# for each replicate
metadata.tmp <- metadata
metadata <- metadata[,c("File.accession", "Output.type","strand",
		      "amount", "id", "Biosample.term.name", "Experiment.accession", "Derived.from")]

		      
binned.arithmfunc <- function(bins, numvar, mcolname, func)
{
  ## count e.g. sum of coverage within bins (windows)
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
		      

names(lincRNAs) = lincRNAs$transcript_id	      
		      
for(i in 1:length(metadata$File.accession)){

  print(i)
  introns <- binned.arithmfunc(introns, covs[[i]], names(covs)[i], viewSums)
  exons <- binned.arithmfunc(exons, covs[[i]], names(covs)[i], viewSums)
  lincRNAs <- binned.arithmfunc(lincRNAs, covs[[i]], names(covs)[i], viewSums)
  
  #normalization
  mcols(introns)[,names(covs)[i]] <- mcols(introns)[,names(covs)[i]]*1. / sum(mcols(introns)[,names(covs)[i]])
  mcols(exons)[,names(covs)[i]] <- mcols(exons)[,names(covs)[i]]*1. / sum(mcols(exons)[,names(covs)[i]])
  mcols(lincRNAs)[,names(covs)[i]] <- mcols(lincRNAs)[,names(covs)[i]]*1. / sum(mcols(lincRNAs)[,names(covs)[i]])

}

#which_covs = which(names(covs[[i]]) %in% seqlevels(lincRNAs))
#covs_lincRNAs <- covs[which_covs]
#lincRNAs <- binned.arithmfunc(lincRNAs, covs[[i]], names(covs)[i], viewSums)


# sum samples on + and - strand
#(introns$transcript_id)
#metadata$Derived.from

introns.tmp = introns
exons.tmp = exons


introns_ave <- data.frame(1:length(introns))
exons_ave <- data.frame(1:length(exons))
lincRNA_ave <- data.frame(1:length(lincRNAs))

aa <- as.character(unique(metadata$Derived.from))
for(i in 1:length(aa)){
   print(i)
   samples_both_strands <- as.character(metadata [ metadata$Derived.from==aa[i] ,]$File.accession)
   #introns
   x1 <- rowSums(as.data.frame(mcols(introns)[samples_both_strands]), na.rm = TRUE)
   y1 <- rowSums(as.data.frame(mcols(exons)[samples_both_strands]), na.rm = TRUE)
   z1 <- rowSums(as.data.frame(mcols(lincRNAs)[samples_both_strands]), na.rm = TRUE)

   introns_ave <- cbind(introns_ave, x1)
   exons_ave <- cbind(exons_ave, y1)
   lincRNA_ave <- cbind(lincRNA_ave, z1)
   #values(introns) <- cbind(values(introns), DataFrame(o=xx))
}

introns$average = rowMeans(introns_ave)
exons$average = rowMeans(exons_ave)
lincRNAs$average = rowMeans(lincRNA_ave)


# take average of replicates within sample
df.introns = data.frame(txid=introns$transcript_id, type="intron", idx=introns$intron_rank, expression=introns$average)
df.exons = data.frame(txid=exons$transcript_id, type="exon", idx=exons$exon_rank, expression=exons$average)
df.lincRNA = data.frame(txid=lincRNAs$transcript_id, type="lincRNA", idx=0, expression=lincRNAs$average)


df <- rbind(df.introns, df.exons, df.lincRNA)
write.table(df, paste0(data_dir, "expression_aortic_smooth_muscle_cell.txt"))
    
    
    
    
