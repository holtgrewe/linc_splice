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
  stopifnot(identical(seqlevels(bins), names(numvar)))
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

# read introns and exons
txdb <- makeTxDbFromGFF(gene_annotation_path, format="gtf", circ_seqs=character())
trak2_txs <- transcriptsBy(txdb, by="gene")
#introns <- intronsByTranscript(txdb, use.names=TRUE)
#exons <- exonsBy(txdb, use.names=TRUE)
saveRDS(trak2_txs, paste0(data_dir, "transcripts.rds"))
#saveRDS(introns, paste0(data_dir, "introns.rds"))
#saveRDS(exons, paste0(data_dir, "exons.rds"))
exons <- readRDS(paste0(data_dir, "exons.rds"))
introns <- readRDS(paste0(data_dir, "introns.rds"))


values(introns[[i]]) <- cbind(values(introns[[i]]), DataFrame( 1:length(introns[[i]]))) 

for(i in 1:length(introns)){

    introns[[i]]$intron_rank <- 1:length(introns[[i]])
}

# read lincRNAs
lincRNAs <- gffToGRanges(lincRNA_annotation_path)

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

# only first replicate
#first_replicate = as.character(metadata$Derived.from)[1]
#metadata = metadata[which(metadata$Derived.from==first_replicate),]



### calculate sum of coverage within introns, exons, lncRNAs

file_acc = as.character(metadata$File.accession)
exp_acc = as.character(metadata$Experiment.accession)
#derived_from = as.character(metadata$Derived.from)
bw_files = paste0(data_dir, file_acc, ".bigWig")


#only for 1st chromosome
lincRNAs = lincRNAs[which(seqnames(lincRNAs)=="chr1"),]
introns = introns[which(seqnames(introns)=="chr1"),]
exons = exons[which(seqnames(exons)=="chr1"),]

chr1.size = seqinfo(BigWigFile(fil))["chr1"]
chr1.gr <- GRanges(seqnames="chr1", ranges=IRanges(start=1, end=seqlengths(chr1.size)), strand="*")

readBigWig.parallel <- function(i, targets=targets, windows=windows){
    readBigWig(targets[i], windows=windows)$chr1
}    
covs <- mclapply(1:nrow(metadata), 
                  readBigWig.parallel,
                  targets=bw_files,
                  windows=chr1.gr,
                  mc.cores=10
                  )
names(covs) <- file_acc

#normalization. Samples from every strands are normalized by sum of coverage from both strands
sum_covs <- sum(sapply(covs, sum)) #sum within derived from
covs <- lapply(1:length(covs), function(x) covs[[x]] / sum_covs)
names(covs) <- file_acc


introns <- unlist(introns)
introns$transcript_id <- names(introns)
exons <- unlist(exons)
exons$transcript_id <- names(exons)
lincRNAs <- unlist(introns)


  i=1
  introns <- binned.arithmfunc(introns, covs[[i]], names(covs)[i], viewSums)
  exons <- binned.arithmfunc(exons, covs[[i]], names(covs)[i], viewSums)
  lincRNAs <- binned.arithmfunc(lincRNAs, covs[[i]], names(covs)[i], viewSums)






    
    
    
    
