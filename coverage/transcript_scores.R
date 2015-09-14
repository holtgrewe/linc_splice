
library(rtracklayer)
library(genomation)
library(GenomicFeatures)

### PATHS
data_dir <- "/home/kasia/lnc_rna_splice_site_hackathon/linc_splice/coverage/data/"
#data_dir = "/home/kwreczy/linc_splice/coverage/data/"
#bw_paths = list.files(path = data_dir, pattern = "*.bigWig")
gene_annotation_path = "~/gencode.v19.annotation.gtf.gz"
lincRNA_annotation_path = paste0(data_dir, "gencode.v19.long_noncoding_RNAs.gtf")
metadata_path = "/home/kasia/lnc_rna_splice_site_hackathon/linc_splice/coverage/data/metadata.tsv"

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

count_coverage <- function(bw_paths, windows=NULL){
    covs = list()    
    for(i in 1:length(bw_paths)){
      print(i)
      if(is.null(windows)){
	
	cov_sample = readBigWig(target=bw_paths[i])
	
      }else{
      
	cov_sample = try(readBigWig(target=bw_paths[i], windows=windows))
	if(is(aa,"try-error")) {
	  windowswithoutchrY=windows[!as.character(seqnames(windows))=="chrY",]
	  seqlevels(windowswithoutchrY) <- c(sapply(1:22, function(x) paste0("chr", x)), "chrX")

	  cov_sample = readBigWig(target=bw_paths[i], windows=windows)
	}else{
	  cov_sample = readBigWig(target=bw_paths[i], windows=windows)
	}
	
      }
      covs[[i]] <- cov_sample
    }
   covs
}


### MAIN

# read introns and exons
txdb <- makeTxDbFromGFF(gene_annotation_path, format="gtf", circ_seqs=character())
trak2_txs <- transcriptsBy(txdb, by="gene")
introns <- intronsByTranscript(txdb, use.names=TRUE)
exons <- exonsBy(txdb, use.names=TRUE)

#saveRDS(introns, "introns.rds")
#saveRDS(exons, "exons.rds")
introns <- readRDS("introns.rds")
exons <- readRDS("exons.rds")

# read lincRNAs
lincRNAs <- gffToGRanges(lincRNA_annotation_path)

metadata <- read.table(metadata_path, 
			sep="\t",
			header=TRUE)	
metadata$strand = "1"
metadata[which(metadata$Output.type=="plus strand signal of unique reads"),]$strand = "+"
metadata[which(metadata$Output.type=="minus strand signal of unique reads"),]$strand = "-"
metadata[which(metadata$Output.type=="plus strand signal of all reads"),]$strand = "+"
metadata[which(metadata$Output.type=="minus strand signal of all reads"),]$strand = "-"
metadata$amount = "1"
metadata[which(metadata$Output.type=="plus strand signal of unique reads"),]$amount = "unique"
metadata[which(metadata$Output.type=="minus strand signal of unique reads"),]$amount = "unique"
metadata[which(metadata$Output.type=="plus strand signal of all reads"),]$amount = "all"
metadata[which(metadata$Output.type=="minus strand signal of all reads"),]$amount = "all"


#only uniquely mapped reads
metadata = metadata[which(metadata$amount=="unique"),]
# only aortic smooth muscle cell
metadata = metadata[which(metadata$Biosample.term.name=="aortic smooth muscle cell"),]
# only 1st chromosome
lincRNAs = lincRNAs[which(seqnames(lincRNAs)=="chr1"),]
introns = introns[which(seqnames(introns)=="chr1"),]
exons = exons[which(seqnames(exons)=="chr1"),]


### calculate sum of coverage within introns, exons, lncRNAs

filenames = as.character(metadata$File.accession)
exp_ids <- unique(as.character(metadata$Experiment.accession))

sapply(filenames, readBigWig(target=paste0(data_dir, filenames[i]), windows=ranges(lincRNAs)))

covs <- count_coverage(paste0(data_dir, filenames[2], ".bigWig"), windows=lincRNAs)

introns <- binned.arithmfunc(introns, cov1, "covmax1", viewSums)
exons <- binned.arithmfunc(introns, cov1, "covmax1", viewSums)
lincRNAs <- binned.arithmfunc(introns, cov1, "covmax1", viewSums)


#normalize = function()



    
    
    
    
