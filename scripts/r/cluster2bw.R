library(GenomicRanges)

################################################################################
#' generate genome-wide bins for counting purpose
#'
#' From: https://github.com/stjude/ChIPseqSpikeInFree
#'
#' Given a chrom.size file, this function allows you to generate a your own sliding windows (bins).
#' @param chromFile chrom.size. Given "hg19","mm10","mm9" or "hg38", will load chrom.size file from package folder.
#' @param binSize size of bins (bp)
#' @param overlap overlaps between two consecutive bins (bp)
#' @param withChr chromosome names in bin File have chr if set withChr to TRUE; FALSE - no chr
#' @return A data.frame of generated bins
#' @export
#' @examples
#' ## 1. generate a mm10 binFile without chr and use a binSize of 1000 bp
#' ## and overlap of 500 bp between two consecutive bins
#' 
#' ## "mm10" will be parsed as system.file("extdata", "mm10.chrom.sizes",
#' ##       package = "ChIPseqSpikeInFree")
#' # binDF <- GenerateBins(chromFile="mm10",binSize=1000, overlap=500,
#' # withChr=FALSE,prefix="mm10")
#' 
#' ## 2. generate a hg19 binFile with chr and use a binSize of 2000 bp
#' 
#' ## "hg19" will be parsed as system.file("extdata", "hg19.chrom.sizes",
#' ##       package = "ChIPseqSpikeInFree")
#' # binDF<- GenerateBins(chromFile="hg19",binSize=2000, overlap=0,
#' # withChr=TRUE,prefix="hg19")
GenerateBins <- function(chromFile, binSize = 1000, overlap = 0, withChr = TRUE) {
  # generate genome-wide bins for counting purpose
  
  options(stringsAsFactors = FALSE)
  if (!file.exists(chromFile)) {
    chromFile <- system.file("extdata", paste0(chromFile, ".chrom.sizes"), package = "ChIPseqSpikeInFree")
    if (!file.exists(chromFile)) {
      stop(paste0("\nchromFile was not found:", chromFile, "\n"))
    }
  }
  if (binSize < 200 && binSize > 10000) {
    stop(paste0("\n**Recommended binSize range 200 ~ 10000 (bp); Your binSize is", binSize, " **\n"))
  }
  cat(paste0("\n\tFollowing file will be used to generate sliding windows: ", basename(chromFile)))
  chromSize <- read.table(chromFile, sep = "\t", header = F, fill = TRUE, quote = "")
  options(scipen = 999) # disable scientific notation when writing out
  counter <- NULL
  excluded <- grepl("random", chromSize[, 1]) | grepl("Un", chromSize[, 1]) | grepl("alt", chromSize[, 1]) | grepl("hap", chromSize[, 1])
  chromSize <- chromSize[!excluded, ]
  chrNotation <- gsub("chrM", "MT", chromSize[, 1])
  chrNotation <- gsub("chr", "", chrNotation)
  if (withChr) {
    chrNotation <- gsub("^", "chr", chrNotation)
    chrNotation <- gsub("chrMT", "chrM", chrNotation)
  }
  chromSize[, 1] <- chrNotation
  beds <- NULL
  for (r in 1:nrow(chromSize)) {
    chrName <- chromSize[r, 1]
    chrLength <- chromSize[r, 2]
    realBinSize <- binSize - overlap
    starts <- seq(0, chrLength, realBinSize) + 1
    ends <- c(seq(binSize, chrLength, realBinSize), chrLength)
    starts <- starts[1:length(ends)]
    bed <- data.frame(chr = chrName, start = starts, end = ends)
    counter <- c(counter, nrow(bed))
    if (r == 1) { # overwrite lines
      beds <- bed
    } else {
      beds <- rbind(beds, bed)
    }
  }
  return(beds)
}


################################################################################
#' Cluster file parser
#' 
#' @param clusters_path File path to clusters file
#' 
#' @return Dataframe with first column being the barcode, second clusters
#' 
#' 
#' 
parse_clusters <- function(path){
  
  clusters <- readr::read_csv(path, col_names = FALSE)
  
  counts <-lapply(clusters$X1, function(x){
    row <-unlist(strsplit(x, "\t"))
    count <- length(row[-1])
    return(c(row[1], paste(row[-1], collapse = "\t"), count))
  })
  
  out <- plyr::ldply(counts, rbind)
  names(out) <- c("Barcode","Clusters", "cluster_size")
  return(out)
}



################################################################################
#'Extract annotation field from cluster read record
#'
#'@param read DNA[mm10]_chr17:39847663-39847713
#'
#'@return mm10
#'
extract_anno <- function(read){
  anno <- unlist(gsubfn::strapply(read, "\\[([^]]*)\\]", back = -1))
  return(anno)
}





#'Using annotation field add counts and percentage of mm10 and hg38 reads within each cluster
#'
#'@param clusters data.frame from parse_clusters function
#'
add_species_counts <- function(clusters){
  
  all_anno <- sapply(as.character(clusters$Clusters), function(x) table(sapply(strsplit(extract_anno(x), ";"), "[", 1)))
  
  all_counts <- plyr::ldply(all_anno, rbind)
  species_counts <- all_counts[,c("mm10", "hg38")]
  species_counts[is.na(species_counts)] <- 0
  
  clusters$mm10 <- species_counts$mm10
  clusters$hg38 <- species_counts$hg38
  
  total <- clusters$mm10 + clusters$hg38

  
  clusters$mm10_perc <- (clusters$mm10/total)*100
  clusters$hg38_perc <- (clusters$hg38/total)*100
  
  return(clusters)
  
}

#'If cells where barcoded by species in first round, use first round barcode to assign identity
#'
#'@param clusters data.frame from parse_clusters function
#'@param mouse_barcodes barcode names corresponding to mouse cells
#'@param human_barcodes barcode names corresponding to human cells
#'@param cell_id_barcode position of cell identity barcode
#'
assign_origin <- function(clusters, mouse_barcodes, human_barcodes, cell_id_barcode=4){
  
  clusters$first_odd <- sapply(strsplit(as.character(clusters$Barcode), "\\."), function(x) unlist(x)[cell_id_barcode])
  
  clusters_mouse <- clusters[clusters$first_odd %in% mouse_barcodes,]
  clusters_human <- clusters[clusters$first_odd %in% human_barcodes,]
  
  clusters$origin <- NA
  clusters$origin[clusters$first_odd %in% mouse] <- "Mouse"
  clusters$origin[clusters$first_odd %in% human] <- "Human"
  
  return(clusters)
}


#' Parse DNA[hg38]_chr4:66020898-66020948
#'
#'@param x a cluster read records: DNA[mm10]_chr17:39847686-39847737
#'
#'@return a GenomicRanges object
parse_coord <- function(x){
  split_coords <- lapply(stringr::str_split(x, "(?<!\\[)\\-|[:\\]\\[]|_"), function(x) x[x != ""])
  out <- plyr::ldply(split_coords, rbind)
  names(out) <- c("molecule", "annotation", "chrom", "start", "end")
  out_gr <- GenomicRanges::makeGRangesFromDataFrame(out, keep.extra.columns = TRUE)
  return(out_gr)
}


#'Convert clusters into bedgraph format
#'
cluster2bw <- function(clusters, resolution, mix_filter, mouse_out, human_out){
  
  clusters_filt <- clusters[abs(clusters$mm10_perc-clusters$hg38_perc) < mix_filter,] 
  print("Clusters being used:", nrow(clusters_filt))
  #remove NA origin
  clusters_filt <- clusters_filt[!is.na(clusters_filt$origin),]
  clusters_filt_mouse <- clusters_filt[clusters_filt$origin == "Mouse",]
  clusters_filt_human <- clusters_filt[clusters_filt$origin == "Human",]
  
  #if original idenity mouse, then only output human mapping reads and vice versa
  
  mouse_gr <- lapply(as.character(clusters_filt_mouse$Clusters), function(x) parse_coord(unlist(strsplit(x, "\t"))))
  human_gr <- lapply(as.character(clusters_filt_human$Clusters), function(x) parse_coord(unlist(strsplit(x, "\t"))))
  #combine all ranges                                      
  mouse_gr <- do.call("c", mouse_gr)
  human_gr <- do.call("c", human_gr)
  
  mouse_export <- mouse_gr[grepl("hg38", mouse_gr$annotation)]
  human_export <- human_gr[grepl("mm10", human_gr$annotation)]
  #Sort seqlevels
  mouse_export <- sort(GenomeInfoDb::sortSeqlevels(mouse_export))
  human_export <- sort(GenomeInfoDb::sortSeqlevels(human_export))
  
  human_seq_levels <- parse_chromsizes("/mnt/data/genomes/GRCh38/hg38.chrom.sizes", "hg38")
  seqinfo(mouse_export, new2old=NULL) <- human_seq_levels
  
  mouse_seq_levels <- parse_chromsizes("/mnt/data/genomes/GRCm38.p6/mm10.chrom.sizes", "mm10")
  human_export <- dropSeqlevels(human_export, c("chr20", "chr21", "chr22"))
  seqinfo(human_export, new2old=NULL) <- mouse_seq_levels
  
  
  human_bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(mouse_export),
                                          tilewidth=resolution,
                                          cut.last.tile.in.chrom=TRUE)
  
  mouse_bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(human_export),
                                          tilewidth=resolution,
                                          cut.last.tile.in.chrom=TRUE)
  
  m_gr_cov <- coverage(mouse_export)
  human_bins_cov <- GenomicRanges::binnedAverage(human_bins, m_gr_cov, "score")
  
  h_gr_cov <- coverage(human_export)
  mouse_bins_cov <- GenomicRanges::binnedAverage(mouse_bins, h_gr_cov, "score")
  
  
  rtracklayer::export.bw(human_bins_cov, mouse_out)
  rtracklayer::export.bw(mouse_bins_cov, human_out)
  
}

cluster2bw(clusters, 1000, 95, "/mnt/data/scRNA/20191205_scrna_mus_hs_2/mouse_human_contaminants.bw", 
           "/mnt/data/scRNA/20191205_scrna_mus_hs_2/human_mouse_contaminants.bw")

#' Make seqlengths named vector
#' 
#' 
parse_chromsizes <- function(chrom_sizes, genome, withChr=TRUE, mouse=FALSE){
  chrom_size <- read.table(chrom_sizes, sep = "\t", header = F, fill = TRUE, quote = "")
  excluded <- grepl("random", chrom_size[, 1]) | grepl("Un", chrom_size[, 1]) | grepl("alt", chrom_size[, 1]) | grepl("hap", chrom_size[, 1])
  chrom_size <- chrom_size[!excluded, ]
  
  names(chrom_size) <- c("seqnames", "seqlengths")
  
  chrom_size$isCircular <- NA
  chrom_size$genome <- genome
  
  seqinfo_out <- Seqinfo(chrom_size$seqnames, chrom_size$seqlengths, chrom_size$isCircular, chrom_size$genome)
  seqinfo_out_sorted <- GenomeInfoDb::sortSeqlevels(seqinfo_out)
  
  return(seqinfo_out_sorted)
}





clusters <- parse_clusters("/mnt/data/scRNA/20191205_scrna_mus_hs_2/workup/clusters/scrna-mus-hs-2_S1_L001.clusters")
clusters <- add_species_counts(clusters)

mouse <- c("Odd2Bo5", "Odd2Bo13", "Odd2Bo21", "Odd2Bo29", "Odd2Bo37", "Odd2Bo45", "Odd2Bo53", "Odd2Bo61", "Odd2Bo69", "Odd2Bo77", "Odd2Bo85", "Odd2Bo93")
human <- c("Odd2Bo6", "Odd2Bo14", "Odd2Bo22", "Odd2Bo30", "Odd2Bo38", "Odd2Bo46", "Odd2Bo54", "Odd2Bo62", "Odd2Bo70", "Odd2Bo78", "Odd2Bo86", "Odd2Bo94")

clusters <- assign_origin(clusters, mouse, human, 4)

cluster2bw

