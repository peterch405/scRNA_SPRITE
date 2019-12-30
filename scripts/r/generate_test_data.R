
library(polyester)
library(Biostrings)
library(Rsubread)
library(ShortRead)
library(readr)

m_cdna <-readDNAStringSet("/mnt/data/genomes/GRCm38.p6/Mus_musculus.GRCm38.cdna.all.fa.gz")
h_cdna <- readDNAStringSet("/mnt/data/genomes/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa.gz")
# subset the FASTA file to first 20 transcripts



#' Removes duplicate sequences from DNAStringSet object.
#'
#' Given a DNAStringSet object, the function dereplicates reads and 
#' adds counts=X to the definition line to indicate replication. 
#'
#' @param dnaSet DNAStringSet object to dereplicate. 
#'
#' @return DNAStringSet object with names describing frequency of repeat.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, 
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples 
#' dnaSet <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC", 
#' "CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC") 
#' dereplicateReads(dnaSet)
dereplicateReads <- function(dnaSet) {
  if(!is(dnaSet,"DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  
  if(is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), sep="-")
  }
  
  dnaSet <- dnaSet[order(dnaSet)]
  counts <- BiocGenerics::table(dnaSet)
  dnaSet <- unique(dnaSet)
  rows <- match(dnaSet, names(counts))
  names(dnaSet) <- paste0(names(dnaSet), 
                          "counts=", as.integer(counts[rows]))
  return(dnaSet)
}





#'Rename and output x number of transcripts
#'
#'
reformat_fa <- function(dnastr, recs_out=1000){
  prot_cod <- dnastr[grepl("transcript_biotype:protein_coding", names(dnastr))]
  #remove identical sequences
  prot_cod <- dereplicateReads(prot_cod)
  new_names <- sapply(strsplit(sapply(strsplit(names(prot_cod), " "), "[", 7), "gene_symbol:"), "[", 2)
  names(prot_cod) <- new_names
  #remove isoforms
  prot_cod_unq <- prot_cod[!duplicated(names(prot_cod))]
  #return x records
  return(prot_cod_unq[1:recs_out])
  
}



mouse_1k <- reformat_fa(m_cdna, 1000)
human_1k <- reformat_fa(h_cdna, 1000)

writeXStringSet(mouse_1k, "/mnt/data/scRNA_SPRITE/mouse_1k_transcripts.fa.gz", compress = TRUE)
writeXStringSet(human_1k, "/mnt/data/scRNA_SPRITE/human_1k_transcripts.fa.gz", compress = TRUE)


fold_changes <- matrix(rep(1,2000), nrow=1000)
head(fold_changes)

#' #'
#' #'
#' #'
#' make_simulated_reads <- function(fasta, out_dir, fold_changes, coverage=20){
#'   # ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
#'   # here all transcripts will have ~equal FPKM
#'   readspertx <- round(coverage * width(mouse_1k) / 100)
#'   
#'   # simulation call:
#'   simulate_experiment(fasta, reads_per_transcript=readspertx, 
#'                       num_reps=c(1,1), fold_changes=fold_changes, outdir=out_dir)   
#' }
#' 
#' make_simulated_reads("D://scRNA_SPRITE/mouse_1k_transcripts.fa", "simulated", fold_changes, 20)


# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta("D://genomes/GRCm38.p6/Mus_musculus.GRCm38.cdna.all.fa.gz")
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)

# Assign a random TPM value to each non-duplicated transcript sequence
TPMs <- rep(0, nrow(transcripts))
TPMs[!transcripts$Duplicate] <- rexp(1000)
TPMs <- rep(2, 1000)

# Generate actual reads.
# The output read file is my-simulated-sample_R1.fastq.gz 
# The true read counts are returned.
true.counts_mouse <- simReads("/mnt/data/scRNA_SPRITE/mouse_1k_transcripts.fa.gz", TPMs, "/mnt/data/scRNA_SPRITE/mouse_1k", paired.end = TRUE)
true.counts_human <- simReads("/mnt/data/scRNA_SPRITE/human_1k_transcripts.fa.gz", TPMs, "/mnt/data/scRNA_SPRITE/human_1k", paired.end = TRUE)


mouse_r1 <- readFastq("/mnt/data/scRNA_SPRITE/mouse_1k_R1.fastq.gz")
human_r1 <- readFastq("/mnt/data/scRNA_SPRITE/human_1k_R1.fastq.gz")

#make barcodes read 2 

# quality, sread, id


#' generate barcode combinations
#' write barcodes into read name so can check accuracy of barcodeID
#' 
barcode_combinations <- function(config_path, num_barcodes){
  config <- read_delim(config_path, 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE, skip = 2)
  #make named vector
  odd_bc <- config$X3[config$X1 == "ODD"]
  names(odd_bc) <- config$X2[config$X1 == "ODD"]
  even_bc <- config$X3[config$X1 == "EVEN"]
  names(even_bc) <- config$X2[config$X1 == "EVEN"]
  y_bc <- config$X3[config$X1 == "Y"]
  names(y_bc) <- config$X2[config$X1 == "Y"]
  
  #all combinations
  all_barcodes <- expand.grid(y_bc, odd_bc, even_bc, odd_bc, stringsAsFactors = FALSE)
  # all_names <- expand.grid(names(odd_bc), names(even_bc), names(odd_bc), names(y_bc), stringsAsFactors = FALSE)
  spacer <- "ATGCAT"
  bc_to_use <- sample(seq(1,nrow(all_barcodes)), num_barcodes)
  barcodes_seqs <- apply(all_barcodes[bc_to_use,], 1, function(x) paste(x, collapse = spacer))
  
  return(barcodes_seqs)
  
}

# sapply(as.character(rnames), function(x) paste(x, )



config_path <- "/mnt/data/scRNA/20191205_scrna_mus_hs_2/config.txt"
barcode_seqs <- barcode_combinations(config_path, 20000)



#'Make perfect barcode R2 fastq
#'
#'@param barcodes_seqs output from barcode_combinations, a vector of barcode sequences
#'@param R1_srq ShortReadQ object with R1 reads
#'
#'@return ShortReadQ object
#'
make_barcode_R2 <- function(barcodes_seqs, R1_srq){
  
  #read in quality score reference
  quality.reference <- system.file("qualf", "ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt", 
                                   package = "Rsubread")
  con <- file(quality.reference, "r")
  qual_ref <- readLines(con)
  close(con)
  
  
  r2_name <- as.character(id(R1_srq))
  r2_seq <- list()
  r2_qual <- list()
  set.seed(1234)
  counter <- 1
  total <- 0
  while(total < length(r2_name)){
    total <- total + 1
    if(total %% 1000 == 0){
      counter <- counter + 1
      
    }
    
    if(total %% 100000 == 0){
      print(total)
    }
    
    r_seq <- unname(barcodes_seqs[counter])
    r2_seq[[total]] <- r_seq
    seq_len <- nchar(r_seq)
    r2_qual[[total]] <- substr(sample(qual_ref, 1), 1, seq_len)
    
  }

  #make DNAStringSet and then ShortReadQ  
  srq <- ShortReadQ(sread = DNAStringSet(as.character(unlist(r2_seq))), 
                   quality = BStringSet(as.character(unlist(r2_qual))), 
                   id = BStringSet(r2_name))
  
  return(srq)
}


mouse_r2 <- make_barcode_R2(barcodes_seqs[1:10000], mouse_r1)
human_r2 <- make_barcode_R2(barcodes_seqs[10000:20000], human_r1)

writeFastq(mouse_r2, "/mnt/data/scRNA_SPRITE/mouse_1k_barcodes_R2.fastq.gz")
writeFastq(human_r2, "/mnt/data/scRNA_SPRITE/human_1k_barcodes_R2.fastq.gz")



# 
# 
# 
# 
# parse_config <- function(config_path){
#   
#   
#   
# }
# 
# library(data.table)
# set.seed(0L)
# 
# output <- data.table(first=odd_bc)[, .(second=even_bc), by=.(first)]
# 

#read configuration 
# con <- file(config_path,"r")
# found <- FALSE
# while(found == FALSE){
#   layout <- readLines(con,n=1)
#   
# }
# close(con)