

library(rtracklayer)


hg38_gtf <- import.gff("/mnt/data/genomes/GRCh38/Homo_sapiens.GRCh38.97.gtf.gz")

hg38_pc <- hg38_gtf[mcols(hg38_gtf)$type == "gene" & mcols(hg38_gtf)$gene_biotype == "protein_coding"]
mcols(hg38_pc)$score <- 0
export.bed(hg38_pc, "/mnt/data/genomes/GRCh38/GRCh38.97_protein_coding.bed")


hg38_exon <- hg38_gtf[mcols(hg38_gtf)$type == "exon" & mcols(hg38_gtf)$gene_biotype == "protein_coding"]
mcols(hg38_exon)$score <- 0
export.bed(hg38_exon, "/mnt/data/genomes/GRCh38/GRCh38.97_protein_coding_exon.bed")


mm10_gtf <- import.gff("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz")

mm10_pc <- mm10_gtf[mcols(mm10_gtf)$type == "gene" & mcols(mm10_gtf)$gene_biotype == "protein_coding"]
mcols(mm10_pc)$score <- 0
export.bed(mm10_pc, "/mnt/data/genomes/GRCm38.p6/GRCm38.95_protein_coding.bed")

mm10_exon <- mm10_gtf[mcols(mm10_gtf)$type == "exon" & mcols(mm10_gtf)$gene_biotype == "protein_coding"]
mcols(mm10_exon)$score <- 0
export.bed(mm10_exon, "/mnt/data/genomes/GRCm38.p6/GRCm38.95_protein_coding_exon.bed")
