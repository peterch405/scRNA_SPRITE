'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process scRNA-seq with the SPRITE barcoding scheme
'''


import os 
import sys
import datetime

#Location of scripts
barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
split_fq = "scripts/python/get_full_barcodes.py"

#Load config.yaml file

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path


#Copy config file into logs
# v = datetime.datetime.now()
# run_date = v.strftime('%Y.%m.%d.')


try:
    email = config['email']
except:
    print("Won't send email on error")
    email = None

try:
    out_dir = config['output_dir']
    print('All data will be written to:', out_dir)
except:
    out_dir = ''
    print('Defaulting to working directory as output directory')

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = out_dir + 'workup/config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "5"
    print('Config "num_tags" not specified, using:', num_tags)

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly'].split('&')
    assert all(i in ['mm10', 'hg38', 'halo'] for i in assembly), 'Only "mm10" or "hg38" currently supported'
    print('Using', ' '.join(assembly))
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'



if config['aligner'] == 'hisat2':
    try:
        anno_gtf = config['anno_gtf']
    except:
        print('Annotation or mask path not specified in config.yaml')
        sys.exit() #no default, exit

    try:
        hisat2_index = config['hisat2_index']
        hisat2_ss = config['hisat2_splice_sites']
    except:
        print('Hisat2 index and splice sites not specified in config.yaml')
        sys.exit()
elif config['aligner'] == 'bowtie2':
    try:
        bowtie2_index = config['bowtie2_index']
    except:
        print('Not bowtie2 index found')
        sys.exit()
elif config['aligner'] == 'star':
    try:
        star_index = config['star_index']
    except:
        print('STAR indexes path not specified in config.yaml')
        sys.exit() #no default, exit


try:
    align_mode = config['align_mode']
    print('Aligning in mode:', align_mode)
except:
    print('Align mode not specified')
    sys.exit()

try:
    samples = config['samples']
    print('Using samples file:', samples)
except:
    samples = './samples.json'
    print('Defaulting to working directory for samples json file')

try:
    trim5_len = config['trim5']
except:
    trim5_len = 0

print("Trimming 5' end:", trim5_len)


#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

# CONFIG = [out_dir + "workup/logs/" + run_date + "config.log"]

#Shared
TRIM = expand(out_dir + "workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
                  sample = ALL_SAMPLES, read = ["R1", "R2"])

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]
MASKED = expand(out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]
BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.barcoded.fastq.gz", sample = ALL_SAMPLES, 
                   read = ["R1", "R2"])
BARCODE_FULL = expand([out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
                       out_dir + "workup/fastqs/{sample}_R2.barcoded_full.fastq.gz"], sample=ALL_SAMPLES)

CUTADAPT = expand([out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz",
                   out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz"],
                   sample = ALL_SAMPLES)

FASTQC = expand([out_dir + "workup/trimmed/{sample}_R1_bfull_trim_fastqc.html",
                 out_dir + "workup/trimmed/{sample}_R2_bfull_trim_fastqc.html"],
                 sample = ALL_SAMPLES)

CLUSTERS = expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
#Hisat2 alignment
Ht2_RNA_ALIGN_SE = expand(out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam", 
                        sample=ALL_SAMPLES, genome=assembly)
Ht2_RNA_ALIGN_PE = expand(out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam", 
                           sample=ALL_SAMPLES, genome=assembly)

Ht2_ANNO_RNA_SE = expand(out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam.featureCounts.bam",
                  sample=ALL_SAMPLES, genome=assembly)
Ht2_ANNO_RNA_PE = expand(out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam.featureCounts.bam",
                  sample=ALL_SAMPLES, genome=assembly)
CLUSTERS_PLOT = [out_dir + "workup/clusters/cluster_sizes.pdf", out_dir + "workup/clusters/cluster_sizes.png"]

BOWTIE2_ALIGN = expand(out_dir + "workup/alignments/{sample}.DNA.b2.mapq20.bam", sample=ALL_SAMPLES)

STAR_RNA = expand(out_dir + "workup/alignments/{sample}.{genome}.Aligned.sortedByCoord.out.mapq20.bam", 
                   sample=ALL_SAMPLES, genome=assembly)

if config['aligner'] == 'hisat2':
    if align_mode == "SE":
        rule all:
            input: CONFIG + ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
                Ht2_RNA_ALIGN_SE + Ht2_ANNO_RNA_SE + 
                CLUSTERS + MULTI_QC + CLUSTERS_PLOT
    elif align_mode == "PE":
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
                CUTADAPT + FASTQC + Ht2_RNA_ALIGN_PE + Ht2_ANNO_RNA_PE +
                CLUSTERS + MULTI_QC + CLUSTERS_PLOT
elif config['aligner'] == 'bowtie2':
        rule all:
            input: CONFIG + ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
                CLUSTERS + MULTI_QC + CLUSTERS_PLOT + BOWTIE2_ALIGN
elif config['aligner'] == 'star':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
                CUTADAPT + FASTQC + STAR_RNA #+
                #CLUSTERS + MULTI_QC + CLUSTERS_PLOT

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')



RNA_star_params = "--runMode alignReads \
--outFilterMultimapNmax 50 \
--outFilterScoreMinOverLread 0.30 \
--outFilterMatchNminOverLread 0.30 \
--outFilterIntronMotifs None \
--alignIntronMax 50000 \
--alignMatesGapMax 1000 \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--alignIntronMin 80 \
--alignSJDBoverhangMin 5 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 250000000 \
--outSAMattributes All \
--readFilesCommand zcat \
--sjdbOverhang 100 \
--twopassMode Basic"

####################################################################################################
#Trimming and barcode identification
####################################################################################################

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
        out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
        out_dir + "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
        out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz",
        out_dir + "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        out_dir + "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}"


# rule log_config:
#     '''Copy config.yaml and place in logs folder with the date run
#     '''
#     input:
#         config_path
#     output:
#         out_dir + "workup/logs/" + run_date + "config.log"
#     shell:
#         "cp {input} {output}"




#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz"
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp(out_dir + "workup/{sample}.ligation_efficiency.txt")
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"
      

rule full_barcode:
    '''
    remove incomplete barcodes
    '''
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2 = out_dir + "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_short.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R2.barcoded_full.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R2.barcoded_short.fastq.gz"
    log:
        r1 = out_dir + "workup/logs/{sample}_DPM_R1.log",
        r2 = out_dir + "workup/logs/{sample}_DPM_R1.log"
    shell:
        '''
        python {split_fq} --r1 {input.r1} &> {log.r1}
        python {split_fq} --r1 {input.r2} &> {log.r2}
        '''


rule cutadapt_pe:
    '''
    Trim barcode sequences for paired end alignment
    --minimum-length 20 (don't filter on read length as single end alignment still valuable)
    -n 3
    '''
    input:
        [out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz", 
        out_dir + "workup/fastqs/{sample}_R2.barcoded_full.fastq.gz"]
    output:
        fastq1=out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz",
        fastq2=out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz",
        qc=out_dir + "workup/trimmed/{sample}.trim.qc.txt"
    threads: 10
    params:
        adapters_r1 = "-a CAAGTCA",
        adapters_r2 = "-G TGACTTGNTTCG",
        others = "-e 0.15"
    log:
        out_dir + "workup/logs/{sample}.cutadapt.log"
    wrapper:
        "0.38.0/bio/cutadapt/pe"


rule fastqc:
    input:
        [out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz",
        out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz"]
    output:
        [out_dir + "workup/trimmed/{sample}_R1_bfull_trim_fastqc.html",
        out_dir + "workup/trimmed/{sample}_R2_bfull_trim_fastqc.html"]
    conda:
        "envs/trim_galore.yaml"
    log:
        out_dir + "workup/logs/{sample}.fastqc.log"
    shell:
        "fastqc {input} &> {log}"


############################################################################################
#RNA alignment
############################################################################################



rule hisat2_align:
    '''
    #from clusterflow pipeline
    # we are currently using a very high penalty score for soft-clipping (--sp 1000,1000)
    #because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
    # we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
    # we are also filtering non-primary alignments (-F 256)
    #filter on mapq score of 20 (Skip alignments with MAPQ smaller than 20)
    #-U FILE Write alignments that are not selected by the various filter options to FILE
    '''
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        mapped=out_dir + "workup/alignments/{sample}.SE.hisat2.mapq20.bam"
    threads: 10
    conda:
        "envs/hisat2.yaml"
    log:
        "workup/logs/{sample}.hisat2.log"
    shell:
        '''
        (hisat2 --sp 1000,1000 \
        -p 10 \
        -t \
        --phred33 \
        --known-splicesite-infile {hisat2_ss} \
        -x {hisat2_index} \
        -U {input.fq} | \
        samtools view -bq 20 -F 4 -F 256 - > {output.mapped}) &> {log}
        '''

 # params:
        # index = lambda wildcards, output, input: if 'hg38'
  
rule hisat2_align_pe:
    '''
    Notes:
        from clusterflow pipeline
        we are currently using a very high penalty score for soft-clipping (--sp 1000,1000)
        because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
        we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
        we are also filtering non-primary alignments (-F 256)
        filter on mapq score of 20 (Skip alignments with MAPQ smaller than 20)
        -U FILE Write alignments that are not selected by the various filter options to FILE
        --dta/--downstream-transcriptome-assembly
            Report alignments tailored for transcript assemblers including StringTie. With this option, 
            HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to 
            fewer alignments with short-anchors, which helps transcript assemblers improve significantly 
            in computation and memory usage.
    '''
    input:
        fq_1=out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
        fq_2=out_dir + "workup/fastqs/{sample}_R2.barcoded_full.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam"
    threads: 10
    params:
        ss = lambda wildcards: hisat2_ss[wildcards.genome],
        index = lambda wildcards: hisat2_index[wildcards.genome]
    conda:
        "envs/hisat2.yaml"
    log:
        out_dir + "workup/logs/{sample}.{genome}.hisat2.log"
    shell:
        ''' 
        (hisat2 \
        --sp 1000,1000 \
        --no-discordant \
        --dta \
        -p {threads} \
        -t \
        --phred33 \
        --trim5 {trim5_len} \
        --known-splicesite-infile {params.ss} \
        -x {params.index} \
        -1 {input.fq_1} -2 {input.fq_2} | \
        samtools view -bq 20 -F 4 -F 256 - > {output}) &> {log}
        '''
# 
# 
# --no-mixed \



rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.DNA.b2.mapq20.bam"
    threads: 10
    params:
        bowtie2_idx = lambda wildcards: bowtie2_index[wildcards.genome]
    log:
        out_dir + "workup/logs/{sample}.bowtie2.log"
    conda:
        "envs/bowtie2.yaml"
    shell:
        "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        --trim5 {trim5_len} \
        -x {params.bowtie2_idx} \
        -U {input.fq} | \
        samtools view -bq 20 -F 4 -F 256 - > {output}) &> {log}"





rule star_align_rna:
    '''
    Align RNA with STAR to the genome first, annotate repeats, 
    anything that does not align, realign with bowtie2 to our repeats reference
    '''
    input:
        fq = out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.{genome}.Aligned.sortedByCoord.out.mapq20.bam",
        out_dir + "workup/alignments/{sample}.{genome}.Aligned.sortedByCoord.out.bam",
        out_dir + "workup/alignments/{sample}.{genome}.Log.final.out",
        out_dir + "workup/alignments/{sample}.{genome}.Log.out",
        out_dir + "workup/alignments/{sample}.{genome}.Log.progress.out",
        out_dir + "workup/alignments/{sample}.{genome}.SJ.out.tab",
        out_dir + "workup/alignments/{sample}.{genome}.unmapped.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{genome}.star.log"
    threads: 10
    params:
        star_idx = lambda wildcards: star_index[wildcards.genome]
    conda:
        'envs/star.yaml'
    shell:
        '''
        STAR {RNA_star_params} \
        --runThreadN {threads} \
        --genomeDir {params.star_idx} \
        --readFilesIn {input.fq} \
        --outFileNamePrefix {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}. &> {log}

        mv {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.Unmapped.out.mate1 \
            {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.unmapped.fastq

        pigz {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.unmapped.fastq

        samtools view -bq 20 {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.Aligned.sortedByCoord.out.bam > \
            {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.Aligned.sortedByCoord.out.mapq20.bam
        '''




rule annotate_rna:
    '''
    -M                  Multi-mapping reads will also be counted. For a multi-
                        mapping read, all its reported alignments will be 
                        counted. The "NH" tag in BAM/SAM input is used to detect 
                        multi-mapping reads.
    -s <int or string>  Perform strand-specific read counting. A single integer
                        value (applied to all input files) or a string of comma-
                        separated values (applied to each corresponding input
                        file) should be provided. Possible values include:
                        0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                        Default value is 0 (ie. unstranded read counting carried
                        out for all input files).
    -t <string>         Specify feature type in GTF annotation. "exon" by 
                        default. Features used for read counting will be 
                        extracted from annotation using the provided value.
    -g <string>         Specify attribute type in GTF annotation. "gene_id" by 
                        default. Meta-features used for read counting will be 
                        extracted from annotation using the provided value.
    '''
    input:
        out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam" if align_mode == "SE" else
        out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam"
    threads: 10
    output:
        bam=out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam.featureCounts.bam" if align_mode == "SE" else
            out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam.featureCounts.bam",
        counts=out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam.featureCounts.txt" if align_mode == "SE" else
               out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam.featureCounts.txt"
    log:
        out_dir + "workup/logs/{sample}.{genome}.anno.log"
    params:
        anno = lambda wildcards: anno_gtf[wildcards.genome]
    conda:
        "envs/annotate_rna.yaml"
    shell:
        '''
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 1 \
        -g gene_name -a {params.anno} -o {output.counts} \
        {input}
        '''


rule add_chr:
    input:
        out_dir + "workup/alignments/{sample}.SE.h2.{genome}.mapq20.bam.featureCounts.bam" if align_mode == "SE" else
        out_dir + "workup/alignments/{sample}.PE.h2.{genome}.mapq20.bam.featureCounts.bam"
    output:
        out_dir + "workup/alignments/{sample}.{genome}.chr.bam"
    log:
        out_dir + "workup/logs/{sample}.{genome}.chr.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {wildcards.genome} &> {log}
        '''


rule make_clusters:
    input:
        expand(out_dir + "workup/alignments/{sample}.{genome}.chr.bam", 
               sample=ALL_SAMPLES, genome=assembly)
    output:
        out_dir + "workup/clusters/{sample}.clusters"
    log:
        out_dir + "workup/clusters/{sample}.make_clusters.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} \
        -i {input} \
        -o {output} \
        -n {num_tags} &> {log}"



rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc workup -o {out_dir}workup/qc"


rule plot_cluster_size:
    input:
        expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/clusters/cluster_sizes.pdf",
        out_dir + "workup/clusters/cluster_sizes.png"
    log:
        out_dir + "workup/logs/cluster_sizes.log"
    conda:
        "envs/r.yaml"
    shell:
        '''
        Rscript scripts/r/get_cluster_size_distribution.r \
            {out_dir}workup/clusters/ \
            .clusters
        '''


# rule gene_body_coverage:
#     input:

#     output:
        
#     conda:
#         "envs/qc.yaml"

# geneBody_coverage.py -r hg19.housekeeping.bed -i test.bam  -o output