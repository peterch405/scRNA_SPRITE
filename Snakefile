'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process scRNA-seq with the SPRITE barcoding scheme
'''


import os 
import sys

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
    bid_config = 'workup/config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "5"
    print('Config "num_tags" not specified, using:', num_tags)

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'


try:
    anno_gtf = config['anno_gtf'][config['assembly']]
except:
    print('Annotation or mask path not specified in config.yaml')
    sys.exit() #no default, exit

try:
    hisat2_index = config['hisat2_index'][config['assembly']]
    hisat2_ss = config['hisat2_splice_sites'][config['assembly']]
except:
    print('Hisat2 index and splice sites not specified in config.yaml')
    sys.exit()



try:
    samples = config['samples']
    print('Using samples file:', samples)
except:
    samples = './samples.json'
    print('Defaulting to working directory for samples json file')

#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])


#Shared
TRIM = expand("workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand("workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
                  sample = ALL_SAMPLES, read = ["R1", "R2"])

LE_LOG_ALL = ["workup/ligation_efficiency.txt"]
MASKED = expand("workup/alignments/{sample}.DNA.chr.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = ["workup/qc/multiqc_report.html"]
BARCODEID = expand("workup/fastqs/{sample}_{read}.barcoded.fastq.gz", sample = ALL_SAMPLES, 
                   read = ["R1", "R2"])
BARCODE_FULL = expand(out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz", sample=ALL_SAMPLES)

CLUSTERS = expand("workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
#Hisat2 alignment
Ht2_RNA_ALIGN = expand("workup/alignments/{sample}.RNA.hisat2.mapq20.bam", 
                        sample=ALL_SAMPLES)
Ht2_ANNO_RNA = expand("workup/alignments/{sample}.RNA.hisat2.mapq20.bam.featureCounts.bam",
                  sample=ALL_SAMPLES)
CLUSTERS_PLOT = [out_dir + "workup/clusters/cluster_sizes.pdf", out_dir + "workup/clusters/cluster_sizes.png"]

rule all:
    input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
           Ht2_RNA_ALIGN + Ht2_ANNO_RNA + CLUSTERS + MULTI_QC + CLUSTERS_PLOT

 

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')


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
         "workup/trimmed/{sample}_R1_val_1.fq.gz",
         "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         "workup/trimmed/{sample}_R2_val_2.fq.gz",
         "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o workup/trimmed/ \
        {input} &> {log}"


#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = "workup/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = "workup/trimmed/{sample}_R2_val_2.fq.gz"
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2_barcoded = "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    log:
        "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp("workup/{sample}.ligation_efficiency.txt")
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand("workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"
      

rule full_barcode:
    '''
    remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}_DPM.log"
    shell:
        "python {split_fq} --r1 {input} &> {log}"


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
        fq="workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        mapped="workup/alignments/{sample}.RNA.hisat2.mapq20.bam"
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


rule annotate_rna:
    '''
    -M                  Multi-mapping reads will also be counted. For a multi-
                        mapping read, all its reported alignments will be 
                        counted. The 'NH' tag in BAM/SAM input is used to detect 
                        multi-mapping reads.
    -s <int or string>  Perform strand-specific read counting. A single integer
                        value (applied to all input files) or a string of comma-
                        separated values (applied to each corresponding input
                        file) should be provided. Possible values include:
                        0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                        Default value is 0 (ie. unstranded read counting carried
                        out for all input files).
    -t <string>         Specify feature type in GTF annotation. 'exon' by 
                        default. Features used for read counting will be 
                        extracted from annotation using the provided value.
    -g <string>         Specify attribute type in GTF annotation. 'gene_id' by 
                        default. Meta-features used for read counting will be 
                        extracted from annotation using the provided value.
        '''
    input:
        "workup/alignments/{sample}.RNA.hisat2.mapq20.bam"
    threads: 10
    output:
        bam="workup/alignments/{sample}.RNA.hisat2.mapq20.bam.featureCounts.bam",
        counts="workup/alignments/{sample}.RNA.hisat2.mapq20.bam.featureCounts.txt"
    log:
        "workup/logs/{sample}.anno.log"
    conda:
        "envs/annotate_rna.yaml"
    shell:
        '''
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 1 \
        -g gene_name -a {anno_gtf} -o {output.counts} \
        {input}
        '''


rule add_chr:
    input:
        "workup/alignments/{sample}.RNA.hisat2.mapq20.bam.featureCounts.bam"
    output:
        "workup/alignments/{sample}.RNA.chr.bam"
    log:
        "workup/logs/{sample}.DNA_chr.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {assembly} &> {log}
        '''


rule make_clusters:
    input:
         "workup/alignments/{sample}.RNA.chr.bam"
    output:
        "workup/clusters/{sample}.clusters"
    log:
        "workup/clusters/{sample}.make_clusters.log"
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
        expand("workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
    output:
        "workup/qc/multiqc_report.html"
    log:
        "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc workup -o workup/qc"


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
            DNA.clusters
        '''
