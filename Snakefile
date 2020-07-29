'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process scRNA-seq with the SPRITE barcoding scheme
'''


import os 
import sys
import datetime
from collections import defaultdict

#Location of scripts
barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
split_fq = "scripts/python/get_full_barcodes.py"
remove_multi = "scripts/python/remove_multialigners.py"

###################################################################################
#Load config.yaml file
###################################################################################

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path


#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')


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


#Bedtools region filtering (e.g. repeat masking)
try:
    in_filter = config['include_filter']
except:
    in_filter = dict()
    print('Include filtering will be skipped')

try:
    ex_filter = config['exclude_filter']
except:
    ex_filter = dict()
    print('Exclude filtering will be skipped')

################################################################################
# Select aligner index
################################################################################

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly'].split('&')
    assert all(i in ['mm10', 'hg38', 'mm10_rep', 'hg38_rep', 'halo'] for i in assembly), 'Only x options currently supported'
    print('Using assembly:', ' '.join(assembly))
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    rep_assembly = config['repeat_assembly'].split('&')
    assert all(i in ['mm10_rep', 'hg38_rep'] for i in rep_assembly), 'Only x options currently supported'
    print('Using repeat assembly:', ' '.join(rep_assembly))
except:
    print('Config "repeat_assembly" not specified, won\'t align to repeats')
    rep_assembly = ''



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




#Additional hisat2 flags
try:
    h2_add_settings = config["hisat2_settings"]
except:
    h2_add_settings = ''

#Additional bowtie2 flags
try:
    b2_add_settings = config["bowtie2_settings"]
except:
    b2_add_settings = ''

#STAR flags
try:
    STAR_RNA_params = config["STAR_RNA_params"]
except:
    STAR_RNA_params = ''


#Mapq filter
if config['aligner'] == 'star':
    mapq = 255
else:
    try:
        mapq = config["mapq_filter"]
    except:
        mapq = 20


try:
    trim5_len = config['trim5']
    trim3_len = config['trim3']
except:
    trim5_len = 0
    trim3_len = 0

print("Trimming 5' end:", trim5_len)
print("Trimming 3' end:", trim3_len)

################################################################################
# UMI settings
################################################################################

try:
    umi = config['umi']
except:
    umi = False
    print('Running without UMI')

if umi == 'True' and trim5_len < 25:
    print('When UMI is present, 5 prime side should be trimmed at least 25bp')
    sys.exit()


#####################################################################################
#Output files, setting up rule all
#####################################################################################

#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(config['samples']))
ALL_SAMPLES = sorted(FILES.keys())


ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]


TRIM = expand(out_dir + "workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
# TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
#                   sample = ALL_SAMPLES, read = ["R1", "R2"])

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]
MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]
BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.barcoded.fastq.gz", 
                   sample = ALL_SAMPLES, read = ["R1", "R2"])
BARCODE_FULL = expand([out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
                       out_dir + "workup/fastqs/{sample}_R2.barcoded_full.fastq.gz"], 
                       sample=ALL_SAMPLES)

CUTADAPT = expand([out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz",
                   out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz"],
                   sample = ALL_SAMPLES)

FASTQC = expand([out_dir + "workup/trimmed/{sample}_R1_bfull_trim_fastqc.html",
                 out_dir + "workup/trimmed/{sample}_R2_bfull_trim_fastqc.html"],
                 sample = ALL_SAMPLES)

#Hisat2 alignment
Ht2_RNA_ALIGN_SE = expand(out_dir + "workup/alignments/{sample}.{genome}.SE.h2.bam", 
                        sample=ALL_SAMPLES, genome=assembly+rep_assembly)
Ht2_RNA_ALIGN_PE = expand(out_dir + "workup/alignments/{sample}.{genome}.PE.h2.bam", 
                           sample=ALL_SAMPLES, genome=assembly+rep_assembly)

Ht2_ANNO_RNA_SE = expand(out_dir + "workup/alignments/{sample}.{genome}.SE.h2.bam.featureCounts.bam",
                  sample=ALL_SAMPLES, genome=assembly+rep_assembly)
Ht2_ANNO_RNA_PE = expand(out_dir + "workup/alignments/{sample}.{genome}.PE.h2.bam.featureCounts.bam",
                  sample=ALL_SAMPLES, genome=assembly+rep_assembly)

ADD_CHR = expand(out_dir + "workup/alignments/{sample}.{genome}.chr.bam", 
                 sample=ALL_SAMPLES, genome=assembly)

UNIQUE =  expand(out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam", 
                 sample=ALL_SAMPLES, genome=assembly)

CLUSTERS_PLOT = [out_dir + "workup/clusters/cluster_sizes.pdf", 
                 out_dir + "workup/clusters/cluster_sizes.png"]

BOWTIE2_ALIGN = expand(out_dir + "workup/alignments/{sample}.{genome}.DNA.b2.bam", 
                       sample=ALL_SAMPLES, genome=assembly+rep_assembly)

STAR_RNA = expand(out_dir + "workup/alignments/{sample}.{genome}.bam", 
                   sample=ALL_SAMPLES, genome=assembly+rep_assembly)

FLAGSTAT = expand([out_dir + "workup/alignments/{sample}.{genome}.chr.bam.flagstat",
                   out_dir + "workup/alignments/{sample}.{genome}.bam.flagstat",
                   out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam.flagstat",
                   out_dir + "workup/alignments/{sample}.{genome}.filt.bam.flagstat"], 
                   sample=ALL_SAMPLES, genome=assembly)

MAPQ = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.chr.unique.mapq{mapq}.bam",
              sample=ALL_SAMPLES, genome=assembly)

CLUSTERS = expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)

UMI_EXTRACT = expand([out_dir + "workup/trimmed/{sample}_R1_umi_ex.fq.gz",
                      out_dir + "workup/trimmed/{sample}_R2_umi_ex.fq.gz"],
                      sample=ALL_SAMPLES)
UMI_DEDUP = expand(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam",
                   sample=ALL_SAMPLES, genome=assembly)


################################################################################
#Merge
################################################################################

if align_mode == "PE":
    MERGE = expand(out_dir + "workup/merged/{sample}_{read}.fastq.gz", 
                            read = ["R1", "R2"],
                            sample=ALL_SAMPLES)
elif align_mode == "SE":
    MERGE = expand(out_dir + "workup/merged/{sample}.fastq.gz",
                        sample=ALL_SAMPLES)

################################################################################
# Filters
################################################################################

def subset_dict(orig_dict, keys_to_keep):
    '''Subset dictionary by keys to keep
    '''
    new_dict = defaultdict()
    for k, v in orig_dict.items(): 
        if k in keys_to_keep:
            new_dict[k] = v
    return new_dict


#Include and exclude reads filter
if len(ex_filter) > 0 and len(in_filter) > 0:
    first_filt = 'in_filt'
    final_filt = 'in_ex_filt'
    in_filt_assembly = [i for i in assembly if i in list(in_filter.keys())]
    ex_filt_assembly = [i for i in assembly if i in list(ex_filter.keys())]
    
    IN_FILTERED = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{first_filt}.bam", 
                  sample=ALL_SAMPLES, genome=in_filt_assembly)
    EX_FILTERED = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{final_filt}.bam", 
                  sample=ALL_SAMPLES, genome=ex_filt_assembly)
    efa = subset_dict(ex_filter, ex_filt_assembly)
    ifa = subset_dict(in_filter, in_filt_assembly)
#Exclude filter only
elif len(in_filter) > 0 and len(ex_filter) == 0:
    first_filt = 'in_filt'
    in_filt_assembly = [i for i in assembly if i in list(in_filter.keys())]

    IN_FILTERED = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{first_filt}.bam", 
                  sample=ALL_SAMPLES, genome=in_filt_assembly)
    EX_FILTERED = []
    ifa = subset_dict(in_filter, in_filt_assembly)
#Exclude filter only
elif len(ex_filter) > 0 and len(in_filter) == 0:
    final_filt = 'ex_filt'
    ex_filt_assembly = [i for i in assembly if i in list(ex_filter.keys())]

    EX_FILTERED = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{final_filt}.bam", 
                  sample=ALL_SAMPLES, genome=ex_filt_assembly)
    IN_FILTERED = []
    efa = subset_dict(ex_filter, ex_filt_assembly)

else:
    first_filt = ''
    final_filt = ''
    IN_FILTERED = []
    EX_FILTERED = []
    in_filt_assembly = ''
    ex_filt_assembly = ''
    # print('Filtering rules not working correctly')
    # sys.exit()

################################################################################
# Select which output files
################################################################################

if config['aligner'] == 'hisat2':
    if align_mode == 'SE':
        VARIABLE = Ht2_RNA_ALIGN_SE #+ Ht2_ANNO_RNA_SE
    elif align_mode == 'PE':
        VARIABLE = Ht2_RNA_ALIGN_PE #+ Ht2_ANNO_RNA_PE
elif config['aligner'] == 'bowtie2': 
    VARIABLE = BOWTIE2_ALIGN 
elif config['aligner'] == 'star':
    VARIABLE = STAR_RNA

if umi == 'True':
    VARIABLE = VARIABLE + UMI_EXTRACT + UMI_DEDUP

########################################################################################################################
#Main rule
########################################################################################################################

rule all:
    input: CONFIG + ALL_FASTQ + TRIM + BARCODEID + LE_LOG_ALL + BARCODE_FULL +
           CUTADAPT + VARIABLE + ADD_CHR + UNIQUE + IN_FILTERED + EX_FILTERED + CLUSTERS + 
           MULTI_QC + CLUSTERS_PLOT + MAPQ #FLAGSTAT

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')

################################################################################
#Merge
################################################################################
#For files that don't need to be merged, just create a symbolic link
rule merge_fastqs_se:
    wildcard_constraints:
        sample='.+(?<!R[12])'
    input:
        lambda wildcards: FILES[wildcards.sample]['R1']
    output: 
        out_dir + "workup/merged/{sample}.fastq.gz"
    shell: 
        '''
        count=$(echo '{input}' | awk -F' ' '{{print NF}}')
        
        if [[ $count -gt 1 ]]
        then
            cat {input} > {output}
        else
            ln -s {input} {output}
        fi
        '''

rule merge_fastqs_pe:
    input: 
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        r1 = out_dir + "workup/merged/{sample}_R1.fastq.gz",
        r2 = out_dir + "workup/merged/{sample}_R2.fastq.gz"
    shell:
        ''' 
        count_1=$(echo '{input.r1}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_1 -gt 1 ]]
        then
            cat {input.r1} > {output.r1}
        else
            ln -s {input.r1} {output.r1}
        fi

        count_2=$(echo '{input.r2}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_2 -gt 1 ]]
        then
            cat {input.r2} > {output.r2}
        else
            ln -s {input.r2} {output.r2}
        fi
        '''

########################################################################################################################
#Trimming and barcode identification
########################################################################################################################

rule adaptor_trimming_pe:
    '''
    Trim adaptors
    multiple cores requires pigz to be installed on the system
    '''
    input:
        [out_dir + "workup/merged/{sample}_R1.fastq.gz", 
        out_dir + "workup/merged/{sample}_R2.fastq.gz"] 
    output:
        out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
        out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz"
    log:
        out_dir + "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        '''
        trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''


rule log_config:
    '''Copy config.yaml and place in logs folder with the date run
    '''
    input:
        config_path
    output:
        out_dir + "workup/logs/config_" + run_date + "yaml"
    shell:
        "cp {input} {output}"




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
    conda:
        "envs/python_dep.yaml"
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
    -n 3  Remove up to COUNT adapters from each read
    -e 0.15 Maximum allowed error rate as value between 0 and 1
            (no. of errors divided by length of matching region).
            Default: 0.1 (=10%)
    Read 1:
    -a 3' adaptor
    -g 5' adaptor
    Read 2:
    -A 3' adaptor
    -G 5' adaptor
    in r2 TGACTTGNTTCG
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
        adapters = "-a CAAGTCA -a AGTTGTC" if align_mode == "SE" else "-a CAAGTCA -a AGTTGTC -G TGACTTG",
        others = "-n 2 --minimum-length 15"
    log:
        out_dir + "workup/logs/{sample}.cutadapt.log"
    wrapper:
        "0.63.0/bio/cutadapt/pe"
# r"(?e)(CAAGTCA){e<=1}|(?e)(AGTTGTC){e<=1}"


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


########################################################################################################################
#UMI extraction
########################################################################################################################

rule umi_extract:
    input:
        fastq1=out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz",
        fastq2=out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz"
    output:
        fastq1=out_dir + "workup/trimmed/{sample}_R1_umi_ex.fq.gz",
        fastq2=out_dir + "workup/trimmed/{sample}_R2_umi_ex.fq.gz"
    log:
        out_dir + "workup/logs/{sample}.umi_extract.log"
    conda:
        "envs/umi_tools.yaml"
    shell:
        '''
        umi_tools extract \
        -I {input.fastq1} \
        --bc-pattern=NNNNNNNNNN \
        --read2-in={input.fastq2} \
        --stdout={output.fastq1} \
        --read2-out={output.fastq2} 2> {log}
        '''
        

rule umi_deduplicate:
    input:
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam"
    output:
        out=out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam",
        idx=temp(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.sorted.bam.bai"),
        sort=temp(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.sorted.bam")
    log:
        out_dir + "workup/logs/{sample}.{genome}.umi_dedup.log"
    params:
        pair = "--paired" if align_mode == "PE" else "",
        stats = out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.dedup.stats"
    threads:
        10
    conda:
        "envs/umi_tools.yaml"
    shell:
        '''
        samtools sort -T {wildcards.sample} -@ {threads} -o {output.sort} {input}
        samtools index {output.sort}

        umi_tools dedup -I {output.sort} {params.pair} --output-stats={params.stats} -S {output.out} 2> {log}
        '''

########################################################################################################################
#RNA alignment
########################################################################################################################



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
       out_dir + "workup/trimmed/{sample}_R1_umi_ex.fq.gz" if umi=='True' else 
       out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz"
    output:
        mapped=out_dir + "workup/alignments/{sample}.{genome}.SE.h2.bam"
    threads: 10
    params:
        ss = lambda wildcards: hisat2_ss[wildcards.genome],
        index = lambda wildcards: hisat2_index[wildcards.genome]
    conda:
        "envs/hisat2.yaml"
    log:
        out_dir + "workup/logs/{sample}.{genome}.SE.hisat2.log"
    shell:
        '''
        (hisat2 \
        -p {threads} \
        -t \
        --phred33 \
        {h2_add_settings} \
        --trim5 {trim5_len} \
        --trim3 {trim3_len} \
        --known-splicesite-infile {params.ss} \
        -x {params.index} \
        -U {input} | \
        samtools view -b -F 4 -F 256 - > {output.mapped}) &> {log}
        '''
# q {params.mq_filter}
 # params:
        # index = lambda wildcards, output, input: if 'hg38'

def hisat2_align_pe_input(wildcards):
    '''Input to hisat2_align_pe
    '''
    if umi == 'True':
        fq_1=out_dir + "workup/trimmed/{sample}_R1_umi_ex.fq.gz"
        fq_2=out_dir + "workup/trimmed/{sample}_R2_umi_ex.fq.gz"
    else:
        fq_1=out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz"
        fq_2=out_dir + "workup/trimmed/{sample}_R2_bfull_trim.fq.gz"
    return {'fq_1':fq_1, 'fq_2':fq_2}

  
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
        hisat2_align_pe_input
    output:
        out_dir + "workup/alignments/{sample}.{genome}.PE.h2.bam"
    threads: 10
    params:
        ss = lambda wildcards: hisat2_ss[wildcards.genome],
        index = lambda wildcards: hisat2_index[wildcards.genome]
    conda:
        "envs/hisat2.yaml"
    log:
        out_dir + "workup/logs/{sample}.{genome}.PE.hisat2.log"
    shell:
        ''' 
        (hisat2 \
        --dta \
        -p {threads} \
        -t \
        --phred33 \
        {h2_add_settings} \
        --trim5 {trim5_len} \
        --trim3 {trim3_len} \
        --known-splicesite-infile {params.ss} \
        -x {params.index} \
        -1 {input['fq_1']} -2 {input['fq_2']} | \
        samtools view -b -F 4 -F 256 - > {output}) &> {log}
        '''
# --known-splicesite-infile {params.ss}
# q {params.mq_filter}
# --sp 1000,1000 \
# --no-mixed \
# --no-discordant \


rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq=out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz"
    output:
        out_dir + "workup/alignments/{sample}.DNA.b2.bam"
    threads: 10
    params:
        bowtie2_idx = lambda wildcards: bowtie2_index[wildcards.genome],
        # mq_filter = mapq
    log:
        out_dir + "workup/logs/{sample}.bowtie2.log"       
    conda:
        "envs/bowtie2.yaml"
    shell:
        "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        {b2_add_settings} \
        --trim5 {trim5_len} \
        --trim3 {trim3_len} \
        -x {params.bowtie2_idx} \
        -U {input.fq} | \
        samtools view -b -F 4 -F 256 - > {output}) &> {log}"





rule star_align_rna:
    '''
    Align RNA with STAR to the genome first, annotate repeats, 
    anything that does not align, realign with bowtie2 to our repeats reference
    '''
    input:
        out_dir + "workup/trimmed/{sample}_R1_umi_ex.fq.gz" if umi=='True' else 
        out_dir + "workup/trimmed/{sample}_R1_bfull_trim.fq.gz"
    output:
        out_dir + "workup/alignments/{sample}.{genome}.bam",
        out_dir + "workup/alignments/{sample}.{genome}.Aligned.sortedByCoord.out.bam",
        out_dir + "workup/alignments/{sample}.{genome}.Log.final.out",
        out_dir + "workup/alignments/{sample}.{genome}.SJ.out.tab",
        out_dir + "workup/alignments/{sample}.{genome}.unmapped.fastq.gz"
    wildcard_constraints:
        genome = ['hg38', 'mm10']
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
        --readFilesIn {input} \
        --outFileNamePrefix {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}. &> {log}

        mv {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.Unmapped.out.mate1 \
            {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.unmapped.fastq

        pigz {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.unmapped.fastq

        '''
   # samtools view -bq 255 {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.Aligned.sortedByCoord.out.bam > \
        #     {out_dir}workup/alignments/{wildcards.sample}.{wildcards.genome}.mapq20.bam

########################################################################################################################
#Read annotation
########################################################################################################################

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
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam"
    threads: 10
    output:
        bam=out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam.featureCounts.bam",
        counts=out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam.featureCounts.txt"
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


########################################################################################################################
#Filtering and cluster making
########################################################################################################################

rule filter_include_reads:
    '''
    Filter, keeping only reads within protein coding genes. (Alternatively, Repeatmask filter)
    -v 	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
    '''
    input:
        out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam" if umi == 'True' else
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam"
    output:
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.{first_filt}.bam"
    wildcard_constraints:
        genome='|'.join(in_filt_assembly)
    params:
        keep_regions = lambda wildcards: ifa[wildcards.genome]
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -a {input} -b {params.keep_regions} > {output}
        '''


def filter_exclude_reads_input(wildcards):
    '''Input to filter_exclude_reads
    '''
    
    if wildcards.genome in list(in_filter.keys()):
        filter_in = [wildcards.genome]
    else:
        fitler_in = []
    if wildcards.genome in list(ex_filter.keys()):
        filter_ex = [wildcards.genome]
    else:
        filter_ex = []
 
    both_filt = set(filter_in) & set(filter_ex) #intersection
    # in_only = set(filter_in) - set(filter_ex) #difference
    ex_only = set(filter_ex) - set(filter_in) #difference
    # no_filt = set(assembly) - (both_filt | in_only | ex_only) #no filter
    
    if len(both_filt) > 0:
        BOTH = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{first_filt}.bam",
                      sample=wildcards.sample, genome=both_filt)
    else:
        BOTH = []
    if len(ex_only) > 0:
        if umi == 'True':
            EX = expand(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam",
                          sample=wildcards.sample, genome=ex_only)
        else:
            EX = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam",
                          sample=wildcards.sample, genome=ex_only)
    else:
        EX = []

    return BOTH + EX


rule filter_exclude_reads:
    '''
    Filter, keeping only reads within protein coding genes. (Alternatively, Repeatmask filter)
    -v 	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
    '''
    input:
        # rules.filter_include_reads.output if len(IN_FILTERED) > 0 else 
        # out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam"
        filter_exclude_reads_input
    output:
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.{final_filt}.bam"
    wildcard_constraints:
        genome='|'.join(ex_filt_assembly)
    params:
        remove_regions = lambda wildcards: efa[wildcards.genome]
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {params.remove_regions} > {output}
        '''



def add_chr_input(wildcards):
    '''Input to add_chr
    If no filtering, use output from aligner
    '''
    
    if wildcards.genome in list(in_filter.keys()):
        filter_in = [wildcards.genome]
    else:
        filter_in = []
    if wildcards.genome in list(ex_filter.keys()):
        filter_ex = [wildcards.genome]
    else:
        filter_ex = []
 
    both_filt = set(filter_in) & set(filter_ex) #intersection
    in_only = set(filter_in) - set(filter_ex) #difference
    ex_only = set(filter_ex) - set(filter_in) #difference
    no_filt = set([wildcards.genome]) - (both_filt | in_only | ex_only) #no filter

    if len(both_filt) > 0:
        BOTH = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{final_filt}.bam",
                      sample=wildcards.sample, genome=both_filt)
    else:
        BOTH = []
    if len(ex_only) > 0:
        EX = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{final_filt}.bam",
                      sample=wildcards.sample, genome=ex_only)
    else:
        EX = []
    if len(in_only) > 0:
        IN =  expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{first_filt}.bam",
                    sample=wildcards.sample, genome=in_only)
    else:
        IN = []
    if len(no_filt) > 0:
        if umi == 'True':
            NO = expand(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam",
                        sample=wildcards.sample, genome=no_filt)
        else:
            NO = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam",
                        sample=wildcards.sample, genome=no_filt)
    else:
        NO = []
    
    return BOTH + EX + IN + NO
    

rule add_chr:
    input:
        add_chr_input
    output:
        out_dir + "workup/alignments/{sample}.{genome}.chr.bam"
    wildcard_constraints:
        genome='|'.join(assembly)
    log:
        out_dir + "workup/logs/{sample}.{genome}.chr.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {wildcards.genome} &> {log}
        '''


def remove_multi_assembly_aligners_input(wildcards):
    '''Input to rule remove_multi_assembly_aligners
    '''
    if umi == 'True':
        if len(rep_assembly) > 0:
            rep = expand(out_dir + f"workup/deduplicated/{{sample}}.{{genome}}.{align_mode}.umi_dedup.bam", 
                         sample=ALL_SAMPLES, genome=rep_assembly)
        else:
            rep = ""

    else:
        if len(rep_assembly) > 0:
            rep = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam", 
                         sample=ALL_SAMPLES, genome=rep_assembly)
        else:
            rep = ""

    main = expand(out_dir + "workup/alignments/{sample}.{genome}.chr.bam", 
               sample=ALL_SAMPLES, genome=assembly)

    return rep + main


rule remove_multi_assembly_aligners:
    input:
        # main = expand(out_dir + "workup/alignments/{sample}.{genome}.chr.bam", 
        #        sample=ALL_SAMPLES, genome=assembly),
        # rep = expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam", 
        #        sample=ALL_SAMPLES, genome=rep_assembly) if len(rep_assembly) > 0 else "" 
        remove_multi_assembly_aligners_input
    output:
        out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {remove_multi} -i {input} &> {out_dir}workup/logs/multi_assembly_align.log
        '''



# rule flagstat_bam:
#     '''
#     Check number of reads mapped with samtools flagstat
#     This is so the alignment stats will be included in multiqc
#     '''
#     input: 
#         bam=out_dir + f"workup/alignments/{{sample}}.{{genome}}.{align_mode}.h2.bam",
#         chr=out_dir + "workup/alignments/{sample}.{genome}.chr.bam",
#         unq=out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam",
#         filt=out_dir + "workup/alignments/{sample}.{genome}.filt.bam"
#     output: 
#         bam=out_dir + "workup/alignments/{sample}.{genome}.bam.flagstat",
#         chr=out_dir + "workup/alignments/{sample}.{genome}.chr.bam.flagstat",
#         unq=out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam.flagstat",
#         filt=out_dir + "workup/alignments/{sample}.{genome}.filt.bam.flagstat"
#     log:
#         bam=out_dir + "workup/logs/{sample}.{genome}.flagstat_bam",    
#         chr=out_dir + "workup/logs/{sample}.{genome}.chr.flagstat_bam",
#         unq=out_dir + "workup/logs/{sample}.{genome}.chr.unique.flagstat_bam",
#         filt=out_dir + "workup/logs/{sample}.{genome}.filt.flagstat_bam"
#     threads: 7
#     message: "flagstat_bam {input}: {threads} threads"
#     shell:
#         '''
#         samtools flagstat -@ {threads} {input.bam} > {output.bam} 2> {log.bam}
#         samtools flagstat -@ {threads} {input.chr} > {output.chr} 2> {log.chr}
#         samtools flagstat -@ {threads} {input.unq} > {output.unq} 2> {log.unq}
#         samtools flagstat -@ {threads} {input.filt} > {output.filt} 2> {log.filt}
#         '''

rule mapq_filter:
    input:
        out_dir + "workup/alignments/{sample}.{genome}.chr.unique.bam"
    output:
        out_dir + f"workup/alignments/{{sample}}.{{genome}}.chr.unique.mapq{mapq}.bam"
    conda:
        "envs/bowtie2.yaml"
    log: 
        out_dir + "workup/logs/{sample}.{genome}.mapq_filt.log"
    params:
        mapq_filt = mapq
    shell:
        '''
        samtools view -bq {params.mapq_filt} {input} > {output} 2> {log}
        '''


rule make_clusters:
    input:
        expand(out_dir + f"workup/alignments/{{sample}}.{{genome}}.chr.unique.mapq{mapq}.bam", 
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
        "multiqc {out_dir}workup -o {out_dir}workup/qc"


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
            "\\.clusters"
        '''


# rule gene_body_coverage:
#     input:

#     output:
        
#     conda:
#         "envs/qc.yaml"

# geneBody_coverage.py -r hg19.housekeeping.bed -i test.bam  -o output
