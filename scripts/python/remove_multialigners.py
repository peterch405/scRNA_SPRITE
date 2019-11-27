import pysam
from collections import defaultdict
import os
import argparse
#if reads align to both human and mouse, exclude these 


def main():

    opts = parse_arguments()

    # hg = '/mnt/data/scRNA/20191122_HumMus/workup/alignments/scRNAHumMus.hg38.chr.bam'
    # mm = '/mnt/data/scRNA/20191122_HumMus/workup/alignments/scRNAHumMus.mm10.chr.bam'
    # files = [hg, mm]

    #if only a single file present, just rename to unique

    if len(opts.input) <= 1:
        out_name = opts.input.replace('.bam', '') + '.unique.bam'
        os.rename(opts.input, out_name)
    else:
        mapped_reads_dict = defaultdict(set)
        for f in opts.input:
            mapped_reads_dict[f] = mapped_reads(f)

        unique_mapped = find_unique(mapped_reads_dict)
        write_unique(unique_mapped)


def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "This program removes reads from a BAM file if they are present " +
            "in another bam.")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',nargs='+',
                        required=True, help = 'Input BAM(s) file')
    return parser.parse_args()


def mapped_reads(bam_path):
    '''Get read name or mapped reads

    Args:
        bam_path(str): Path to bam file
    '''
    mapped = set()
    with pysam.AlignmentFile(bam_path, 'rb') as input_file:
        for read in input_file.fetch(until_eof = True):
            if not read.is_secondary and not read.is_unmapped:
                mapped.add(read.query_name)

    return(mapped)


def find_unique(sets_of_reads):
    '''Find unique reads for each dataset

    Args:
        sets_of_reads(dict): A file (k) dict of read names (v)
    '''
    datasets = set(sets_of_reads.keys())
    unique_sets_of_reads = defaultdict(set)

    for k, v in sets_of_reads.items():
        the_rest = datasets.difference(set([k]))
        other_reads = set()
            
        for k2,v2 in sets_of_reads.items():
            if k2 in the_rest:
                other_reads.update(v2)

        unique_sets_of_reads[k] = v.difference(other_reads)
        
    return(unique_sets_of_reads)


# test_find_unique = {'data1':set(['a','b','c']), 'data2':set(['a','d','e']), 
#                     'data3':set(['a', 'e','f'])}
# find_unique(test_find_unique)


def write_unique(unique_sets_of_reads):
    '''Write out bam files with only uniquely mapped reads
    '''
    for k, v in unique_sets_of_reads.items():
        out_file = k.replace('.bam', '') + '.unique.bam'
        total = 0
        out_reads = 0
        with pysam.AlignmentFile(k, 'rb') as input_file:
            output_file = pysam.AlignmentFile(out_file, "wb", template = input_file)
            for read in input_file.fetch(until_eof = True):
                total += 1
                if read.query_name in v:
                    output_file.write(read)
                    out_reads += 1

        print(out_reads, 'written for', k, 'out of a total of', total)


if __name__ == "__main__":
    main()
