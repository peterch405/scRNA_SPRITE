import argparse
import cluster as c

def main():
    args = parse_arguments()
    clusters = c.parse_cluster(args.input)
    c.write_bam(clusters, args.num_tags, args.bam, args.output)
    print("done")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a clusters file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The input cluster file.")
    parser.add_argument('-b', '--bam',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The BAM file used to create cluster file, \
                                or with reads of interest.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The output BAM file.")
    parser.add_argument('-n', '--num_tags',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        required=True,
                        help = "The number of tags contained in the barcode " +
                               "of each BAM record.")

    return parser.parse_args()

if __name__ == "__main__":
    main()
