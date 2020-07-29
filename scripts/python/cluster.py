
import pysam
import re
import gzip
import os
import sys
from collections import defaultdict

class Position:
    """This class represents a genomic position, with type of nucleic acid (DNA)

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "DNA[feature]_chrX:1000-1500"
    """

    def __init__(self, type, feature, chromosome, start_coordinate, end_coordinate):
        self._type = type
        self._feature = feature
        self._chromosome = chromosome
        self._start_coordinate = start_coordinate
        self._end_coordinate = end_coordinate

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return (self._type == other._type and
                self._feature == other._feature and
                self._chromosome == other._chromosome and
                self._start_coordinate == other._start_coordinate and
                self._end_coordinate == other._end_coordinate)

    def __hash__(self):
        return hash((self._type, self._feature, self._chromosome, 
                     self._start_coordinate, self._end_coordinate))

    def to_string(self):
        try:
            out = self._type + "[" + self._feature + "]" + "_" + \
                  self._chromosome + ":" + \
                  str(self._start_coordinate) + "-" + str(self._end_coordinate)
        except:
            print(self._type, self._feature, self._chromosome)
            print('Elements are not as expect!')
            sys.exit()
        return out

class Cluster:
    """This class represents a barcoding cluster as a collection of genomic
    positions.

    The underlying data structure is a set, so duplicate positions are
    discarded.

    Methods:
    - add_position(position): Adds a genomic position to this cluster

    - size(): Returns the number of reads or positions in this cluster

    - to_string(): Returns a string representation of this cluster as a
      tab-delimtited series of positions. See Position#to_string for how
      positions are represented as strings.
    """

    def __init__(self):
        self._positions = set()

    def add_position(self, position):
        self._positions.add(position)

    def size(self):
        return len(self._positions)

    def to_string(self):
        position_strings = [position.to_string() for position in self._positions]
        return "\t".join(position_strings)
    
    def to_list(self):
        position_strings = [position.to_string() for position in self._positions]
        return position_strings



class Clusters:
    """This class represents a collection of barcoding clusters.

    Methods:
    - get_cluster(barcode): Returns the cluster that corresponds to the given
      barcode. If the cluster does not exist, it is initialized (with zero
      positions), and this empty cluster is returned.

    - add_position(barcode, position): Adds the position to the cluster
      that corresponds with the given barcodes

    - to_strings(): Returns an iterator over the string representations of all
      of the contained clusters.

    - remove_cluster(barcode): Removes a cluster with the specified barcode

    - unique(): 
    """
    def __init__(self):
        self._clusters = {}

    def get_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]

    def add_position(self, barcode, position):
        self.get_cluster(barcode).add_position(position)

    def to_strings(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.to_string()

    def remove_cluster(self, barcode):
        del self._clusters[barcode]

    def make_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            lookup[barcode].update(cluster.to_list())
        return lookup 



def get_clusters(bamfile, num_tags):
    """Parses a BAM file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """
    #strip RPM DPM from barcode
    #TODO add file name as an additional barcode
    
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    for bam in bamfile:
        #get sample name from bamfile
        file_name = os.path.basename(bam)
        #get genome build
        if 'hg38' in file_name:
            assembly = 'hg38'
        elif 'mm10' in file_name:
            assembly = 'mm10'

        sample_name = file_name.split('.')[0]
        try:
            with pysam.AlignmentFile(bam, "rb") as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    strand = '+' if not read.is_reverse else '-'
                    umi = name.split('_')[-1]

                    if read.has_tag('XT'):
                        gene_anno = read.get_tag('XT')
                    elif read.has_tag('XS'):
                        gene_anno = read.get_tag('XS')
                    else:
                        gene_anno = ''
                    
                    
                    anno = ';'.join(filter(None, [assembly, gene_anno, strand, umi]))

                    position = Position('DNA', anno, read.reference_name,
                                        read.reference_start, read.reference_end)
                    barcode.append(sample_name)
                    barcode_str = ".".join(barcode)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print('BAM file provided is not a BAM or is empty!')

    return clusters



def write_clusters_to_file(clusters, outfile, unique=False):
    """Writes a Clusters object to a file"""

    count = 0
    with open(outfile, 'w') as f:

        if unique:
            for cluster_string in clusters.unique():
                f.write(cluster_string)
                f.write("\n")
                count += 1
        else:
            for cluster_string in clusters.to_strings():
                f.write(cluster_string)
                f.write("\n")
                count += 1
    print('Number of clusters written:',count)



def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f


def parse_cluster(c_file):
    '''
    Parse cluster file

    Args:
        c_file(str): input path of cluster file
    '''

    total_reads = 0
    clusters = Clusters()
    pattern = re.compile('([a-zA-Z0-9]+)\[([a-zA-Z0-9_;\-\+]+)\]_([a-zA-Z0-9_\-]+):([0-9]+)\-([0-9]+)')
    # match = pattern.search('DNA[hg38;+]_chrX:78818752-78819946')
    with file_open(c_file) as c:
        for line in c:

            barcode, *reads = line.decode('utf-8').rstrip('\n').split('\t')

            for read in reads:
                total_reads += 1
                try:
                    match = pattern.search(read)
                    n_type, anno, chrom, start, end = match.groups()
                    position = Position(n_type, anno, chrom, start, end)
                    clusters.add_position(barcode, position)
                except:
                    print(read)
                    raise Exception('Pattern did not match above printed string')
    print('Total cluster reads:', total_reads)
    return(clusters)


# clusters_test = parse_cluster('/mnt/data/scRNA/20191205_scrna_mus_hs_2/workup/clusters/scrna-mus-hs-2_S1_L001.clusters')




def write_bam(cluster, num_tags, original_bam, output_bam):
    '''From a cluster make/subset bam file
    If barcode, chrom, start and end coordinates match, write out read into new BAM

    Args:
        cluster(Clusters): 
    '''

    #get sample name from bamfile
    file_name = os.path.basename(original_bam)
    sample_name = file_name.split('.')[0]

    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    #get genome build
    if 'hg38' in file_name:
        assembly = 'hg38'
    elif 'mm10' in file_name:
        assembly = 'mm10'

    read_lookup = cluster.make_lookup()
    out_reads = 0
    with pysam.AlignmentFile(original_bam, "rb") as in_bam:
        out_bam = pysam.AlignmentFile(output_bam, "wb", template = in_bam)
        for read in in_bam.fetch(until_eof = True):
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())
            strand = '+' if not read.is_reverse else '-'
                        

            if read.has_tag('XT'):
                gene_anno = read.get_tag('XT')
            elif read.has_tag('XS'):
                gene_anno = read.get_tag('XS')
            else:
                gene_anno = ''
                    
                    
            anno = ';'.join(filter(None, [assembly, gene_anno, strand]))
            barcode.append(sample_name)
            barcode_str = ".".join(barcode)

            position = Position('DNA', anno, read.reference_name,
                                read.reference_start, read.reference_end)

            barcode_reads = read_lookup.get(barcode_str, "Not present")
         
            if position.to_string() in barcode_reads:
                out_bam.write(read)
                out_reads += 1 

    out_bam.close()
    print('Total reads written:', out_reads)