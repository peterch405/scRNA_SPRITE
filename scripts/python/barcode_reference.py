
from collections import defaultdict

#create a fasta with all barcodes and spacers for hisat2 alignment

config = '/mnt/data/scRNA/20191122_HumMus/config.txt'


def parse_config(path):
    '''Parse barcodeID config

    Args:
        path(str): Input config.txt path
    '''

    barcode = ''
    header = []
    barcodes_dct = defaultdict(set)
    present_bc = defaultdict(set)
    with open(config, 'r') as cf:
        for _ in range(10):
            line = cf.readline()
            if line.startswith(('EVEN', 'ODD', 'DPM')):
                tag, name, barcode, errors = line.rstrip().split('\t')
                barcodes_dct[tag].add(barcode)
                break
            else:
                header.append(line)

        for line in cf:
            tag, name, barcode, errors = line.rstrip().split('\t')
            barcodes_dct[tag].add(barcode)

    for item in header:
        if item.startswith('READ2'):
            barcode = item.rstrip('\n').split('= ')[-1]

    barcodes = barcode.split('|')
    for k, v in barcodes_dct.items():
        if k in barcodes:
            present_bc[k] = v

    return((present_bc, barcode))


#96*11 = 1056

bc_dict, barcode = parse_config(config)



def make_fasta(seqs, out, barcodes):
    '''Make a fasta file from a dictionary of sequences

    Args:
        seqs(dict): a dictionary of barcode type (k) and sequences (v)
        out(str): out path for fasta file
        barcodes(str): string of barcoding scheme
    '''

    total_len = 0
    seq = ''
    for v in seqs.values():
        seq += ''.join(v)

    with open(out, 'w') as fa_out:
        fa_out.write('>SPRITE_' + barcodes + '\n')
        len_out = 0
        initial_len = len(seq)
        total_len += initial_len
        while len_out < initial_len:
            fa_out.write(seq[:80] + '\n')
            seq = seq[80:]
            len_out += 80

make_fasta(bc_dict, '/mnt/data/scRNA/barcodes.fasta', barcode)

#make splice site file