# pypy script <in_fasta> <out_fasta>

import re
import sys
from operator import itemgetter
from libs.bio_file_parsers import write_fasta, fasta_parser

def main(args):
    update_sizes(args['in_fasta'], args['in_clstr'])
    

def update_sizes(fasta, clstr):
    
    # Load fasta into memory
    with open(fasta, 'r') as in_fasta:
        fasta_records = [x for x in fasta_parser(in_fasta)]
    
    # Parse cluster sizes
    clus_sizes = parse_clus_sizes(clstr)
    
    # Get total cluster size
    total_clusters_size = 0
    for key in clus_sizes:
        total_clusters_size += clus_sizes[key]
    
    # Zip records and sizes together
    fasta_records = [(x[0], x[1], clus_sizes[x[0]]) for x in fasta_records]
    
    # Write new fasta
    with open(fasta, 'w') as out_handle:
        for record in sorted(fasta_records, key=itemgetter(2), reverse=True):
            # Replace new size in header
            new_header = '{0};size={1}'.format(record[0].split(';size=')[0],
                                               record[2])
            write_fasta(out_handle, new_header, record[1])
    
    return total_clusters_size
    
def parse_clus_sizes(clstr):
    """ Makes a dictionary with centroid name as key and cluster size as value.
    """
    name_pat = re.compile(r'.*>(.*)\.{3}.*')
    size_pat = re.compile('.*size=([0-9]+).*')
    clus_sizes = {}
    
    with open(clstr, 'r') as in_handle:
        for clus_num, info in clstr_parser(in_handle):
            clus_size = 0
            for entry in info:
                # Mark as centroid if ends in *
                if entry.endswith('*'):
                    centroid = name_pat.search(entry).group(1)
                # Add read size to total cluster size
                clus_size += int(size_pat.search(entry).group(1))
            # Save to dictionary
            clus_sizes[centroid] = clus_size
    
    return clus_sizes

def clstr_parser(handle):
    """ Parses a cd-hit-est clstr file. Very similar to parsing a fasta file.
    """

    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        yield title, lines

        if not line:
            return  # StopIteration

if __name__ == '__main__':
    
    args = {}
    args['in_fasta'] = sys.argv[1]
    args['in_clstr'] = sys.argv[2]
    
    main(args)
