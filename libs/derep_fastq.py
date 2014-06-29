# pypy script <in_fastq> <out_fastq>

import sys
import gzip
from libs.bio_file_parsers import fastq_parser
#~ from bio_file_parsers import fastq_parser
from operator import itemgetter

def main(args):
    
    derep_fastq(args['in'], args['out'])

def derep_fastq(in_file, out_file):
    
    seqs = {}
    dup = 0
    
    # Open correct handle for file compression
    if in_file.split('.')[-1] == 'gz':
        in_handle = gzip.open(in_file, 'r')
    else:
        in_handle = open(in_file, 'r')
    
    # Save records to a dictionary, using the seq as a key
    for record in fastq_parser(in_handle):
        if record[1] not in seqs:
            seqs[record[1]] = {'title':record[0],
                               'qual':record[2],
                               'size':1}
        else:
            seqs[record[1]]['size'] += 1
            dup += 1
    
    # Close file
    in_handle.close()
    
    # Write one copy of each to a new file
    with open(out_file, 'w') as in_handle:
        for key in sorted_keys(seqs):
            title = seqs[key]['title'].split(' ')[0]
            title = '@{0};size={1}'.format(title, seqs[key]['size'])
            entry = [title, key, '+', seqs[key]['qual']]
            in_handle.write('\n'.join(entry) + '\n')
    
    del seqs, dup
    
    return 1

def sorted_keys(seq_dict):
    pairs = []
    for key in seq_dict.keys():
        pairs.append((key, seq_dict[key]['size']))
    for key, size in sorted(pairs, key=itemgetter(1), reverse=True):
        yield key

if __name__ == '__main__':
    
    args = {}
    
    args['in'] = sys.argv[1]
    args['out'] = sys.argv[2]
    
    main(args)
