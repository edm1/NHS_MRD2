# pypy script <in_fastq> <out_fastq>

import sys
from libs.bio_file_parsers import fasta_parser
from libs.bio_file_parsers import write_fasta

def main(args):
    
    seqs = {}
    dup = 0
    
    for record in fasta_parser(open(args['in'], 'r')):
        if record[1] not in seqs:
            seqs[record[1]] = {'title':record[0],
                               'size':1}
        else:
            seqs[record[1]]['size'] += 1
            dup += 1
    
    with open(args['out'], 'w') as out_handle:
        for key in sorted(seqs.keys(), key=lambda x: seqs[x]['size'], reverse=True):
            parts = seqs[key]['title'].split(';size=')
            title = parts[0]
            size = int(parts[1]) + seqs[key]['size']
            title = '{0};size={1}'.format(title, size)
            write_fasta(out_handle, title, key)
    
    print len(seqs)
    print len(seqs) + dup

if __name__ == '__main__':
    
    args = {}
    
    args['in'] = sys.argv[1]
    args['out'] = sys.argv[2]
    
    main(args)
