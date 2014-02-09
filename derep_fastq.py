# pypy script <in_fastq> <out_fastq>

import sys
import gzip
from libs.bio_file_parsers import fastq_parser

def main(args):
    
    seqs = {}
    dup = 0
    
    for record in fastq_parser(gzip.open(args['in'], 'r')):
        if record[1] not in seqs:
            seqs[record[1]] = {'title':record[0],
                               'qual':record[2],
                               'size':1}
        else:
            seqs[record[1]]['size'] += 1
            dup += 1
    
    with open(args['out'], 'w') as in_handle:
        for key in sorted(seqs.keys(), key=lambda x: seqs[x]['size'], reverse=True):
            title = seqs[key]['title'].split(' ')[0]
            title = '@{0};size={1}'.format(title, seqs[key]['size'])
            entry = [title, key, '+', seqs[key]['qual']]
            in_handle.write('\n'.join(entry) + '\n')
    
    print len(seqs)
    print len(seqs) + dup

if __name__ == '__main__':
    
    args = {}
    
    args['in'] = sys.argv[1]
    args['out'] = sys.argv[2]
    
    main(args)
