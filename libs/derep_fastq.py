# pypy script <in_fastq> <out_fastq>

import sys
import gzip
from libs.bio_file_parsers import fastq_parser
#~ from bio_file_parsers import fastq_parser
from operator import itemgetter

def main(args):
    
    derep_fastq(args['in'], args['out'])

def derep_fastq(in_file, out_file):
    
    phred_dict = phred_score_dict(33)
    
    
    # Open correct handle for file compression
    if in_file.split('.')[-1] == 'gz':
        in_handle = gzip.open(in_file, 'r')
    else:
        in_handle = open(in_file, 'r')
    
    seqs = {}
    total_reads = 0
    dup = 0
    
    # Save records to a dictionary, using the seq as a key
    for q_header, q_seq, q_qual in fastq_parser(in_handle):
    
        total_reads += 1
        
        # Calc average quality
        q_ave_qual = average_phred_score(q_qual, phred_dict)
        
        # Check if seq exists
        if q_seq not in seqs:
            seqs[q_seq] = [q_header, q_ave_qual, 1]
        else:
            dup += 1
            # Add 1 to the count
            seqs[q_seq][2] += 1
            # If the quality is higher then replace
            if q_ave_qual > seqs[q_seq][1]:
                seqs[q_seq][0] = q_header
                seqs[q_seq][1] = q_ave_qual
    
    in_handle.close()
    
    # Make a dict of the uniqe reads headers and their cluster sizes
    to_keep = {}
    for seq in seqs:
        to_keep[seqs[seq][0]] = seqs[seq][2]
    
    ### Iterate over the fastq again and keep unique records ##
    
    # Open correct handle for file compression
    if in_file.split('.')[-1] == 'gz':
        in_handle = gzip.open(in_file, 'r')
    else:
        in_handle = open(in_file, 'r')
    
    # Write new fastq
    with open(out_file, 'w') as out_handle:
        for header, seq, qual in fastq_parser(in_handle):
            if not header in to_keep:
                continue
            title = header.split(' ')[0]
            title = '@{0};size={1}'.format(title, to_keep[header])
            entry = [title, seq, '+', qual]
            out_handle.write('\n'.join(entry) + '\n')

    return 0

def average_phred_score(qual_string, phred_dict):
    """ Given a qual string, calculates the average phred score
    """
    scores = [phred_dict[x] for x in qual_string]
    return sum(scores)/len(scores)

def phred_score_dict(offset):
    """ Creates a dict of phred score values
    """
    ascii_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    offset_string = ascii_string[offset - 33:]
    
    phred_dict = {}
    for i in range(len(offset_string)):
        phred_dict[offset_string[i]] = float(i)
    
    return phred_dict


if __name__ == '__main__':
    
    args = {}
    
    args['in'] = sys.argv[1]
    args['out'] = sys.argv[2]
    
    main(args)
