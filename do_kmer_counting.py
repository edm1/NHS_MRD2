# pypy script <in_fasta> <out_fasta>

import sys
from libs.bio_file_parsers import fasta_parser

def main(args):
    
    kmers = {}
    read_names = []
    kmer_size = 5
    
    # First pass to get read names
    with open(args['in_fasta'], 'r') as in_handle:
        for header, seq in fasta_parser(in_handle):
            read_names.append(header)

    # Second pass to count kmers
    with open(args['in_fasta'], 'r') as in_handle:
        
        read_idx = 0
        for header, seq in fasta_parser(in_handle):
            # Count kmers using sliding window
            for i in range(len(seq) - kmer_size + 1):
                
                kmer = seq[i:i + kmer_size]
                
                if not kmer in kmers:
                    kmers[kmer] = [0] * len(read_names)
                    
                    kmers[kmer][read_idx] += 1
            
            print read_idx
            read_idx += 1
            if read_idx == 100:
                break
                
    #~ # Second pass to count kmers
    #~ with open(args['in_fasta'], 'r') as in_handle:
        #~ 
        #~ read_idx = 0
        #~ for header, seq in fasta_parser(in_handle):
            #~ # Count kmers using sliding window
            #~ for i in range(len(seq) - kmer_size + 1):
                #~ 
                #~ kmer = seq[i:i + kmer_size]
                #~ 
                #~ if not kmer in kmers:
                    #~ kmers[kmer] = {}
                #~ 
                #~ try:
                    #~ kmers[kmer][read_idx] += 1
                #~ except KeyError:
                    #~ kmers[kmer][read_idx] = 1
            #~ 
            #~ print read_idx
            #~ read_idx += 1
            #~ if read_idx == 1000000:
                #~ break
    
    #~ print kmers
    print len(kmers)
    #~ sys.exit()
    
    # Populate pandas dataframe
    


if __name__ == '__main__':
    
    args = {}
    args['in_fasta'] = sys.argv[1]
    
    main(args)
