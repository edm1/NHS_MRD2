# pypy script <in_fasta> <out_fasta>

import pyximport; pyximport.install()
import probabilistic_detection
import sys

def main():
    
    # Parse args
    targets_in = sys.argv[1]
    followup_reads = sys.argv[2]
    
    # Run cython script
    probabilistic_detection.main(targets_in, followup_reads)
    
    return 0

if __name__ == '__main__':
    
    main()
