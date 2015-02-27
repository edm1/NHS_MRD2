#!/usr/bin/env pypy
# -*- coding: utf-8 -*-

"""NHS MRD.

Usage:
  nhs_mrd.py detect  <fastq> <id> [--identity <int>] [--top <int>]
                     [--r <regex>] [--germ <regex>] [--ab <int>] [--bb <int>]
                     [--m <int>] [--t <int>] [--e <dec>] [--tree <int>]
                     [--min_dup <int>] [--quiet]
  nhs_mrd.py compare <sample1.dat> <sample2.dat> <id> [--queries <int>]
                     [--matches <int>]
  nhs_mrd.py (-h | --help)
  nhs_mrd.py --version

Options for clone detection:
  --identity <int>  Minimum seq identity to be considered a cluster [default: 98].
  --top <int>      Number of top clusters to report [default: 100].
  --r <regex>      J region conserved nucleotide motif [default: CCCCAG].
  --germ <regex>   Regular expressions that match germ line specific sequences
                   in the J regions. See documentation for defaults.
  --ab <int>       Bases to include after end of motif [default: 125].
  --bb <int>       Bases to include before start of motif [default: 10].
  --m <int>        Max distance from seq start that motif can be [default: 50].
  --t <int>        Number of hits to retrieve in annotation [default: 3].
  --e <dec>        Minimum E-value for BLAST results [default: 0.05].
  --min_dup <int>  The minimum number of reads a duplication group must have.
                   [default: 2].

Options for sample comparison:
  --queries <int>  Number of top clusters to compare between samples
                   [default: 10].
  --matches <int>  Number of matches to report for each query [default: 3].

Other options:
  --help, -h       Show this screen.
  --version        Show version.
  --quiet          Do not print progress

"""

# Import libraries
import os
import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages') # Location of docopt
from docopt import docopt
import platform

def main():
    
    # Check if PyPy is being used and warn user if not
    if not platform.python_implementation() == 'PyPy':
        print ('\nWarning: using pypy will speed things up considerably. It '
               'is recommended that you run the pipeline using the following '
               'command:\n')
        print '\tpypy nhs_mrd.py [options]\n'

    #
    # Get arguments
    #
    
    # Parse arguments from docstring
    args = docopt(__doc__, version='NHS MRD 1.5.3')

    # Need to add default germline regex separately
    if not args['--germ']:
        args['--germ'] = "".join(["CACGGTG[ATGCN]{22,24}AGAAACCCA|",
                                  "CACAGTC[ATGCN]{22,24}ACAAAAACA|",
                                  "CACACAG[ATGCN]{24,26}ACAAAAACC|",
                                  "CACACAG[ATGCN]{24,26}ACACAAACC|",
                                  "CACATTG[ATGCN]{22,24}ACAAAAACC|",
                                  "CACATTG[ATGCN]{20,22}ACAAAGAAC|",
                                  "CACAATG[ATGCN]{21,23}ACAAAAACC"])
    
    #
    # Validate arguments
    #
    
    # Integers
    for int_arg in ['--ab', '--bb', '--m', '--matches',
                    '--queries', '--t', '--top', '--min_dup']:
        if not args[int_arg] == None:
            try:
                args[int_arg] = int(args[int_arg])
            except ValueError:
                raise ValueError, '{0} must be an integer'.format(int_arg)
    
    # Floats
    for float_arg in ['--e', '--identity']:
        if not args[float_arg] == None:
            try:
                args[float_arg] = float(args[float_arg])
            except ValueError:
                raise ValueError, '{0} must be a decimal'.format(float_arg)
    
    # Files
    for file_arg in ['<fastq>', '<sample1.dat>', '<sample2.dat>']:
        if not args[file_arg] == None:
            if not os.path.exists(args[file_arg]):
                raise ValueError, '{0} is not an existing file'.format(file_arg)
    
    # Check identity is between 1 and 100
    err = 'Identity must be between 1-100 (default 98).'
    assert args['--identity'] > 1 and args['--identity'] <= 100, err

    #
    # Run clone detection stage
    #

    if args['detect']:
        
        # Import clone detection
        import stage1
        
        stage1.main(args)
        
    #
    # Run sample comparison stage
    #
    
    elif args['compare']:
        
        # Import sample comparison
        import stage2
        
        # Need to be available as top-level modules so that pickle can import
        # the .dat file
        from libs.stage1_classes import Records
        from libs.stage1_classes import Cluster
        from libs.stage1_classes import Segment
        
        # Run
        stage2.main(args)
        
    return 0

if __name__ == '__main__':
    main()
