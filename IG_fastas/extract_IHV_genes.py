# pypy script <in_fasta> <out_fasta>

import sys
from math import ceil

def main(args):
    
    with open(args['out'], 'w') as out_handle:
        
        for record in fasta_parser(open(args['in'], 'r')):
            if 'IGHV' in record[0]:
                write_fasta(out_handle, record[0], record[1])

def fasta_parser(handle):
    """Generator function to iterator over Fasta records (as string tuples).

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> for values in SimpleFastaParser(open("Fasta/dups.fasta")):
    ...     print(values)
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

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

        #Remove trailing whitespace, and any internal spaces
        #(and any embedded \r which are possible in mangled files
        #when not opened in universal read lines mode)
        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return  # StopIteration

def write_fasta(handle, header, seq, max_line=79):
    """ Will write the fasta sequence to the handle.
    """
    
    # Calculate how many columns per line
    num_lines = int(len(seq) / max_line) + 1
    len_line = int(ceil(float(len(seq))/num_lines))
    
    # Write header
    handle.write('>')
    handle.write(header)
    handle.write('\n')
    
    # Write wrapped lines
    for seq_part in wrap(seq, len_line):
        handle.write(seq_part)
        handle.write('\n')
        
    return 0

def wrap(string, length):
    """ Yield successive length-sized chunks from string.
    """
    for i in xrange(0, len(string), length):
        yield string[i:i + length]

if __name__ == '__main__':
    
    args = {}
    
    args['in'] = sys.argv[1]
    args['out'] = sys.argv[2]
    
    main(args)
