#!/usr/bin/env python
# -*- coding: utf-8 -*-
#


import sys
import pylab as pl
import os

def main(args):
    
    # Get lengths
    lengths = get_lengths(args['in_ndn_derep'])
    
    # Get plot parameters
    title = os.path.split(os.path.split(args['in_ndn_derep'])[0])[-1]
    plrange = max(lengths) - min(lengths)
    
    # Make histogram
    pl.figure()
    pl.hist(lengths, bins=(plrange/2))
    pl.title(title)
    pl.show()
    
    return 0

def get_lengths(in_file):
    lengths = []
    with open(in_file, 'r') as in_handle:
        for header, seq in fasta_parser(in_handle):
            length = len(seq)
            size = int(header.split(';size=')[1])
            for i in range(size):
                lengths.append(length)
    return lengths

def fasta_parser(handle):
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

if __name__ == '__main__':
    
    args = {}
    args['in_ndn_derep'] = sys.argv[1]
    
    main(args)
