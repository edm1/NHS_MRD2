# pypy script <in_fasta> <out_fasta>

import re
import sys
from bio_file_parsers import write_fasta, fasta_parser
from random import choice

def main(args):
    make_consensus(args['in_fasta'], args['in_clstr'], args['out_file'])
    

def make_consensus(ndn_fasta, clstr_meta, out_fasta):
    """ Main function for creating consensus sequences.
    """
    # Load fasta into a dictionary
    fasta_dict = make_fasta_dict(ndn_fasta)
    num_of_clusters = 0
    
    # Make a consensus for each cluster
    with open(clstr_meta, 'r') as in_handle:
        with open(out_fasta, 'w') as out_handle:
            pattern = re.compile(r">(.*)\.\.\..*((([0-9]+):([0-9]+):([0-9]+):([0-9]+))|\*)")
            for title, lines in clstr_parser(in_handle):
                num_of_clusters += 1
                if len(lines) > 1:
                    header, cons_seq = make_consensus_seq(lines, fasta_dict, pattern)
                else:
                    header = pattern.search(lines[0]).group(1)
                    cons_seq = fasta_dict[header]
                write_fasta(out_handle, header, cons_seq)
    return num_of_clusters

def make_consensus_seq(members, fasta_dict, re_pattern):
    """ Takes the cd-hit output clstrs file and produces a consensus seq for
        each cluster.
    """
    
    # Align sequences using clustering meta data
    align = []
    for member in members:
        match_obj = re_pattern.search(member)
        header = match_obj.group(1)
        if match_obj.group(2) == "*":
            title = header
            align.append(fasta_dict[header])
        else:
            # Add spaces to the beginning
            seq = " " * (int(match_obj.group(6)) - 1)
            # Remove bases from the start
            seq += fasta_dict[header][(int(match_obj.group(4)) - 1):]
            align.append(seq)
    
    # Calc consensus
    max_len = max([len(x) for x in align])
    cons_seq = ""
    for i in range(max_len):
        used = []
        for seq in align:
            try:
                used.append(seq[i])
            except IndexError:
                used.append(" ")
        # Add choice to cons seq
        nucl = mode(used)
        cons_seq += nucl
    
    cons_seq = get_longest_contig(cons_seq)
    return title, cons_seq

def get_longest_contig(cons_seq):
    """ Splits the seq by whitespace and returns the longest part.
    """
    parts = cons_seq.split(" ")
    return sorted(parts, key=len, reverse=True)[0]

def mode(nuc_list):
    """ Returns the most common element of a list, picks randomly if its a
        draw.
    """
    # Count nucs in list
    d = {}
    for nuc in nuc_list:
        try:
            d[nuc] += 1
        except KeyError:
            d[nuc] = 1
    # Get nuc choices
    max_count = max([d[x] for x in d])
    choices = []
    for key in d:
        if d[key] == max_count:
            choices.append(key)
    # Return random choice
    if len(choices) == 1:
        return choices[0]
    else:
        return choice(choices)
    

def make_fasta_dict(fasta):
    """ Parses the fasta file into a dictionary
    """
    fasta_dict = {}
    with open(fasta, 'r') as in_handle:
        for header, seq in fasta_parser(in_handle):
            fasta_dict[header] = seq
    return fasta_dict

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
    args['out_file'] = sys.argv[3]
    
    main(args)
