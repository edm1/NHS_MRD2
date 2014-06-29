# pypy script <in_fasta> <out_fasta>

import re
import sys
from bio_file_parsers import write_fasta, fastq_parser, clstr_parser
from random import choice

def main(args):
    make_consensus(args['in_fastq'], args['in_clstr'], args['out_file'])
    

def make_consensus(ndn_fastq, clstr_meta, out_fasta):
    """ Main function for creating consensus sequences.
    """
    # Load fasta into a dictionary
    fastq_dict = make_fastq_dict(ndn_fastq)
    num_of_clusters = 0
    
    # Pattern to pick out header
    pattern = re.compile(r">(.*)\.\.\..*((([0-9]+):([0-9]+):([0-9]+):([0-9]+))|\*)")
    
    # Check members of each cluster are a good match
    with open(clstr_meta, 'r') as in_handle:
        for _, lines in clstr_parser(in_handle):
            split_clusters = check_cluster_membership(lines, fastq_dict, pattern)
    
    # Make a consensus for each cluster
    #~ with open(clstr_meta, 'r') as in_handle:
        #~ with open(out_fasta, 'w') as out_handle:
            #~ for title, lines in clstr_parser(in_handle):
                #~ num_of_clusters += 1
                #~ if len(lines) > 1:
                    #~ print header # DEBUG
                    #~ header, cons_seq = make_consensus_seq(lines, fasta_dict, pattern)
                #~ else:
                    #~ header = pattern.search(lines[0]).group(1)
                    #~ cons_seq = fasta_dict[header]
                #~ write_fasta(out_handle, header, cons_seq)
    #~ return num_of_clusters

def check_cluster_membership(members, fastq_dict, pattern):
    """ Checks each cd-hit cluster to see if members are a good match.
        Returns list of separate clusters if not.
    """
    
    # Return if there is only 1 member
    if len(members) == 1:
        header = pattern.search(members[0]).group(1)
        seq, qual = fastq_dict[header]
        return [(header, seq, qual)]
    
    # Find centroid and others
    centroid = None
    others = []
    for member in members:
        match_obj = pattern.search(member)
        header = match_obj.group(1)
        if match_obj.group(2) == "*":
            centroid = header
        else:
            others.append(header)
    
    # Find optimal alignment for each other against centroid
    alignments = {}
    centroid_seq = fastq_dict[centroid][0]
    max_front, max_back = 0, 0
    for other in others:
        other_seq = fastq_dict[other][0]
        front, back = align_seq(centroid_seq, other_seq)
        alignments[other] = (front, back)
        if front > max_front:
            max_front = front
        if back > max_back:
            max_back = back
    
    # Represent alignment
    aligned = []
    # Add centroid at top
    aligned.append(max_front * '-' + centroid_seq + max_back * '-')
    # Add each of the others
    for other in others:
        other_seq = fastq_dict[other][0]
        add_front = max_front - alignments[other][0]
        add_back = max_back - alignments[other][1]
        aligned.append(add_front * '-' + other_seq + add_back * '-')
    
    for align in aligned:
        print align
    print
    
    
    #~ other_seq = fastq_dict[others[0]][0]
    #~ centroid_seq =  'ACATINAHATDOGDOG'
    #~ other_seq = 'FACEACATINAHATLYLKY'
    #~ front, back = align_seq(centroid_seq, other_seq)
    #~ 
    #~ if front > 0 :
        #~ centroid_seq = '-' * front + centroid_seq
    #~ elif front < 0:
        #~ other_seq = '-' * -front + other_seq
    #~ 
    #~ if back > 0 :
        #~ centroid_seq = centroid_seq + '-' * back
    #~ elif back < 0:
        #~ other_seq = other_seq + '-' * -back
         #~ 
    #~ print centroid_seq
    #~ print other_seq
    
    
    #~ sys.exit()

def align_seq(centroid, other):
    
    cen_len = len(centroid)
    other_len = len(other)
    
    # Slide each seq against each other and measure number of matches
    best_score = 0
    best_index = None
    for i in range(1, cen_len + other_len):
        cen_sub = centroid[max(0, i - other_len):min(i, cen_len)]
        other_sub = other[max(0, other_len - i):min(other_len, cen_len + other_len - i)]
        # Skip if not at least 50% of the centroid length
        if not len(cen_sub) > cen_len / 2:
            continue
        score = count_matches(cen_sub, other_sub)
        if score > best_score:
            best_score = score
            best_index = i
    
    # Convert into shorthand, ins/del at front and back of other
    front =  other_len - best_index 
    back = best_index - cen_len
    return (front, back)
        
def count_matches(seqA, seqB):
    score = 0
    for i in range(len(seqA)):
        if seqA[i] == seqB[i] and not seqA[i] == 'N':
            score += 1
    return score

def yield_indexes(other_len, cen_len):
    """ Yields the indexes to slide along 
    """
    pass
    #~ for i in range(other_len):
        #~ yield  (other_len -1 - i, None)
    #~ 
    #~ for i in range(1, other_len):
        #~ yield (None, -i)

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
            # Cluster centroid seq
            title = header
            #~ align.append(fasta_dict[header])
            align.append(fasta_dict[header] + ' *')
        else:
            # Add spaces to the beginning
            seq = "N" * (int(match_obj.group(6)) - 1)
            # Remove bases from the start
            seq += fasta_dict[header][(int(match_obj.group(4)) - 1):]
            align.append(seq)
    
    # Debug
    #~ for al in align:
        #~ print al
    #~ print
    
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
    

def make_fastq_dict(fasta):
    """ Parses the fasta file into a dictionary
    """
    fastq_dict = {}
    with open(fasta, 'r') as in_handle:
        for header, seq, qual in fastq_parser(in_handle):
            fastq_dict[header] = (seq, qual)
    return fastq_dict

if __name__ == '__main__':
    
    args = {}
    args['in_fastq'] = sys.argv[1]
    args['in_clstr'] = sys.argv[2]
    args['out_file'] = sys.argv[3]
    
    main(args)
