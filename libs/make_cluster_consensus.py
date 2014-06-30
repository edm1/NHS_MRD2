# pypy script <in_fasta> <out_fasta>

import re
import sys
from bio_file_parsers import fastq_parser, clstr_parser, phred_score_dict
from random import choice
from probabilistic_seq_match import sequences_match_prob
from operator import itemgetter

def main(args):
    make_consensus(args['in_fastq'], args['in_clstr'], args['out_file'])
    

def make_consensus(ndn_fastq, clstr_meta, out_fastq):
    """ Main function for creating consensus sequences.
    """
    # Load fasta into a dictionary
    fastq_dict = make_fastq_dict(ndn_fastq)
    num_of_clusters = 0
    total_clus_size = 0
    
    # Pattern to pick out header
    pattern = re.compile(r">(.*)\.\.\..*((([0-9]+):([0-9]+):([0-9]+):([0-9]+))|\*)")
    
    # Phred score dict
    phred_dict = phred_score_dict(33)
    phred_dict_inv = {int(v):k for k, v in phred_dict.items()}
    
    out_clusters = []
    
    with open(clstr_meta, 'r') as in_handle:
        for _, lines in clstr_parser(in_handle):
                
                # Check members of each cluster are a good match
                separate_clusters = check_cluster_membership(lines, fastq_dict, pattern, phred_dict)
                
                # For each separate cluster, create a new consensus seq and qual string
                for cluster_list in separate_clusters:
                    num_of_clusters += 1
                    
                    # Append to a list
                    header, cons_seq, cons_qual, clus_size = clusters_consensus(cluster_list, phred_dict, phred_dict_inv)
                    total_clus_size += clus_size
                    out_clusters.append((header, cons_seq, cons_qual, clus_size))
    
    # Write fastq entries in size order
    with open(out_fastq, 'w') as out_handle:
        for header, cons_seq, cons_qual, clus_size in sorted(out_clusters, reverse=True, key=itemgetter(3)):
            out = ['@' + header, cons_seq, '+', cons_qual]
            out_handle.write('\n'.join(out) + '\n')
    
    return num_of_clusters, total_clus_size

def clusters_consensus(cluster_list, phred_dict, phred_dict_inv):
    """ Creates a consensus sequence by taking the highest quality base
        in each position.
    """
    # Unzip
    headers, seqs, quals = zip(*cluster_list)
    seq_len = len(seqs[0])
    
    # Convert qual string to phred scores
    quals = list(quals)
    for i in range(len(quals)):
        quals[i] = [phred_dict[x] for x in quals[i]]
    
    # Create consensus seq and qual
    cons_seq = [None] * seq_len
    cons_qual = [None] * seq_len
    for i in range(seq_len):
        # Check that more that 50% of reads have a base here and that at least
        # 1 base has quailty >= 20
        if not check_50perc_bases(seqs, quals, i, 20):
            cons_seq[i] = '-'
            cons_qual[i] = '-'
        else:
            # Find the highest quality base
            d = {'A':0, 'G':0, 'C':0, 'T':0, 'N':0}
            for j in range(len(seqs)):
                base = seqs[j][i]
                qual = int(quals[j][i])
                try:
                    d[base] = max(d[base], qual)
                except KeyError:
                    # Catches '-'
                    pass
            best_base, best_qual = sorted(d.iteritems(), key=itemgetter(1), reverse=True)[0]
            # Build consensus seq
            cons_seq[i] = best_base
            cons_qual[i] = best_qual
    
    # Convert to strings
    cons_seq = ''.join(cons_seq)
    cons_qual = ''.join([phred_dict_inv[x] if not x == '-' else '-' for x  in cons_qual  ])
    
    # Calc new cluster size
    clus_size = sum_cluster_sizes(headers)
    # Make new header
    cent_header = headers[0].split(';size=')[0] + ';size={0}'.format(clus_size)
    
    # Strip padding
    cons_seq_unpadded, cons_qual_unpadded = strip_padding(cons_seq, cons_qual)
    
    # Split into contig
    cons_seq_contig, cons_qual_contig = split_into_contiguous(cons_seq_unpadded, cons_qual_unpadded)
    
    return cent_header, cons_seq_contig, cons_qual_contig, clus_size

def split_into_contiguous(seq, qual):
    """ Splits seq and qual by '-' then returns the longest
    """
    # Split seq and qual into contigs
    seq_qual = zip(list(seq), list(qual))
    contigs = []
    contig = []
    for i in range(len(seq_qual)):
        if seq_qual[i][0] != '-':
            contig.append(seq_qual[i])
        else:
            contigs.append(contig)
            contig = []
    contigs.append(contig)
    
    # Take the longest contigs
    longest = sorted(contigs, key=len, reverse=True)[0]
    seq, qual = zip(*longest)
    seq = ''.join(seq)
    qual = ''.join(qual)

    return seq, qual    

def strip_padding(seq, qual):

    # Strip front
    i = 0
    while seq[i] == '-':
        i += 1
    seq = seq[i:]
    qual = qual[i:]
    
    # Strip back
    i = len(seq) - 1
    while seq[i] == '-':
        i = i - 1
    seq = seq[:i+1]
    qual = qual[:i+1]
    
    return seq, qual
    
def sum_cluster_sizes(header_list):
    summed = 0
    for header in header_list:
        size = int(header.split(';size=')[1])
        summed += size
    return summed

def check_50perc_bases(seq_list, qual_list, i, qual_thres):
    """ Returns True if 50% of the bases are not '-'
    """
    num_present = 0
    for seq in seq_list:
        if seq[i] != '-':
            num_present += 1
    
    #~ max_qual = max([x[i] for x in qual_list])
    if float(num_present) / len(seq_list) >= 0.5: # and max_qual >= qual_thres:
        return True
    else:
        return False

def check_cluster_membership(members, fastq_dict, pattern, phred_dict):
    """ Checks each cd-hit cluster to see if members are a good match.
        Returns list of separate clusters if not.
    """
    
    # Return if there is only 1 member
    if len(members) == 1:
        header = pattern.search(members[0]).group(1)
        seq, qual = fastq_dict[header]
        return [[(header, seq, qual)]]
    
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
    max_front, max_back = 0, 0 # Max number of insertions at front and back
    for other in others:
        other_seq = fastq_dict[other][0]
        front, back = align_seq(centroid_seq, other_seq)
        alignments[other] = (front, back)
        # Update max insertions
        if front > max_front:
            max_front = front
        if back > max_back:
            max_back = back
    
    # Represent alignment seqs and qual strings
    aligned = []
    # Add centroid at top
    seq = pad_alignment(fastq_dict[centroid][0], max_front, max_back, 0, 0)
    qual = pad_alignment(fastq_dict[centroid][1], max_front, max_back, 0, 0)
    aligned.append((centroid, seq, qual))
    # Add each of the others
    for other in others:
        seq = pad_alignment(fastq_dict[other][0], max_front, max_back, alignments[other][0], alignments[other][1])
        qual = pad_alignment(fastq_dict[other][1], max_front, max_back, alignments[other][0], alignments[other][1])
        aligned.append((other, seq, qual))
    
    # Re-cluster members using their qual strings for information
    separate_clusters = recluster_members(aligned, phred_dict)
    
    return separate_clusters

def recluster_members(aligned_list, phred_dict):
    
    # New clusters
    separate_clusters = []
    
    # Do re-clustering
    for query in aligned_list:
        new_centroid = True
        q_header, q_seq, q_qual = query
        # For each of the new cluster groups
        for group in separate_clusters:
            added = False
            # For each member in the group
            for target in group:
                # Take seq and qual from the first member in group
                t_header, t_seq, t_qual = target
                # Calculate the probability of a match
                prob = sequences_match_prob(q_seq, q_qual, t_seq, t_qual, phred_dict, 0)
                # Arbitarily selected 1e-20 as threshold
                if prob > 1e-20:
                    group.append(query)
                    new_centroid = False
                    added = True
                    break
            if added == True:
                break
        # If it didn't match any then add as new cluster
        if new_centroid == True:
            separate_clusters.append([query])
    
    # For viewing clusters
    #~ for clus in separate_clusters:
        #~ for align in clus:
            #~ print align[1]
        #~ print
    
    return separate_clusters
    

def pad_alignment(seq, max_front, max_back, front, back):
    
    pad = '-'
    add_front = max_front - front
    add_back = max_back - back
    return add_front * pad + seq + add_back * pad
    
def align_seq(centroid, other):
    """ Returns the (front, back) front/back are the numer of insertions or del
        of the other compared to the centroid
    """
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
