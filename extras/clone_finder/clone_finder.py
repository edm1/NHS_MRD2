#!/usr/bin/env python
# -*- coding: utf-8 -*-
#


def main():
    
    # Constants
    global ascii_string
    ascii_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    
    # Try for whole fastq
    from bio_file_parsers import fastq_parser
    import gzip
    
    fastq_iter = fastq_parser(gzip.open('../../test_data/272-Diagnostic-mutliplexFR1_S25_L001_R1_001.fastq.gz', 'r'))
    a = fastq_iter.next()
    a = fastq_iter.next()
    a = fastq_iter.next()
    c = 1
    thresh = 0.01
    for b in fastq_iter:
        prob = sequences_match_prob(a, b, thresh)
        if prob > thresh:
            print c, prob
        c += 1
        #~ print c
    
    # Test data
    #~ seq1 = ('AGGGTTCCCTGGCCCCAGGGGTCGAACCAGGAGATAGCAGCTGGTACTACTACAATAAAGGCCCCCCGTCTCTCGCACAGTAATACACGGCCGTGTCCTCAGATCTCAGGCTGCTCAGCTCCATGTAGACTGTGCTCGTGGACGTGTCCCTGGTCATGGTGACTCTGCCCTGGAACTTCTGTGCGTAGCTTGTGCTACCACCACTAGGGTTGATTATTCCCATCCACTCAAGCCCTTGTCCAG',
            #~ '11>AAFFFFFFAFGFE?ECCCGGGCFHGCFGGGHHHHHGHHHHGHHHHHHHHHHHHHGHGHHHGGGGGGGHHHHGGGGGHHHHHHHGHGGGGGGGHHHHHHGHHHHHHHHHHHGHHHHHHHHHHHHHGHHHHHHHHHGGHGGHGGGHHHHHHHHHGHHHHHGHHHHHHGHHHGHEHHHHHGHHHGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGFFFFFFFFFFFFFFFFFEFFFFFFFFF')
    #~ seq2 = ('AGGGTTCCCTGGCCCCAGGGGTCGAACCAGGAGATAGCAGCTGGTACTACTACAATATCCTAGTCCGGGGGTCGCACAGTAATACACAGCCGTGTCTGCGGCGGTCACAGAGCTCAGCTTCAGGGAGAACTGGTTCTTGGACGTGTCTACGGATATGGTGACTCGACTCTTGAGGGACGGGTTGTAGTAGGTGCTCCCACTATAATAGATACTCCCAATCCACTCCAGCCCCTTCCCTGGGGG',
            #~ 'ABCCCFFFFFFCGGGGGGGGGGGGGGHGHHGHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGGGGGGGGGGHHHHHHHHHHHHHGGGHHHHHGGGGGFGAG1>FF/DGCCBGGGHGHH<C.CCCBG0<::C0G0GG.ABFFFFFG?AFFFG9FFFGGGA?EGGG0:FFD.ADFF-ADA:B;BBF9;/:9BFBF./BFFF///BFFFFFEFFFFFEBF/FFFFFFFEFFFFFF.AAF')
    #~ seq2sub = ('AACCAGGAGATAGCAGCTGGTACTACTAC',
               #~ 'GGHGHHGHGHHHHHHHHHHHHHHHHHHHH')
               #~ 
    #~ sub = (0, 60)
    #~ seq3 = (seq1[0][sub[0]:sub[1]], seq1[1][sub[0]:sub[1]])
    #~ seq4 = (seq2[0][sub[0]:sub[1]], seq2[1][sub[0]:sub[1]])
    #~ 
    #~ print seq3[0]
    #~ print seq4[0]
    #~ print sequences_match_prob(seq3, seq4)
    
    return 0

def optimal_alignment(A, B, min_overlap=5):
    
    # Find longest and shortest
    if len(seqA[0]) >= len(seqB[0]):
        seqLong, qualLong = A
        seqShort, qualShort = B
    else:
        seqLong, qualLong = B
        seqShort, qualShort = A
        
    # Return if shortest is less than min
    if len(seqShort) < min_overlap:
        return 0, None
    
    #

def sequences_match_prob(A, B, stop_thresh=None):
    """ Given two sequences and their quality scores
    """
    
    # Unpack values
    headerA, seqA, qualA = A
    headerB, seqB, qualB = B
    
    # Calc prob
    prob = 1.0
    
    # For each base
    for i in range(len(seqA)):
        # If Either is N, then prob is 0.25
        if seqA[i] == 'N' or seqB[i] == 'N':
            match_prob = 0.25
        else:
            # Calculate the base probabilities
            prob_A = base_prob(phred_score(qualA[i]))
            prob_B = base_prob(phred_score(qualB[i]))
            # Calc probability of a match
            if seqA[i] == seqB[i]:
                match_prob = match_given_match_prob(prob_A, prob_B)
            elif seqA[i] != seqB[i]:
                match_prob = match_given_mismatch_prob(prob_A, prob_B)
        # Combine with overall prob
        prob = prob * match_prob
        #~ cumul.append(prob)
        # Return if less than threshold
        if stop_thresh:
            if prob < stop_thresh:
                #~ print cumul
                return prob
                
    
    return prob
    

def match_given_mismatch_prob(x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (  (1.0/3) * (1 - x_prob) * y_prob
            + (1.0/3) * (1 - y_prob) * x_prob
            + (2.0/9) * x_prob * y_prob )
    
def match_given_match_prob(x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (1 - x_prob) * (1 - y_prob) + (x_prob * y_prob) / 3

def base_prob(phred_score):
    """ Returns the probabilty that a base is incorrect, given its
        Phred score.
    """
    return 10.0**(-float(phred_score)/10)

def phred_score(char, offset=33):
    """ Returns the Phred score of a fastq quality character.
    """
    
    # Offset the string
    offset_string = ascii_string[offset - 33:]
    # Find the characters score
    score = 0
    for i in range(len(offset_string)):
        if char == offset_string[i]:
            return score
        score += 1
    # If no score is found then there must be an error
    raise ValueError, 'Invalid fastq quality character'

if __name__ == '__main__':
    
    args = {}
    main()
