# cython script <in_fasta> <out_fasta>

import sys

def main(bytes targets_in, followup_fastq):
    

    cdef bytes target, title, seq, qual
    cdef char* ascii
    ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    
    # For each target
    for clone_seq, clone_qual in target_parser('temp'):
#~         print target
        # For each record check prob match
        for target_title, target_seq, target_qual in fastq_parser(followup_fastq):
#~                 print len(target_seq)
#~                 print len(clone_seq)
#~                 break
           prob = sequences_match_prob(target_seq, target_qual, clone_seq, clone_qual, ascii, 0.0)
#~            if prob > 0.8:
#~                print target_title
#~                print prob
    
    return 0 

def sequences_match_prob(bytes a_seq, bytes a_qual, bytes b_seq,
                         bytes b_qual, char* ascii, double stop_thresh):
    """ Given two sequences and their quality scores
    """
    
    cdef double prob, prob_A, prob_B, match_prob
    
    # Calc prob
    prob = 1.0
    
    # For each base
    for i in range(len(a_seq)):
        # If Either is N, then prob is 0.25
        if a_seq[i] == 'N' or b_seq == 'N':
            match_prob = 0.25
        else:
            # Calculate the base probabilities
            prob_A = base_prob(phred_score(a_qual[i], 33, ascii))
            prob_B = base_prob(phred_score(b_qual[i], 33, ascii))
            # Calc probability of a match
            if a_seq[i] == b_seq[i]:
                match_prob = match_given_match_prob(prob_A, prob_B)
            elif a_seq[i] != b_seq[i]:
                match_prob = match_given_mismatch_prob(prob_A, prob_B)
        # Combine with overall prob
        prob = prob * match_prob
        #~ cumul.append(prob)
        # Return if less than threshold
        if prob < stop_thresh:
            #~ print cumul
            return prob
    return prob

def match_given_mismatch_prob(double x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (  (1.0/3) * (1 - x_prob) * y_prob
            + (1.0/3) * (1 - y_prob) * x_prob
            + (2.0/9) * x_prob * y_prob )

def match_given_match_prob(double x_prob, y_prob):
    """ Gives the prob of true match given two bases match.
    """
    return (1 - x_prob) * (1 - y_prob) + (x_prob * y_prob) / 3

def base_prob(int phred_score):
    """ Returns the probabilty that a base is incorrect, given its
        Phred score.
    """
    cdef double prob
    prob = 10.0**((-<double>phred_score)/10)
    return prob

def phred_score(bytes letter, int offset, char* ascii):
    """ Returns the Phred score of a fastq quality character.
    """
    cdef int score
    
    # Offset the string
    offset_string = ascii[offset - 33:]
    # Find the characters score
    score = 0
    while score < len(offset_string):
        if letter == offset_string[score]:
            return score
        score += 1
    # If no score is found then there must be an error
    raise ValueError, 'Invalid fastq quality character'

def target_parser(bytes target_file):
    """ Temp placeholder
    """
    seqs = [('AGGGTTCCCTGGCCCCAGGGGTCGAACCAGGAGATAGCAGCTGGTACTACTACAATAAAGGCCCCCCGTCTCTCGCACAGTAATACACGGCCGTGTCCTCAGATCTCAGGCTGCTCAGCTCCATGTAGACTGTGCTCGTGGACGTGTCCCTGGTCATGGTGACTCTGCCCTGGAACTTCTGTGCGTAGCTTGTGCTACCACCACTAGGGTTGATTATTCCCATCCACTCAAGCCCTTGTCCAGGGGCCTGTCGCACCCAGTGCATATAGTAGCTGGTGAAGGTGTATCCAGATGCCTTGC',
             '11>AAFFFFFFAFGFE?ECCCGGGCFHGCFGGGHHHHHGHHHHGHHHHHHHHHHHHHGHGHHHGGGGGGGHHHHGGGGGHHHHHHHGHGGGGGGGHHHHHHGHHHHHHHHHHHGHHHHHHHHHHHHHGHHHHHHHHHGGHGGHGGGHHHHHHHHHGHHHHHGHHHHHHGHHHGHEHHHHHGHHHGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGFFFFFFFFFFFFFFFFFEFFFFFFFFFFB@@@FFFFFE@@@FF?FBFBFFFF/FFFFFFBFFFFFBFFFFFFFFEBBFEFFFF9')]
    for seq in seqs:
        yield seq

def fastq_parser(bytes fastq_file):
    """Function from BioPython that parses a fastq
    """
    cdef bytes line, title_line, seq_string, quality_string
    
    with open(fastq_file, 'r') as handle:
        
        #We need to call handle.readline() at least four times per record,
        #so we'll save a property look up each time:
        handle_readline = handle.readline
    
        #Skip any text before the first record (e.g. blank lines, comments?)
        while True:
            line = handle_readline()
            if not line:
                return  # Premature end of file, or just empty?
            if line[0] == "@":
                break
            if isinstance(line[0], int):
                raise ValueError("Is this handle in binary mode not text mode?")
    
        while line:
            if line[0] != "@":
                raise ValueError(
                    "Records in Fastq files should start with '@' character")
            title_line = line[1:].rstrip()
            #Will now be at least one line of quality data - in most FASTQ files
            #just one line! We therefore use string concatenation (if needed)
            #rather using than the "".join(...) trick just in case it is multiline:
            seq_string = handle_readline().rstrip()
            #There may now be more sequence lines, or the "+" quality marker line:
            while True:
                line = handle_readline()
                if not line:
                    raise ValueError("End of file without quality information.")
                if line[0] == "+":
                    #The title here is optional, but if present must match!
                    second_title = line[1:].rstrip()
                    if second_title and second_title != title_line:
                        raise ValueError("Sequence and quality captions differ.")
                    break
                seq_string += line.rstrip()  # removes trailing newlines
            #This is going to slow things down a little, but assuming
            #this isn't allowed we should try and catch it here:
            if " " in seq_string or "\t" in seq_string:
                raise ValueError("Whitespace is not allowed in the sequence.")
            seq_len = len(seq_string)
    
            #Will now be at least one line of quality data...
            quality_string = handle_readline().rstrip()
            #There may now be more quality data, or another sequence, or EOF
            while True:
                line = handle_readline()
                if not line:
                    break  # end of file
                if line[0] == "@":
                    #This COULD be the start of a new sequence. However, it MAY just
                    #be a line of quality data which starts with a "@" character.  We
                    #should be able to check this by looking at the sequence length
                    #and the amount of quality data found so far.
                    if len(quality_string) >= seq_len:
                        #We expect it to be equal if this is the start of a new record.
                        #If the quality data is longer, we'll raise an error below.
                        break
                    #Continue - its just some (more) quality data.
                quality_string += line.rstrip()
    
            if seq_len != len(quality_string):
                raise ValueError("Lengths of sequence and quality values differs "
                                 " for %s (%i and %i)."
                                 % (title_line, seq_len, len(quality_string)))
    
            #Return the record and then continue...
            yield (title_line, seq_string, quality_string)
        raise StopIteration
#~     
