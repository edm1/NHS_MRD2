#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Contains simple parsers for fasta and fastq files, taken directly from the
# biopython source code. Also, includes a homemade BLAST xml output parser that
# is written specifically with NHS MRD needs in mind.
#

from xml.etree.ElementTree import ElementTree
from math import ceil # To round up line wrap

def wrap(string, length):
    """ Yield successive length-sized chunks from string.
    """
    for i in xrange(0, len(string), length):
        yield string[i:i + length]

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
    if len(seq) > 0:
        for seq_part in wrap(seq, len_line):
            handle.write(seq_part)
            handle.write('\n')
    else:
        handle.write(seq)
        handle.write('\n')
        
    return 0

def blast_xml_parser(xml_handle):
    """ Will parse the BLAST xml output file and for each query yield a dict
        containing the following values:
        
        dict = {'query_id': str
                'hits': [hit1, hit2, hit3, ...]
                          |
                         hit1 = {'hit_id': str
                                 'hit_len': int
                                 'e'
                                 'qstart': int
                                 'qend': int
                                 'sstart': int
                                 'send': int
                                 'align_len': int
                                 'identity': float
    """
    
    # Parse XML file
    tree = ElementTree()
    tree.parse(xml_handle)
    
    # Split into queries
    iteration = tree.iter('Iteration')
    for query in iteration:
        
        # Initiate dict
        query_dict = {}
        query_dict['query_id'] = query.findtext('Iteration_query-def')
        query_dict['hits'] = []
        
        # Iterate through each hit
        for hit in query.iter('Hit'):
            hit_dict = {}
            hit_dict['hit_id'] = hit.findtext('Hit_id')
            hit_dict['hit_len'] = int(hit.findtext('Hit_len'))
            
            # Move to top Hsp
            top_hsp = hit.iter('Hsp').next()
            
            # Extract everything needed
            hit_dict['e'] = float(top_hsp.findtext('Hsp_evalue'))
            hit_dict['qstart'] = int(top_hsp.findtext('Hsp_query-from'))
            hit_dict['qend'] = int(top_hsp.findtext('Hsp_query-to'))
            hit_dict['sstart'] = int(top_hsp.findtext('Hsp_hit-from'))
            hit_dict['send'] = int(top_hsp.findtext('Hsp_hit-to'))
            hit_dict['align_len'] = int(top_hsp.findtext('Hsp_align-len'))
            hit_dict['identity'] = int(top_hsp.findtext('Hsp_identity'))
            
            # Add this hit to the query hit list
            query_dict['hits'].append(hit_dict)
        
        yield query_dict

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

def fastq_parser(handle):
    """Iterate over Fastq records as string tuples (not as SeqRecord objects).

    This code does not try to interpret the quality string numerically.  It
    just returns tuples of the title, sequence and quality as strings.  For
    the sequence and quality, any whitespace (such as new lines) is removed.

    Our SeqRecord based FASTQ iterators call this function internally, and then
    turn the strings into a SeqRecord objects, mapping the quality string into
    a list of numerical scores.  If you want to do a custom quality mapping,
    then you might consider calling this function directly.

    For parsing FASTQ files, the title string from the "@" line at the start
    of each record can optionally be omitted on the "+" lines.  If it is
    repeated, it must be identical.

    The sequence string and the quality string can optionally be split over
    multiple lines, although several sources discourage this.  In comparison,
    for the FASTA file format line breaks between 60 and 80 characters are
    the norm.

    WARNING - Because the "@" character can appear in the quality string,
    this can cause problems as this is also the marker for the start of
    a new sequence.  In fact, the "+" sign can also appear as well.  Some
    sources recommended having no line breaks in the  quality to avoid this,
    but even that is not enough, consider this example::

        @071113_EAS56_0053:1:1:998:236
        TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA
        +071113_EAS56_0053:1:1:998:236
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
        @071113_EAS56_0053:1:1:182:712
        ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG
        +
        @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
        @071113_EAS56_0053:1:1:153:10
        TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT
        +
        IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
        @071113_EAS56_0053:1:3:990:501
        TGGGAGGTTTTATGTGGA
        AAGCAGCAATGTACAAGA
        +
        IIIIIII.IIIIII1@44
        @-7.%<&+/$/%4(++(%

    This is four PHRED encoded FASTQ entries originally from an NCBI source
    (given the read length of 36, these are probably Solexa Illumna reads where
    the quality has been mapped onto the PHRED values).

    This example has been edited to illustrate some of the nasty things allowed
    in the FASTQ format.  Firstly, on the "+" lines most but not all of the
    (redundant) identifiers are omitted.  In real files it is likely that all or
    none of these extra identifiers will be present.

    Secondly, while the first three sequences have been shown without line
    breaks, the last has been split over multiple lines.  In real files any line
    breaks are likely to be consistent.

    Thirdly, some of the quality string lines start with an "@" character.  For
    the second record this is unavoidable.  However for the fourth sequence this
    only happens because its quality string is split over two lines.  A naive
    parser could wrongly treat any line starting with an "@" as the beginning of
    a new sequence!  This code copes with this possible ambiguity by keeping
    track of the length of the sequence which gives the expected length of the
    quality string.

    Using this tricky example file as input, this short bit of code demonstrates
    what this parsing function would return:

    >>> with open("Quality/tricky.fastq", "rU") as handle:
    ...     for (title, sequence, quality) in FastqGeneralIterator(handle):
    ...         print(title)
    ...         print("%s %s" % (sequence, quality))
    ... 
    071113_EAS56_0053:1:1:998:236
    TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
    071113_EAS56_0053:1:1:182:712
    ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
    071113_EAS56_0053:1:1:153:10
    TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
    071113_EAS56_0053:1:3:990:501
    TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%

    Finally we note that some sources state that the quality string should
    start with "!" (which using the PHRED mapping means the first letter always
    has a quality score of zero).  This rather restrictive rule is not widely
    observed, so is therefore ignored here.  One plus point about this "!" rule
    is that (provided there are no line breaks in the quality sequence) it
    would prevent the above problem with the "@" character.
    """
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

def main():
    """ Used to test the xml parser.
    """
    
    
    import sys
    
    test_xml = sys.argv[1]
    
    for d in blast_xml_parser(open(test_xml, 'r')):
        if d['query_id'] == 'M01996:13:000000000-A5U9G:1:1101:11910:4467':
            print d['hits'][0]['e']
    
    return 0

if __name__ == '__main__':
	main()

 
