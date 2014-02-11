# python script <in_fasta> <out_fasta>

import os
import sys
import pysam
from libs.bio_file_parsers import write_fasta

class Read:
    
    def __init__(self, query_name):
        self.query_name = str(query_name)
        self.j_ref = None
        self.v_ref = None
        self.seq = None
        self.insert_start = 0
        self.insert_end = -1
    
    def parse_J_attr(self, alignment):
        # Add reference match to list
        self.j_ref = int(alignment.rname)
        # J deletion size
        self.j_del = int(alignment.pos)
        # 0 based index of the first insert base
        self.insert_end = int(alignment.qstart - 1)
        # Add the seq
        if not self.seq:
            self.seq = alignment.seq

    def parse_V_attr(self, alignment, ref_len):
        # Add reference match to list
        self.v_ref = int(alignment.rname)
        # V deletion size
        self.v_del = int(ref_len - (alignment.aend + 1))
        # 0 based index of the last insert base
        self.insert_start = int(alignment.qend + 1)
        # Add the seq
        if not self.seq:
            self.seq = alignment.seq        

class Group:
    
    def __init__(self, parent_read):
        self.j = parent_read.j_ref
        self.v = parent_read.v_ref
        self.members = [parent_read]
    
    def check_membership(self, new_read):
        # If read shares the same V and J, add to group and return True
        if new_read.j_ref == self.j and new_read.v_ref == self.v:
            self.members.append(new_read)
            return True
        else:
            return False

def main(args):
    
    # Check out_dir exists and make it if not
    if not os.path.exists(args['out_dir']):
        os.mkdir(args['out_dir'])
    
    print 'Parsing SAMs...'
    read_dict, ref_names = parse_sams(args['J_in'], args['V_in'])

    print 'Grouping reads...'
    groups = group_reads(read_dict)
    print '-Number of groups: {0}'.format(len(groups))
    
    print 'Writing group N-D-N subsequences...'
    write_group_subsequences(groups, read_dict, args['out_dir'])
    
    #~ print 'Writing subsequence fasta...'
    #~ write_sub_sequences('results/clone_specific_seqs.fasta', args['J_in'],
                        #~ read_dict)

def write_group_subsequences(groups, read_dict, out_dir):
    """ Takes the list of groups and writes the N-D-N sequences to a separate
        fasta for each group.
    """
    group_num = 1
    for group in groups:
        fasta_file = os.path.join(out_dir, 'group_{0}.fasta'.format(group_num))
        write_sub_sequences(fasta_file, read_dict, group.members)
        group_num += 1

def write_sub_sequences(out_fasta, read_dict, read_list):
    """ Writes a fasta file containing the subsequences between the J and V
        regions.
    """
    read_name_list = [read.query_name for read in read_list]
    
    with open(out_fasta, 'w') as out_handle:
        for read_name in read_name_list:
            # Get read record
            record = read_dict[read_name]
            # Write fasta
            seq = record.seq[record.insert_start:record.insert_end]
            print len(seq)
            write_fasta(out_handle, read_name, seq)

def group_reads(read_dict):
    """ Takes the dictionary of reads and groups them by V/J usage
    """
    i = 0
    
    # Group reads by J and V
    groups = []
    for read_name in sorted(read_dict):
        read = read_dict[read_name]
        make_new = True
        for group in groups:
            if group.check_membership(read):
                make_new = False
                break
        if make_new:
            groups.append(Group(read))
        i += 1
    
    return groups

def parse_sams(j_sam, v_sam):
    
    # Input files with segment as key
    in_files = {'J': j_sam,
                'V': v_sam}
    
    # Dictionary to hold all read records
    read_dict = {}
    ref_names = {}
    
    for segment in in_files:
    #~ for segment in ['V']:
        with pysam.Samfile(in_files[segment], 'r') as sam_file:
            
            # Get names for the ref numbers
            ref_names[segment] = get_ref_name_dict(sam_file)
            
            # Go through each read in the alignment
            for alignedread in sam_file.fetch():
                
                # Skip if read is unmapped
                if alignedread.is_unmapped:
                    continue
                
                # Reads should be in reverse orientation (J -> V)
                if not alignedread.is_reverse:
                    print 'Warning! Mapping incorrect orientation: {0}'.format(alignedread.qname)
                    continue
                
                #~ print 'List of pairs:\n{0}'.format(alignedread.aligned_pairs)
                #~ print 'Mapping quality: {0}'.format(alignedread.mapq)
                #~ print 'Ref start: {0}'.format(alignedread.pos)
                #~ print 'Ref end + 1: {0}'.format(alignedread.aend)
                #~ print 'Query name: {0}'.format(alignedread.qname)
                #~ print 'Query start: {0}'.format(alignedread.qstart)
                #~ print 'Query end from pairs: {0}'.format(alignedread.aligned_pairs[-1][0])
                #~ print 'Query end: {0}'.format(alignedread.qend)
                #~ print 'Query length: {0}'.format(alignedread.qlen)
                #~ print 'Read length: {0}'.format(alignedread.rlen)
                #~ print 'Inferred length: {0}'.format(alignedread.inferred_length)
                #~ print 'T length: {0}'.format(alignedread.tlen)
                #~ print 'Ref number: {0}'.format(alignedread.rname)
                #~ print 'Ref name: {0}'.format(sam_file.getrname(alignedread.rname))
                #~ print 'CIGAR str: {0}'.format(alignedread.cigarstring)
                #~ print 'Ref lenght: {0}'.format(sam_file.lengths[alignedread.rname])
                #~ print
                #~ sys.exit()
                
                # Get length of reference
                ref_len = sam_file.lengths[alignedread.rname]
                
                # Make an instance for the read and link to it with record
                if not alignedread.qname in read_dict:
                    read_dict[alignedread.qname] = Read(alignedread.qname)
                record = read_dict[alignedread.qname]
                
                # Parse alignment info
                if segment == 'J':
                    record.parse_J_attr(alignedread)
                elif segment == 'V':
                    record.parse_V_attr(alignedread, ref_len)
                
    return read_dict, ref_names

def get_ref_name_dict(sam_file):
    """ Given a sam file stream, it will return a dictionary of the reference
        name for each reference number.
    """
    name_dict = {}
    name_dict[None] = None
    
    for i in range(0, sam_file.nreferences):
        name_dict[i] = sam_file.getrname(i)
    
    return name_dict

if __name__ == '__main__':
    
    args = {}
    
    args['J_in'] = sys.argv[1]
    args['V_in'] = sys.argv[2]
    args['out_dir'] = sys.argv[3]
    
    main(args)
