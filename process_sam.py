# pypy script <in_fasta> <out_fasta>

import sys
import pysam
from libs.bio_file_parsers import write_fasta

class Read:
    
    def __init__(self, query_name):
        self.query_name = str(query_name)
        self.j_ref = None
        self.v_ref = None
    
    def parse_J_attr(self, alignment):
        # Add reference match to list
        self.j_ref = int(alignment.rname)
        # J deletion size
        self.j_del = int(alignment.pos)
        # 0 based index of the first insert base
        self.insert_end = int(alignment.qstart - 1)

    def parse_V_attr(self, alignment, ref_len):
        # Add reference match to list
        self.v_ref = int(alignment.rname)
        # V deletion size
        self.v_del = int(ref_len - (alignment.aend + 1))
        # 0 based index of the last insert base
        self.insert_start = int(alignment.qend + 1)

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
    
    print 'Parsing SAMs...'
    read_dict, ref_names = parse_sams(args['J_in'], args['V_in'])
    
    print 'Writing subsequence fasta...'
    write_sub_sequences('results/clone_specific_seqs.fasta', args['J_in'],
                        read_dict)

def write_sub_sequences(out_fasta, j_sam, read_dict):
    """ Writes a fasta file containing the subsequences between the J and V
        regions.
    """
    with open(out_fasta, 'w') as out_handle:
        with pysam.Samfile(j_sam, 'r') as sam_file:
            # Go through each read in the alignment
            for alignedread in sam_file.fetch():
                # Skip if unmapped
                if alignedread.is_unmapped:
                    continue
                # Get read record
                record = read_dict[alignedread.qname]
                #~ try:
                    #~ ins_start = record.insert_start
                #~ except AttributeError:
                    #~ ins_start = 0
                #~ try:
                    #~ ins_end = record.insert_end + 1
                #~ except AttributeError:
                    #~ ins_end = -1
                #~ 
                #~ # Write fasta
                #~ seq = alignedread.seq[ins_start:ins_end]
                #~ write_fasta(out_handle, record.query_name, seq)  
                try:
                    ins_start = record.insert_start
                    ins_end = record.insert_end + 1
                    seq = alignedread.seq[ins_start:ins_end]
                    write_fasta(out_handle, record.query_name, seq)
                except AttributeError:
                    pass

def group_reads():
    """ Not yet implemented
    """
    #~ i = 0
    #~ tot_reads = len(read_dict)
    #~ # Group reads by J and V
    #~ print 'Grouping reads...'
    #~ groups = []
    #~ for read_name in sorted(read_dict):
        #~ read = read_dict[read_name]
        #~ make_new = True
        #~ for group in groups:
            #~ if group.check_membership(read):
                #~ make_new = False
                #~ break
        #~ if make_new:
            #~ groups.append(Group(read))
        #~ 
        #~ print 'Progress: {0}%'.format(int(float(i*100)/tot_reads)) 
        #~ i += 1
    #~ 
    #~ with open('groups.txt', 'w') as out_handle:
        #~ for group in groups:
            #~ try:
                #~ out_handle.write('j_{0} v_{1} size={2}\n'.format(
                     #~ ref_names['J'][group.j], ref_names['V'][group.v], len(group.members)))
            #~ except:
                #~ print group.j, group.v    

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
                
                # Add additional attributes
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
    
    main(args)
