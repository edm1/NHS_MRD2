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
        self.j_del = None
        self.v_del = None
        self.seq = None
        self.insert_start = 0
        self.insert_end = -1
        self.is_phix = False
    
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
        self.v_del = int(ref_len - (alignment.aend))
        # 0 based index of the last insert base
        self.insert_start = int(alignment.qend + 1)
        # Add the seq
        if not self.seq:
            self.seq = alignment.seq
    
    def get_header(self):
        return '{0}_{1},{2}:{3},{4}'.format(self.query_name,
                                            self.j_ref,
                                            self.j_del,
                                            self.v_ref,
                                            self.v_del)
    
    def get_insert(self):
        return self.seq[self.insert_start:self.insert_end]

class Group:
    
    def __init__(self, parent_read):
        self.j = parent_read.j_ref
        self.v = parent_read.v_ref
    
    def check_membership(self, new_read):
        # If read shares the same V and J, add to group and return True
        if new_read.j_ref == self.j and new_read.v_ref == self.v:
            return True
        else:
            return False

def main(args):
    
    # Check out_dir exists and make it if not
    for folder in [args['out_dir'], os.path.join(args['out_dir'], 'groups')]:
        if not os.path.exists(folder):
            os.mkdir(folder)
    
    print 'Parsing SAMs...'
    read_dict, ref_names = parse_sams(args['J_in'], args['V_in'],
                                      args['out_dir'])
    

def parse_sams(j_sam, v_sam, out_dir):
    
    # Open sam iterators
    j_handle = pysam.Samfile(j_sam, 'r')
    v_handle = pysam.Samfile(v_sam, 'r')
    j_iter = j_handle.fetch()
    v_iter = v_handle.fetch()
        
    # Get refname dicts
    ref_names = {}
    ref_names['J'] = get_ref_name_dict(j_handle)
    ref_names['V'] = get_ref_name_dict(v_handle)
    
    # File names
    phix_fasta = os.path.join(out_dir, 'phiX174_reads.fasta')
    unmapped_fasta = os.path.join(out_dir, 'unmapped_reads.fasta')
    group_dir = os.path.join(out_dir, 'groups')
    
    # Counters
    total_count = 0
    unmapped_count = 0
    phiX_count = 0
    mapped_count = 0
    
    # Initiate list of groups
    groups = []
    
    # Iterate over J and V sams
    for j_align in j_iter:
        v_align = v_iter.next()
        
        # Get read size
        read_size = int(str(j_align.qname).split('size=')[-1])
        
        total_count += 1 * read_size
        
        # Skip if both are not mapped
        if j_align.is_unmapped and v_align.is_unmapped:
            with open(unmapped_fasta, 'a') as out_handle:
                write_fasta(out_handle, j_align.qname, j_align.seq)
            unmapped_count += 1 * read_size
            continue
        
        # Skip if read is phix
        if not j_align.rname == -1:
            if ref_names['J'][j_align.rname] == 'phiX174':
                with open(phix_fasta, 'a') as out_handle:
                    write_fasta(out_handle, j_align.qname, j_align.seq)
                phiX_count += 1 * read_size
                continue
        
        mapped_count += 1 * read_size
        
        # Get V length of reference
        v_ref_len = v_handle.lengths[v_align.rname]
        
        # Add details to read_record
        read_record = Read(j_align.qname)
        if not j_align.is_unmapped:
            read_record.parse_J_attr(j_align)
        if not v_align.is_unmapped:
            read_record.parse_V_attr(v_align, v_ref_len)
        
        # Find read's group number
        group_num = 0
        make_new = True
        for group in groups:
            if group.check_membership(read_record):
                make_new = False
                break
            group_num += 1
        if make_new:
            groups.append(Group(read_record))
            
        # Write the reads N-D-N seq to groups fasta file
        fasta_name = os.path.join(group_dir, 'group_{0}.fasta'.format(group_num))
        with open(fasta_name, 'a') as out_handle:
            write_fasta(out_handle, read_record.get_header(),
                        read_record.get_insert())
    
    print 'PhiX174 reads: {0}'. format(phiX_count)
    print 'Reads unmapped: {0}'. format(unmapped_count)
    print 'Reads mapped: {0}'. format(mapped_count)
    print 'Total reads: {0}'. format(total_count)

    
    # Close sam iterators
    j_handle.close()
    v_handle.close()
    
    return 1,1

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
