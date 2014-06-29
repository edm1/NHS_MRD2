# python script <in_fasta> <out_fasta>

import os
import sys
import pysam
#~ from libs.bio_file_parsers import write_fasta
from bio_file_parsers import write_fasta, fastq_parser

class Read:
    
    def __init__(self, query_name):
        self.query_name = str(query_name)
        self.j_ref = None
        self.v_ref = None
        self.j_del = None
        self.v_del = None
        self.insert_start = 0
        self.insert_end = -1
        self.is_phix = False
    
    def parse_J_attr(self, alignment):
        # Add reference match to list
        self.j_ref = int(alignment.rname)
        # J deletion size
        self.j_del = int(alignment.pos)
        # 0 based index of the first insert base
        self.insert_end = int(alignment.qstart)

    def parse_V_attr(self, alignment, ref_len):
        # Add reference match to list
        self.v_ref = int(alignment.rname)
        # V deletion size
        self.v_del = int(ref_len - (alignment.aend))
        # 0 based index of the last insert base
        self.insert_start = int(alignment.qend + 1)
    
    def get_header(self):
        parts = str(self.query_name).split(';size=')
        title = parts[0]
        size = int(parts[1])
        
        return '{0}_{1},{2}:{3},{4};size={5}'.format(title,
                                                    self.j_ref,
                                                    self.j_del,
                                                    self.v_ref,
                                                    self.v_del,
                                                    size)

def main(args):
    
    # Check out_dir exists and make it if not
    for folder in [args['out_dir']]:
        if not os.path.exists(folder):
            os.mkdir(folder)
    
    print 'Parsing SAMs...'
    ref_names, metrics, ndn_fasta = parse_sams(args['J_in'], args['V_in'],
                                    args['out_dir'], args['fastq_in'])
    

def parse_sams(j_sam, v_sam, out_dir, fastq_in):
    
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
    ndn_fastq = os.path.join(out_dir, 'NDN_reads.fastq')
    phix_fasta = os.path.join(out_dir, 'phiX174_reads.fasta')
    pUPATrap_fasta = os.path.join(out_dir, 'pUPATrap.fasta')
    unmapped_fasta = os.path.join(out_dir, 'unmapped_reads.fasta')
    
    # Counters
    metrics = {}
    metrics['total_count'] = 0
    metrics['unmapped_count'] = 0
    metrics['phiX_count'] = 0
    metrics['pUPATrap_count'] = 0
    metrics['mapped_count'] = 0
    
    # Dict to hold insert start and ends
    insert_pos = {}
    
    # Iterate over J and V sams
    for j_align in j_iter:
        v_align = v_iter.next()
        
        # Get read size
        read_size = int(str(j_align.qname).split('size=')[-1])
        
        metrics['total_count'] += 1 * read_size
        
        # Skip if both are not mapped
        if j_align.is_unmapped and v_align.is_unmapped:
            with open(unmapped_fasta, 'a') as out_handle:
                write_fasta(out_handle, j_align.qname, j_align.seq)
            metrics['unmapped_count'] += 1 * read_size
            continue
        
        # Skip if read is phix
        if not j_align.rname == -1:
            if ref_names['J'][j_align.rname] == 'phiX174':
                with open(phix_fasta, 'a') as out_handle:
                    write_fasta(out_handle, j_align.qname, j_align.seq)
                metrics['phiX_count'] += 1 * read_size
                continue
            elif ref_names['J'][j_align.rname] == 'pUPATrap-CRV2_vector':
                with open(pUPATrap_fasta, 'a') as out_handle:
                    write_fasta(out_handle, j_align.qname, j_align.seq)
                metrics['pUPATrap_count'] += 1 * read_size
                continue
        
        metrics['mapped_count'] += 1 * read_size
        
        # Get V length of reference
        v_ref_len = v_handle.lengths[v_align.rname]
        
        # Add details to read_record
        read_record = Read(j_align.qname)
        if not j_align.is_unmapped:
            read_record.parse_J_attr(j_align)
        if not v_align.is_unmapped:
            read_record.parse_V_attr(v_align, v_ref_len)
                        
        # Save header, insert start and end to dict
        insert_pos[read_record.query_name] = (read_record.get_header(),
                                              read_record.insert_start,
                                              read_record.insert_end)
    
    # Iterate over the input fastq and output NDN regions
    rev_comp_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    with open(ndn_fastq, 'w') as out_h:
        with open(fastq_in, 'r') as in_h:
            for title, seq, qual in fastq_parser(in_h):
                if title in insert_pos:
                    start = insert_pos[title][1]
                    end = insert_pos[title][2]
                    rev_comp = reverse_complement(seq, rev_comp_dict)[start:end]
                    rev_qual = qual[::-1][start:end]
                    # Seq len must be > 0
                    if not len(rev_comp) > 0:
                        continue
                    # Make output
                    out = ['@{0}'.format(insert_pos[title][0]),
                           rev_comp,
                           '+',
                           rev_qual]
                    out_h.write('\n'.join(out) + '\n')
                
    
    # Print metrics
    for key, value in metrics.iteritems():
        print '{0}: {1}'.format(key, value)
    
    # Close sam iterators
    j_handle.close()
    v_handle.close()
    
    return ref_names, metrics, ndn_fastq

def reverse_complement(seq, rev_comp_dict):
    rev_list = list(seq.upper())
    return ''.join([rev_comp_dict[x] for x in rev_list])[::-1]

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
    args['fastq_in'] = sys.argv[3]
    args['out_dir'] = sys.argv[4]
    
    main(args)
