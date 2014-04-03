# python script <in_fasta> <out_fasta>

import os
import sys
import pysam

def main(args):
            

    print 'Parsing SAMs...'
    parse_sam(args['sam_in'])
    

def parse_sam(sam_in):
    
    with pysam.Samfile(sam_in, 'r') as in_handle:

        # Get sam iter     
        ref_names = get_ref_name_dict(in_handle)
        
        # Metric counters
        metrics = {}
        metrics['total_count'] = 0
        metrics['unmapped_count'] = 0
        metrics['mapped_count'] = 0
        
        for align in in_handle.fetch():
            
            # Get read size
            read_size = int(str(align.qname).split('size=')[-1])
            metrics['total_count'] += 1 * read_size
            
            # Skip if unmapped
            if align.is_unmapped:
                metrics['unmapped_count'] += 1 * read_size
                continue
            
            # See if mapping passes criteria
            metrics['mapped_count'] += 1 * read_size
            
            ref_id = align.tid
            ref_len = in_handle.lengths[ref_id]
            
            read_name = align.qname
            
            align_len = align.alen
            
            #~ print align.aligned_pairs
            #~ print align.positions
            print read_name, in_handle.getrname(ref_id)
            print align.seq
        
        print metrics['mapped_count']
        print metrics['unmapped_count']
    
    sys.exit()
    
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
    ndn_fasta = os.path.join(out_dir, 'NDN_reads.fasta')
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
    
    # Iterate over J and V sams
    for j_align in j_iter:
        v_align = v_iter.next()
        
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
                        
        # Write the reads N-D-N seq to fasta file
        with open(ndn_fasta, 'a') as out_handle:
            write_fasta(out_handle, read_record.get_header(),
                        read_record.get_insert())
    
    # Print metrics
    for key, value in metrics.iteritems():
        print '{0}: {1}'.format(key, value)
    
    # Close sam iterators
    j_handle.close()
    v_handle.close()
    
    return ref_names, metrics, ndn_fasta

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
    args['sam_in'] = sys.argv[1]
    
    main(args)
