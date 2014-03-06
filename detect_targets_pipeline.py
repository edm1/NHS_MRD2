# python script <in_fasta> <out_fasta>

import os
import sys
from shutil import rmtree
import subprocess

def main(args):
    
    # Set initial variables
    out_dir = args['out_dir']
    
    # Creat out directory
    #~ if not os.path.exists(out_dir):
        #~ os.mkdir(out_dir)
    #~ else:
        #~ resp = raw_input('Out directory already exists. Overwrite? [y/n]')
        #~ if resp == 'y':
            #~ rmtree(out_dir)
            #~ os.mkdir(out_dir)
        #~ else:
            #~ print 'Exiting.'
            #~ return 1
    
    # Derep fastq
    from libs.derep_fastq import derep_fastq
    print 'Dereplicating raw fastq...'
    fastq_file = os.path.join(out_dir, 'derep_reads.fastq')
    #~ derep_fastq(args['in_reads'], fastq_file)
    
    # Run bowtie
    print 'Running bowtie2 alignment...'
    j_sam = os.path.join(out_dir, 'J_align.sam')
    v_sam = os.path.join(out_dir, 'V_align.sam')
    j_bowtie_cmd = 'bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -U {0} -S {1}'
    v_bowtie_cmd = 'bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -U {0} -S {1}'
    j_bowtie_cmd = j_bowtie_cmd.format(fastq_file, j_sam)
    v_bowtie_cmd = v_bowtie_cmd.format(fastq_file, v_sam)
    #~ subprocess.call(j_bowtie_cmd, shell=True)
    #~ subprocess.call(v_bowtie_cmd, shell=True)
    
    # Process sams
    from libs.process_sam import parse_sams
    print 'Processing SAM files...'
    #~ ref_names, metrics, ndn_fasta = parse_sams(j_sam, v_sam, out_dir)
    ndn_fasta = os.path.join(out_dir, 'NDN_reads.fasta') # DEBUG
    
    # Derep fasta
    from libs.derep_fasta import derep_fasta
    print 'Dereplicating N-D-N sequences...'
    ndn_derep = os.path.join(out_dir, 'derep_ndn.fasta')
    derep_fasta(ndn_fasta, ndn_derep)
    
    # Run cd-hit
    print 'Clustering N-D-N sequences...'
    clstr_out = os.path.join(out_dir, 'NDN_clusters.fasta')
    clstr_meta = clstr_out + '.clstr'
    cdhit_cmd = 'cd-hit-est -i {0} -o {1} -c 0.95 -G 1 -d 0 -s 0.9 -r 0 -T 0 -M 2000'
    cdhit_cmd = cdhit_cmd.format(ndn_derep, clstr_out)
    subprocess.call(cdhit_cmd, shell=True)
    
    # Update cd-hit cluster sizes
    from libs.update_cdhit_sizes import update_sizes
    print 'Updating cluster sizes...'
    update_sizes(clstr_out, clstr_meta)
    
    # TODO
    
    return 0


if __name__ == '__main__':
    
    args = {}
    args['in_reads'] = sys.argv[1]
    args['out_dir'] = sys.argv[2]
    
    main(args)
