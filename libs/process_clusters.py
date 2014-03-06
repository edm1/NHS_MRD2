# pypy script <in_fasta> <out_fasta>

import re
import sys
from libs.bio_file_parsers import fasta_parser

def main(args):
    
    return 0
    
def process_clusters(clus_file, ref_dict, metrics, out_file):
        
        pattern = re.compile('.*_([0-9]+|None),([0-9]+|None):([0-9]+|None),([0-9]+|None);size=([0-9]+)')
        
        with open(out_file, 'w') as out_handle:
            # Write header
            line = ['clus num', 'size', 'proportion', 'v name', 'v del', 'j name', 'j del', 'centroid name', 'centroid seq']
            out_handle.write('\t'.join(line) + '\n')
            
            # Write each cluster
            cluster_num = 1
            for header, seq in fasta_parser(open(clus_file, 'r')):
                # Get J/V and size info from header
                match = pattern.search(header)
                j_num = convert_match(match.group(1))
                j_del = convert_match(match.group(2))
                v_num = convert_match(match.group(3))
                v_del = convert_match(match.group(4))
                size = convert_match(match.group(5))
                # Get V/J name
                j_name = ref_dict['J'].get(j_num, None)
                v_name = ref_dict['V'].get(v_num, None)
                # Convert size to proportion
                proportion = float(size) / metrics['mapped_count']
                # Build line
                line = [cluster_num,
                        size,
                        '{0:.4f}'.format(proportion),
                        v_name,
                        v_del,
                        j_name,
                        j_del,
                        header,
                        seq]
                line = [str(x) for x in line]
                out_handle.write('\t'.join(line) + '\n')
                
                cluster_num += 1
                
        return 0



def convert_match(value):
    """ Tries to convert regex match to an int and returns None if it cant.
    """
    try:
        return int(value)
    except ValueError:
        return None

if __name__ == '__main__':
    
    main()
