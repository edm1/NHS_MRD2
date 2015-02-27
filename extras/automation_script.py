#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  copy this script to the parent folder (NHS_MRD). Change the the location of
#  the input_files. Locations can be obtained using:
#     ls ~/<location_of_data>/*.fastq.gz
#  and run using:
#     python automation

import os,sys

def main():
    
    input_files = \
'''/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-2-PATIENT-166_S30_L001_R1_001.fastq.gz        
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-2-PATIENT-187_S16_L001_R1_001.fastq.gz        
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-4-PATIENT-166_S31_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-4-PATIENT-187_S17_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-6-PATIENT-166_S32_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-1E-6-PATIENT-187_S18_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-D28-PATIENT-166_S33_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-D28-PATIENT-187_S19_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-DIAGNOSTIC-PATIENT-166_S29_L001_R1_001.fastq.gz
/home/em13383/NHS_MRD/test_data/jack_comp/FR2-MULTIPLEX-DIAGNOSTIC-PATIENT-187_S15_L001_R1_001.fastq.gz
'''
                    
    input_files = input_files.split('\n')
    for in_file in input_files:
        in_reads = in_file.strip()
        out_id = in_reads.split('/')[-1].split('.')[0]
        print out_id
        if os.path.exists('./results/{0}/{0}_bar.pdf'.format(out_id)):
            continue
        arg = ' '.join(['pypy nhs_mrd.py detect',
                        '{0}'.format(in_reads),
                        '{0}'.format(out_id)])
        os.system(arg)
        
    return 0

if __name__ == '__main__':
	main()

