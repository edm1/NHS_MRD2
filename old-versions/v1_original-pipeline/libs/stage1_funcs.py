#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from shutil import rmtree
import os
import time
import subprocess

def prepare_folders(args):
    """ Will check if a folder with the same ID name already exists and offer
        to overwrite it. Then will make some new project folders.
    """
    
    # Check if project folder already exists
    if os.path.exists('./tmp/{0}/'.format(args['<id>'])):
        resp = raw_input('This project folder already exists, would you ' + 
                         'like to overwrite it? (y/n): ')
        if resp == 'y':
            try:
                rmtree('./tmp/{0}/'.format(args['<id>']))
                rmtree('./results/{0}/'.format(args['<id>']))
            except OSError:
                pass
        else:
            print 'Run again with different --id parameter'
            sys.exit()
                
    # Make new project folders
    os.mkdir('./tmp/{0}/'.format(args['<id>']))
    os.mkdir('./results/{0}/'.format(args['<id>']))

def log_parameters(log_file, arg_dict):
    """Appends all the run parameters to the log file.
    """
    with open(log_file, 'a') as out_handle:
        out_handle.write('\n--Run parameters--\n')
        for option in arg_dict:
            out_handle.write('{0}: {1}\n'.format(option, arg_dict[option]))

def print_timestamp(message):
    """Prints a timestamp and message for user.
    """
    print "{0} {1}".format(time.strftime('%H:%M:%S',time.localtime()),
                           message)

def run_D_blast(in_file, top_seg):
    """Function to run BLAST. Required to multiprocess.
    """
    subprocess.call(' '.join(['blastn',
                              '-db ./database/human_gl_D',
                              '-query {0}'.format(in_file),
                              '-out {0}_D.blast'.format(in_file),
                              '-evalue 100000',
                              '-dust no',
                              '-word_size 5',
                              '-penalty -4',
                              '-outfmt 5',
                              '-max_target_seqs {0}'.format(top_seg)
                              ]), shell=True)

def run_J_blast(in_file, top_seg):
    """Function to run BLAST. Required to multiprocess.
    """
    subprocess.call(' '.join(['blastn',
                              '-db ./database/human_gl_J',
                              '-query {0}'.format(in_file),
                              '-out {0}_J.blast'.format(in_file),
                              '-evalue 1000',
                              '-dust no',
                              '-word_size 7',
                              '-penalty -3',
                              '-outfmt 5',
                              '-max_target_seqs {0}'.format(top_seg)
                              ]), shell=True)

def run_V_blast(in_file, top_seg):
    """Function to run BLAST. Required to multiprocess.
    """
    subprocess.call(' '.join(['blastn',
                              '-db ./database/human_gl_V',
                              '-query {0}'.format(in_file),
                              '-out {0}_V.blast'.format(in_file),
                              '-evalue 20',
                              '-dust no',
                              '-word_size 9',
                              '-penalty -2',
                              '-outfmt 5',
                              '-max_target_seqs {0}'.format(top_seg)
                              ]), shell=True)
                              
def first_entry_name(top_list):
    """Given a list of the top hits it will return the name of the first hit.
    """
    if len(top_list) > 0:
        return top_list[0].hit
    else:
        return None

def detailed_string(entry, hit_num):
    """Produces VDJ string for the detailed output.
    """
    align_len = entry.q_end - entry.q_start + 1
    return ' '.join(['\t{0}:'.format(hit_num),
                      '(-{0}){1}(-{2})'.format(entry.del5, entry.hit,
                                               entry.del3),
                      'e={0:.2e}'.format(entry.e),
                      'qalign={0}..{1}'.format(entry.q_start, entry.q_end),
                      'len={0}'.format(align_len),
                      'mismatch={0}'.format(align_len - entry.identity)
                      ]) + '\n'

def reverse_complement_DNA(seq, reverse=True):
    """ Will return the (reverse) complement of a DNA seq
    """
    
    map_d = {'A': 'T',
             'T': 'A',
             'C': 'G',
             'G': 'C',
             'N': 'N'}
    
    # Get complimentary bases
    comp_seq = [ map_d[nuc] for nuc in seq ]
    
    # Reverse if True
    if reverse:
        comp_seq.reverse()
    
    return ''.join(comp_seq)
    
    
