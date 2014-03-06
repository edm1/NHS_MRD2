#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Will install necessary dependancies, set up folders and download databases
# necessay to use the pipeline.
#
# Should be run using:
# 
# python setup.py ( install | update | updatedb )
#

import os
import sys
import shutil
import stat
import subprocess

def main():
    
    # Check for sudo access
    if not os.getenv("USER") == 'root':
        print ('\nSuperuser access is required to install dependancies. ' +
               'Run again with:')
        print '  sudo python install.py\n'
        sys.exit()
    
    # Write project folders
    folders = ['results']
    exist = []
    for folder in folders:
        if os.path.exists(folder):
            exist.append(folder)
    if len(exist) > 0:
        resp = raw_input('This will overwrite the folders containing ' +
                         'any existing results. Would you like to ' +
                         'continue? (y/n)')
        if resp == 'y':
            for folder in exist:
                shutil.rmtree(folder)
        else:
            print 'Exiting.'
            return 1
    for folder in folders:
        os.mkdir(folder)
        # Change folder permissions to 777
        os.chmod(folder, stat.S_IRWXO)
    
    # Try to install dependancies
    ret = subprocess.call(' '.join([
            'apt-get install',
            'git',
            'python-pip']), shell=True)
            
    if not ret == 0:
        print ('Dependancies were not succesfully installed. They may ' +
               'need to be installed manually.')
        return 1
    
    # Install pysam using pip
    ret = subprocess.call(' '.join([
            'pip install pysam']), shell=True)        
    if not ret == 0:
        print ('pysam was not succesfully installed. Exiting.')
        return 1
    
    # Install bowtie2 and cd-hit-est
    subprocess.call('bash libs/install_programs.sh', shell=True)
    
    
    print 'Done.'

    return 0

if __name__ == '__main__':
	main()

