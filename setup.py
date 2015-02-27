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
    
    # Extract argument
    try:
        possible_arg = ['install', 'update', 'updatedb']
        task = sys.argv[1].lower()
        assert task in possible_arg
    except:
        print '\nArgument is invalid. Command must be one of the following:'
        print '  sudo python setup.py install'
        print '  python setup.py update'
        print '  python setup.py updatedb\n'
        sys.exit()
    
    # Ig Database URL
    db_url = ('ftp://ftp.ncbi.nih.gov/blast/executables/'
              'igblast/release/database/human*')
    
    #
    # Install
    #
    
    if task == 'install':
        # Check for sudo access
        if not os.getenv("USER") == 'root':
            print ('\nSuperuser access is required to install dependancies. ' +
                   'Run again with:')
            print '  sudo python setup.py install\n'
            sys.exit()
        
        # Write project folders
        folders = ['database', 'tmp', 'results', 'comparisons']
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
            
        # Download databases
        ret = subprocess.call('wget -P ./database {0}'.format(db_url),
                              shell=True)
        if not ret == 0:
            print 'Failed to download Ig databases'
            return 1
        
        # Try to install dependancies
        ret = subprocess.call(' '.join([
                'apt-get install',
                'git',
                'ncbi-blast+',
                'pypy',
                'python-pip']), shell=True)
                
        if not ret == 0:
            print ('Dependancies were not succesfully installed. They may ' +
                   'need to be installed manually.')
            return 1
        
        # Install docopt using pip
        ret = subprocess.call(' '.join([
                'pip install docopt==0.6.1']), shell=True)        
        
        if not ret == 0:
            print ('Docopt was not succesfully installed. Exiting.')
            return 1
        
        print 'Done.'
    
    elif task == 'update':

        # Try to git pull
        ret = subprocess.call('git pull', shell=True)
        if not ret == 0:
            ret2 = subprocess.call('git fetch --all', shell=True)
            ret3 = subprocess.call('git reset --hard origin/master',shell=True)
            if not ret2 == 0 and ret3 == 0:
                print 'Script was unable to update pipeline.'
                return 1
        print '\nPipeline updated.'
        
    elif task == 'updatedb':
        ret = subprocess.call('wget -P ./database {0}'.format(db_url),
                              shell=True)
        if not ret == 0:
            print 'Failed to update Ig databases'
            return 1
        print 'Done'    

    return 0

if __name__ == '__main__':
	main()

