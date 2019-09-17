#!/usr/bin/python

import sys
import re
import os
import glob

#
# Parsing of the input file
#
def rdinp(filename):

    # Set defaults
    dirfile=''
    normcut=1.
    dthresh=1e-6
    refsta=[]

    # Read in the input file
    with open(filename, 'r') as infile:
        inp=infile.readlines()

    # Remove comment lines marked with '#' or '!'
    inp = [re.split('#|!', ln)[0] for ln in inp]

    # Remove blank lines
    filtered_inp = filter(blankline, inp)

    # Parse the input
    for line in filtered_inp:

        if '$dirfile' in line:
            # Directory file
            if not '=' in line:
                input_error('$dirfile','no argument given')
            else:
                dirfile=(line.split('=')[1])
                dirfile=clean_keyword(dirfile)
                
        elif '$refstates' in line:
            # Reference geometry states
            if not '=' in line:
                input_error('$refstates','no argument given')
            else:
                k=line.index('=')+1
                refsta=(line[k:].split(','))
                refsta=list(map(int, refsta))
                
        elif '$dthresh' in line:
            # Hadamard screening threshold
            if not '=' in line:
                input_error('$dthresh','no argument given')
            else:
                dthesh=float((line.split('=')[1]))

        elif '$norm_cutoff' in line:
            # Norm cutoff
            if not '=' in line:
                input_error('$norm_cutoff','no argument given')
            else:
                normcut=float((line.split('=')[1]))

        else:
            # Unknown keyword
            print('\n','Error parsing line: '+line,'\n')
            sys.exit()
            
    return dirfile,normcut,dthresh,refsta

#
# blankline
#
def blankline(line):
    
    blank=['','\n']

    if (line in blank):
        return False
    else:
        return True

#
# Input error handling
#
def input_error(keyword,reason):

    print('\n','Error parsing keyword: '+keyword)
    print('\n','Reason: '+reason,'\n')
    sys.exit()

#
# Clean keywords up
#
def clean_keyword(keyword):

    # Remove any trailing newline characters
    clean_keyword=keyword.rstrip('\n')

    # Remove any blank spaces at the start
    clean_keyword=clean_keyword.strip()

    return clean_keyword
    
#
# Parsing of the directory file
#
def rddirfile(dirfile):

    # Initialisation
    dirlist=[]
    
    # Make sure that the directory file exists
    try:
        f=open(dirfile)
        f.close()
    except:
        print('File does not exist: '+dirfile)

    # Read in the directories file
    with open(dirfile, 'r') as infile:
        inp=infile.readlines()
        
    # Remove comment lines marked with '#' or '!'
    inp = [re.split('#|!', ln)[0] for ln in inp]

    # Clean-up the directory names
    for i in range(len(inp)):
        inp[i]=clean_keyword(inp[i])
        
    # Remove blank lines
    filtered_inp = filter(blankline,inp)
    
    # Parse the directory names
    for line in filtered_inp:
        dirlist.append(line)

    return dirlist

#
# Extract determinant files
#
def extract_dets(detdir,outdir):
        
    # Extract the determinant files in the directory
    # detdir to the directory outdir
    targz=glob.glob(detdir+'/*.bin.tar.gz')
    for i in range(len(targz)):
        os.system('tar -zxvf '+targz[i]+' -C'+outdir)

#
# Last directory in a path
#
def lastdir(path):

    filtered=filter(blankline, path.split('/'))
    lbl=str(filtered[-1:])[2:-2]    

    return lbl

#
# Write a blockdiag input file
#
def wrbdinp(filename):

    print('FINISH WRITING WRBDINP!')
    sys.exit()
    
#
# Main routine
#
# Read the name of the input file
if len(sys.argv) < 2:
    print('Error: no input file given')
    sys.exit()
infile=str(sys.argv[1])

# Parse the input file
dirfile,normcut,dthresh,refsta=rdinp(infile)

# Parse the directory file
dirlist=rddirfile(dirfile)

# Set the reference geometry directory
refdir=dirlist[0]

# Make sure that the directory outdir exists and is
# empty
if os.path.isdir('ref'):
    os.system('rm -r ref')
os.system('mkdir ref')

# Extract and tar.gz determinant files in the reference directory
extract_dets(refdir,'ref')

# Perform the chain of blockdiag calculations
for i in range(1,len(dirlist)):

    # Current directory
    currdir=dirlist[i]

    # Ouput our progress
    lbl=lastdir(currdir)
    print('\n'+lbl)
    
    # Extract any tar.gz determinant files
    extract_dets(currdir,'./')

    # Write the blockdiag input file
    bdinpfile=lbl+'_diab.inp'
    wrbdinp(bdinpfile)
