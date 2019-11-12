#!/usr/bin/python

#######################################################################
# Runchain: a simple, if badly written, code to automate the running
#           of a chain of blockdiag calculations
#######################################################################

import sys
import re
import os
import glob
import copy

#
# Parsing of the input file
#
def rdinp(filename):

    # Set defaults
    dirfile=''
    normcut=1.
    dthresh=1e-6
    refsta=[]
    dmattrans=False
    
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
                dthresh=float((line.split('=')[1]))
                
        elif '$norm_cutoff' in line:
            # Norm cutoff
            if not '=' in line:
                input_error('$norm_cutoff','no argument given')
            else:
                normcut=float((line.split('=')[1]))
        
        elif '$dmat_trans' in line:
            # DFT/MRCI dmat transformation file output
            dmattrans=True
            
        else:
            # Unknown keyword
            print('\n','Error parsing line: '+line,'\n')
            sys.exit()
            
    return dirfile,normcut,dthresh,refsta,dmattrans

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
# Determination of the quantum chemistry code used
#
def getqctype(directory):

    qctype=None

    # DFT/MRCI?
    for f in ('auxbasis','out3','control'):
        check=glob.glob(directory+'/'+f)
        if len(check)!=0:
            qctype='dftmrci'

    # Columbus MRCI?
    for f in ('ciudgin.drt1','ciudgsm.drt1.sp'):
        check=glob.glob(directory+'/'+f)
        if len(check)!=0:
            qctype='colmrci'

    # Exit if we could not determine the calculation type
    if qctype==None:
        print('\n Error: the quantum chemistry code used could not be determined\n')
        sys.exit()
    
    return qctype
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
# Determine the names of the disp. determinant files
#
def getdispdets():

    dispdets=glob.glob('disp/*.bin')
    dispdets.sort(key=natural_keys)

    return dispdets

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

#
# Determine the ref. states from the blockdiag output at
# the previous geometry
#
def getrefstates(lastlbl):

    statelist=[]
    
    filename=lastlbl+'.log'
    
    with open(filename,"r") as logfile:
        log=logfile.readlines()

    for line in log:
        if 'Selected state' in line:
            state=line.split()[-1:]
            statelist.append(int(state[0]))

    return statelist

#
# Read the disp. state energies
#
def rddispen(qctype,directory):
    
    ener=[]

    # DFT/MRCI
    if qctype=='dftmrci':
        filename=directory+'/out3'
        if not os.path.exists(filename):
            print('\nCould not find the DFT/MRCI output file'+filename+'\n')
            sys.exit()
        with open(filename,"r") as outfile:
            lines=outfile.readlines()
        for line in lines:
            if 'DFTCI  ' in line:
                string=line.split()
                ener.append(string[string.index('DFTCI')+1])

    # Columbus MRCI
    elif qctype=='colmrci':
        # Check to see if ciudgsm.drt1.sp is present in
        # directory/. or directory/LISTINGS
        found=False
        filename=directory+'/ciudgsm.drt1.sp'
        if os.path.exists(filename):
            found=True
            with open(filename,"r") as outfile:
                lines=outfile.readlines()
        else:
            filename=directory+'/LISTINGS/ciudgsm.drt1.sp'
            if os.path.exists(filename):
                found=True
                with open(filename,"r") as outfile:
                    lines=outfile.readlines()

        # Exit if ciudgsm.drt1.sp could not be found
        if not found:
            print('\nCoule not find the MRCI output file ciudgsm.drt1.sp\n')
            sys.exit()

        # Parse ciudgsm.drt1.sp
        for line in lines:
            if 'eci       =' in line:
                string=line.split()
                ener.append(string[string.index('=')+1])

    return ener
    
#
# Write a blockdiag input file
#
def wrbdinp(filename,i,dispdets,refcurr,lastlbl,dthresh,
            normcut,dispen,dmattrans):

    # Open the blockdiag input file
    f=open(filename,"w+")
    
    # MO files
    f.write('$mos_ref=ref/mos.dat \n')
    f.write('$mos_disp=disp/mos.dat \n')
    
    # Ref. determinants
    f.write('\n$dets_ref')
    for k in range(len(refcurr)):
        f.write('\n ref/det.1.'+str(refcurr[k])+'.bin '+str(k+1))
    f.write('\n$end\n')

    # Disp. determinants
    f.write('\n$dets_disp')
    for k in range(len(dispdets)):
        f.write('\n '+dispdets[k]+' '+str(k+1))
    f.write('\n$end\n')

    # Ref. ADT matrix file
    if i!=1:
        f.write('\n$ref_trans='+lastlbl+'.log\n')

    # Disp. state energies
    f.write('\n$energies')
    for k in range(len(dispen)):
        f.write('\n'+str(dispen[k])+' '+str(k+1))
    f.write('\n$end\n')
    
    # Hadamard screening threshold
    f.write('\n$dthresh='+str(dthresh)+'\n')

    # Norm cutoff
    f.write('\n$norm_cutoff='+str(normcut)+'\n')

    # DFT/MRCI dmat transformation file output
    if dmattrans==True:
        f.write('\n$dmat_trans\n')
    
    # Close the blockdiag input file
    f.close()

#
# Check whether or not a blockdiag calculation ran to completion
#
def blockdiag_status(filename):

    completed=False

    with open(filename,"r") as outfile:
        lines=outfile.readlines()

    for line in lines:
        if 'CPU Time:' in line:
            completed=True
        
    return completed

#
# Main routine
#
# Read the name of the input file
if len(sys.argv) < 2:
    print('Error: no input file given')
    sys.exit()
infile=str(sys.argv[1])

# Parse the input file
dirfile,normcut,dthresh,refsta,dmattrans=rdinp(infile)

# Parse the directory file
dirlist=rddirfile(dirfile)

# Try and determine the quantum chemistry code used
qctype=getqctype(dirlist[0])

# Make sure that the determinant directories exists and are
# empty
if os.path.isdir('ref'):
    os.system('rm -r ref')
os.system('mkdir ref')

if os.path.isdir('disp'):
    os.system('rm -r disp')
os.system('mkdir disp')

# Initialise lastlbl
lastlbl=''

# Perform the chain of blockdiag calculations
for i in range(1,len(dirlist)):
    
    # Ouput our progress
    lbl=lastdir(dirlist[i])+'_diab'
    print('\n'+lbl)
    
    # Extract any tar.gz determinant files
    extract_dets(dirlist[i-1],'ref')
    extract_dets(dirlist[i],'disp')
    
    # Copy over the MO files
    os.system('cp '+dirlist[i-1]+'/mos.dat ref/')
    os.system('cp '+dirlist[i]+'/mos.dat disp/')
    
    # Determine the ref. states to include
    if i==1:
        refcurr=refsta
    else:
        refcurr=getrefstates(lastlbl)

    # Get the list of all disp. determinant files
    dispdets=getdispdets()

    # Get the disp. geometry adiabatic energies
    dispen=rddispen(qctype,dirlist[i])
    
    # Write the blockdiag input file
    bdinpfile=lbl+'.inp'
    wrbdinp(bdinpfile,i,dispdets,refcurr,lastlbl,dthresh,
            normcut,dispen,dmattrans)

    # Run the blockdiag calculation
    inputfile=lbl+'.inp'
    os.system('blockdiag.x '+inputfile)

    # Exit here if the blockdiag calculation failed
    completed=blockdiag_status(lbl+'.log')
    if not completed:
        print('\n Blockdiag failed...\n')
        sys.exit()
    
    # Update lastlbl
    lastlbl=copy.deepcopy(lbl)
    
