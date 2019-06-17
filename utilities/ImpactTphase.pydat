#!/usr/bin/python

import os, sys, string, pdb

#
# Note: python's index start from 0
# On Windows, this script has to be run under Cygwin
# The Cygwin can downloaded from: http://cygwin.com/install.html
#---------------------
# (c) Copyright, 2012 by the Regents of the University of California.
#Author: Ji Qiang.
#Description: Energy scan of driven phase (output file engout, col. 5 vs. 1).
#           

#new_row is the line number -1 of that phase variable in the ImpactT.in file.
#new_col is the location of phase variable in that line - 1. 
new_row = 43 #44 row
new_col = 7 #8th column

#the following 3 lines define the range of phase in the scan.
nsteps = 36
inival = 0
delval = 10

val = inival
newValueIndex = 0
while newValueIndex < nsteps:
    test5file = open('ImpactT.in','r')
    lines = test5file.readlines()
    test5file.close()

    # modifiy line 44, number 8
    new_line = string.split(lines[new_row])
    new_line[new_col] = str(val)
    new_line = string.join(new_line) + '\n'

    lines[new_row] = new_line
    test5file2 = open('ImpactT.in', 'w')
    test5file2.writelines(lines)
    test5file2.close()
       
    os.system('./ImpactT > a')
    os.system('tail -1 fort.18 >> tmpeng')
    os.system('echo '+str(val)+' >> tmpphase')
    os.system('tail -1 fort.18')
    newValueIndex = newValueIndex+1
    val = val + delval
    print 'index: ',val

os.system('paste tmpphase tmpeng > engout')
os.system('rm tmpeng tmpphase')
