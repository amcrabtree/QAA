#!/usr/bin/env python

#################### MODULES ######################
import sys
import re 

#################### ARGPARSE ######################
infile=sys.argv[1]

#################### MAIN ######################
## check if current read is mapped
fh=open(infile, "r")
m_reads=0
u_reads=0
for line in fh:
    if line.startswith("@") == False:
        # store bitwise flag (2nd position) in SAM header
        flag=int(line.strip("\n").split("\t")[1])
        if((flag & 256) != 256): # bit 256 is not true, therefore NO 2ndary alignment
            if ((flag & 4) != 4): # bit 4 is false, therefore read mapped
                    m_reads+=1
            else: # read unmapped
                    u_reads+=1
     
print("mapped reads: ", m_reads)
print("umapped reads: ", u_reads)