#!/usr/bin/env python

import numpy as np
import Bioinfo
import argparse
import gzip

# define command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Plots average qscore for all reads in input FASTQ file")
	parser.add_argument("-f", "--fq_in", help="input fastq file", required=True)
	parser.add_argument("-r", "--read_len", help="read length", required=True, type=int)
	parser.add_argument("-o", "--output", help="output file prefix", required=True)
	return parser.parse_args()

# Assign variables from arguments in the command line
args = get_args()
fq_in = args.fq_in
rlen = args.read_len
out_f = args.output

## Use this cell to populate your numpy array 
#numreads=int(sum(1 for line in gzip.open(fq_in, "r"))/4)
def mean_qscores(fq_in, readlen) -> np.ndarray:
    av_npa = np.zeros((readlen), dtype=int)
    with gzip.open(fq_in,"rt") as fh:
        i, r = 0, 0
        for line in fh:
            i+=1
            if i%4 == 0:
                r+=1
                line=line.rstrip()
                for basepos in range(len(line)):
                    p=int(Bioinfo.convert_phred(line[basepos]))
                    av_npa[basepos]+=p
        av_npa=av_npa/r
    return av_npa

mean=mean_qscores(fq_in, rlen)

f = open(out_f, "a")
f.truncate(0)
f.write("position\tmean\n")
for i in range(len(mean)):
    f.write(str(i+1)+"\t"+str(mean[i]) + "\n")
f.close()