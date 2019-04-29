#!/usr/bin/env python

################################################
# Script: fasta2gap_bed.py
# Author: Andy Yuan
# Email:  yuxuan.yuan@outlook.com
# Created date: 05/06/2018
# Last modified date: 05/06/2018
################################################

from Bio import SeqIO
import sys
import os
from argparse import ArgumentParser
  

def fasta2gap_bed(fasta, path, name):
    i=1
    j=1
    with open ('%s/%s.gap.bed'%(path, name), 'w') as gap:
        for seqs in SeqIO.parse(fasta,'fasta'):
            Seq=seqs.seq.upper()
            Len=len(Seq)
            start=0
            seq=Seq
            if Seq[0]!='N':
                while len(seq)>0:
                    gap_s=seq.find('N')
                    start=start+gap_s+1
                    seq=Seq[start-1:]
                    pool = [seq.find('A'),seq.find('T'),seq.find('G'),seq.find('C')]
                    if max(pool)!=-1:
                        nul=[ n for n in pool if n>=0]
                        end= start+ min(nul)-1
                        gap.write('%s\t%s\t%s\tGap_%s\n'%(i, start,end,j))
                        start=end
                        seq=Seq[start:]
                        j+=1    
                        if seq.find('N') == -1:
                            break                     
                    else:
                        end=Len
                        gap.write('%s\t%s\t%s\tGap_%s\n'%(i, start,end,j))
                        j+=1
                        break
            else:
                start=1
                while len(seq)>0:
                    pool = [seq.find('A'),seq.find('T'),seq.find('G'),seq.find('C')]
                    if max(pool)!=-1:
                        nul=[ n for n in pool if n>=0]
                        end=start+min(nul)-1
                        gap.write('%s\t%s\t%s\tGap_%s\n'%(i, start,end,j))
                        j+=1
                        start=end
                        seq=Seq[start:]
                        gap_s=seq.find('N')
                        if gap_s != -1:
                            start=start+gap_s+1
                            seq=Seq[start-1:]
                        else:
                            break
                    else:
                        end=Len
                        gap.write('%s\t%s\t%s\tGap_%s\n'%(i, start,end,j))
                        j+=1
                        break
            i+=1

parser = ArgumentParser(description='covert an input fasta file into a bionano required gap bed file')
parser.add_argument('-v', '--version', action='version', version='1.0')
parser.add_argument('-f', dest='fasta', help="the fasta format file")
parser.add_argument('-p', dest='prefix', help='prefix of outputs. Default: prefix of the input fasta file')
parser.add_argument('-o', dest='output', help='Output directory')
args = parser.parse_args()

if None not in [args.fasta, args.output]:
    args.output = os.path.abspath(args.output)
    if not os.path.isdir(args.output):
        print >> sys.stderr, '\nOops! It seems the path to output directory is not existent. Please check!\n'
        sys.exit(1)
        
    if args.prefix == None:
        args.prefix = args.fasta.rsplit('.',1)[0].split('/')[-1]
    fasta2gap_bed(args.fasta, args.output, args.prefix)
else:
    print 
    parser.print_help()
    print 
    sys.exit(1)    
