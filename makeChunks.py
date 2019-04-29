#!/usr/bin/env python2

# Written by: Andy Yuan
# Email: yuxuan.yuan@outlook.com
# Created date: 01/09/2017
# Last modification date: 03/09/2017

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from math import trunc

def makeChunksByPieces(fasta_file,nreads,ns,path):
    """
    This script is used to split a given fasta file into tiny pieces.
    """
    fasta_file = os.path.abspath(fasta_file)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print '\nSomething is wrong with your output directory! Please check!'
        sys.exit()
    nreads=float(nreads)
    ns=float(ns)
    if nreads%ns == 0:
        step=int(nreads/ns)
    else:
        step=int(trunc(nreads/ns)+1)
    i=1
    j=1

    for seqs in SeqIO.parse(fasta_file, 'fasta'):
        ID=seqs.id
        if ID[0]=="@":
            ID=ID[1:]
        seq=str(seqs.seq)
        start=(i-1)*step+1
        end=i*step
        if end >nreads:
            end = int(nreads)
        fd=open('%s/chunk_%s_%s_%s.fa'%(path,i,start,end),'a')
        if j%step!=0:
            fd.write('>%s\n'%ID)
            fd.write('%s\n'%seq)
        else:
            fd.write('>%s\n'%ID)
            fd.write('%s\n'%seq)
            i+=1
        j+=1

def makeChunksBySize(fasta, size, path):
    """
    This script is used to split a given fasta file into tiny pieces.
    """
    fasta = os.path.abspath(fasta)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print '\nSomething is wrong with your output directory! Please check!'
        sys.exit()
    size = 1000*float(size)
    i=1
    j=0
    for seqs in SeqIO.parse (fasta, 'fasta'):
        ID = seqs.id
        if ID[0]=="@":
            ID=ID[1:]
        seq = str(seqs.seq)
        j+=len(seq)
        if j <= size:
            with open ("%s/chunk_%s.fa"%(path, i), 'a') as fh:
                fh.write('>%s\n'%ID)
                fh.write('%s\n'%seq)
        else:
            j=0
            i+=1
            with open ("%s/chunk_%s.fa"%(path, i), 'a') as fh:
                fh.write('>%s\n'%ID)
                fh.write('%s\n'%seq)

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description="\nSplit assembly into chunks by pieces or size (kbp)")
    parser.add_argument('-v', action='version', version='1.01')
    parser.add_argument('-f', dest='fasta', help='a fasta format file that contains all genome sequences', type = str)
    parser.add_argument('-t', dest='nreads', help='total number of reads in the fasta file', type = int)
    parser.add_argument('-p', dest='npieces', help='number of chunks to split to ', type=int)
    parser.add_argument('-s', dest='size', help='rough size of the chunks you want to split (kbp)', type = float)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.output]:
        if None not in [args.nreads, args.npieces] and None in [args.size]:
            try:
                makeChunksByPieces(args.fasta, args.nreads, args.npieces, args.output)
            except:
                print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
        if None in [args.nreads, args.npieces] and None not in [args.size]:
            try:
                makeChunksBySize(args.fasta, args.size, args.output)
            except:
                print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
    else:
        print
        parser.print_help()
        print
        sys.exit()
