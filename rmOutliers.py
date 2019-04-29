#!/usr/bin/env python2

# Written by: Andy Yuan	
# Email: yuxuan.yuan@outlook.com


from __future__ import print_function, division
import numpy as np
import matplotlib.pylab as plt
from PyAstronomy import pyasl
from math import log, exp
import pandas as pd
import argparse
import os
import sys


def rmOutliers(myfile, chrName, method, pvalue, miniSize, path):
    
    pvalue=float(pvalue)
    miniSize = int(miniSize)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print ('\nSomething is wrong with your output directory! Please check!\n')
        sys.exit(1)

    fd = pd.read_csv(myfile, sep = '\t', names=['chr', 'start', 'end', 'length'])
    
    if chrName!='all':
        sub=fd[fd['chr']==chrName].reset_index(drop=True)
    else:
        sub=fd	 
    if len(sub)<1:
        print('\nOpps! It seems that there is no SV in the input chromosome or the input chrName does not match with the chr in the file!\n')
        sys.exit(1)
	
    x = np.array([log(float(x)) for x in sub['length']])
    
    Q3=np.percentile(x, 75)
    Q1=np.percentile(x, 25)
    upBound = 2.5*Q3 -1.5*Q1
    lowBound = 2.5*Q1 -1.5*Q3
    sum=0
    for j in range(len(sub)):
        if x[j] <lowBound:
            sum+=1
        if x [j] > upBound:
            sum+=1
    #print (sum)
    
    #plt.boxplot(x)
    # Apply the generalized ESD
    if sum < 2:
        r = pyasl.generalizedESD(x, 2, pvalue, fullOutput=True)
    else:
        r = pyasl.generalizedESD(x, sum, pvalue, fullOutput=True)
    print("The checked chromosome is: ", chrName)
    print("Number of outliers found by Tukey's method: ", sum)
    print("Number of outliers found by generalized ESD: ", r[0])
    # print("Indices of outliers: ", r[1])
    # Plot the "data"
    # plt.plot(x, 'b.')
    
    # and mark the outliers.
    # for i in range(r[0]):
    #    plt.plot(r[1][i], x[r[1][i]], 'rp')
    # plt.show()
    if r[0]>0:
        smallest = exp(x[r[1][-1]])
    
        print("The samllest outlier is ", smallest)
    
        with open ('%s/%s.%s.rmOutliers.txt' % (path, chrName, method),'w') as fh:
            for j in range(len(sub)):
                if smallest > lowBound:
                    if sub['length'][j]<smallest and sub['length'][j] >= miniSize:
                        fh.write('%s\t%s\t%s\t%s\n' % (sub['chr'][j], sub['start'][j], sub['end'][j], sub['length'][j]))
                else:
                    if sub['length'][j] >= miniSize:
                        fh.write('%s\t%s\t%s\t%s\n' % (sub['chr'][j], sub['start'][j], sub['end'][j], sub['length'][j]))
    else:
        with open ('%s/%s.%s.rmOutliers.txt' % (path, chrName, method),'w') as fh:
            for j in range(len(sub)):
                if sub['length'][j]>=miniSize:
                    fh.write('%s\t%s\t%s\t%s\n' % (sub['chr'][j], sub['start'][j], sub['end'][j], sub['length'][j]))            


if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='remove outliers by using generalized ESD method.')
    parser.add_argument('-v', '--version',action='version', version='1.0')
    parser.add_argument('-i', dest='info', help='SV information file contains chr, start, end and length. seprated by tab', type = str)
    parser.add_argument('-c', dest='chrName', help='name of the targeted chromosome. Can be all or particular chromosome', type = str)
    parser.add_argument('-m', dest='method', help='method used to call SVs. Such as Pindel', type= str)
    parser.add_argument('-p', dest='pvalue', help='pvalue used to calculate outliers. Default: 0.05', default=0.05, type=float)
    parser.add_argument('-s', dest='size', help='SV cutoff size (bp). Default: 1.', default=1, type=int)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.info, args.chrName, args.method, args.output]:
        run = rmOutliers(args.info, args.chrName, args.method, args.pvalue, args.size, args.output)
    else:
        print
        parser.print_help()
        print
        sys.exit()
