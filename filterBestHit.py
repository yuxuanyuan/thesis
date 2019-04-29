#!/usr/bin/env python


# author: Andy Yuan
# email: yuxuan.yuan@outlook.com

import pandas as pd
import sys

def filterbestHits(genelist, blastp, database, n, output):
    n=int(n)
    with open (output, 'w') as fd:
        fd.write('qseqid\tpident\tlength\tevalue\tbitscore\t%s\n'%database)
        hints = pd.read_csv(blastp, sep='\t')
        with open (genelist, 'r') as Genes:
            for i in Genes:
                i = i.strip().split()[0]
                sub = hints[hints['qseqid']==i].reset_index(drop=True)
                if len(sub) == 0:
                    fd.write('%s\tNone\tNone\tNone\tNone\tNone\n' % i)
                else:
                    for j in range(n):
                        try:
                            sortedSub = sub.sort_values(by=['bitscore', 'evalue'],ascending=[False, True]).reset_index(drop=True)
                            info=sortedSub['sseqid'][j]+'_'+sortedSub['subject_description'][j]             
                            fd.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (i, sortedSub['pident'][j], sortedSub['length'][j], sortedSub['evalue'][j], sortedSub['bitscore'][j], info))
                        except KeyError:
                            pass

filterbestHits(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
