#!/usr/bin/env python

to_get = {} # key: subject, value: list of queries

import sys

ins = sys.argv[2:]
goa = sys.argv[1]

for l in ins:
    seen = set()
    for line in open(l):
        ll = line.split()
        query, subject1 = ll[0], ll[1]
        subject=subject1.split('|')[1]
        if query in seen: continue
        seen.add(query)
        if subject in to_get:
            to_get[subject].append(query)
        else:
            to_get[subject] = [query]
go_dict = {}

with open (goa, 'r') as GOA:
    for line in GOA:
        if line.startswith('!'): continue
        ll = line.split('\t')
        gene = ll[1]
        go = ll[4]
        if gene in to_get:
            if gene in go_dict:
                go_dict[gene].append(go)
            else:
                go_dict[gene] = [go]

for g in sorted(go_dict):
    queries = to_get[g]
    for q in queries:
        print q + '\t' + '\t'.join(go_dict[g])
