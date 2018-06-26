#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
import itertools

def powerset(l0):
    return itertools.chain.from_iterable(itertools.combinations(l0, r) for r in range(len(l0)+1))

A = [1,2,3,4,5]
B = [3,2,4]
C = [3,4,5]
titles = [partition for partition in powerset(['A', 'B', 'C']) if any(map(None, partition))]
source = [partition for partition in powerset([A, B, C]) if any(map(None, partition))]
dic = {}
for elt in (zip(titles, source)):
    title = "∩".join(str(x) for x in list(elt[0]))
    inter = reduce(lambda x, y: set(x).intersection(y), list(elt[1]))
    dic[title] = inter
    print elt
import pprint 
pprint.pprint(dic)
