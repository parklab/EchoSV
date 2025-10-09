#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2023/11/09
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import gzip
from collections import defaultdict 
from intervaltree import IntervalTree

def loadBed(bedFile, tree=None, max_size=None, goodonly=False):   # dict of IntervalTrees, each interval: chrom/start/end
    if tree is None:
        tree = defaultdict(IntervalTree)
        
    if bedFile.endswith('.gz'):
        f = gzip.open(bedFile, 'rt')
    else:
        f = open(bedFile, 'r')
    if "centromere" in bedFile or "telomere" in bedFile:
        tolerance = 500
    else:
        tolerance = 0
    for line in f:
        if line.startswith('#') or line.startswith("chrom"):
            continue
        info = line.strip().split()
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        if (goodonly and info[3] != "1"):
            continue
        # tree[chrom].addi(start, end + 1)
        if max_size is None or (max_size is not None and end - start <= max_size):
            tree[chrom].addi(max(0, start - tolerance), end + tolerance + 1)
    f.close()
    return tree # keys: chrom, values: [(start1, en1), (start2, end2), ...]