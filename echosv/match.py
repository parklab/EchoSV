#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/01/22
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import gzip, timeit, argparse, os
from pysam import AlignmentFile
from echosv.collect_suppread_lrs import AnnSV as AnnSV_lrs
from echosv.collect_suppread_srs import AnnSV as AnnSV_srs

def match_main(args=None):
    import argparse
    parser = argparse.ArgumentParser(description="Annotate SVs from given dataset.")
    parser.add_argument("--longread", action="store_true", help="Use genotyping for long read data.")
    parser.add_argument("--shortread", action="store_true", help="Use genotyping for short read data.")
    parser.add_argument("-i", "--input", help="Input SV vcf file.", default=None, required=True)
    parser.add_argument("-b", "--bam", nargs="+", help="Bam file.", default=None, required=True)
    parser.add_argument("-o", "--output", help="Output vcf file.")
    parser.add_argument("-q", "--min_mapq", help="Minimum mapping quality, default 1", default=1, type=int)
    parser.add_argument("-d", "--dist_threshold", help="Distance threshold for read alignment, default 500bp", default=500, type=int)
    parser.add_argument("-s", "--ratio_of_seqsize", help="Ratio of sequence size for SVs, default 0.5", default=0.5, type=float)
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode.")
    args = parser.parse_args()
    
    if args.longread:
        anno = AnnSV_lrs(sv_vcf=args.input, bam_list=args.bam, out_file=args.output, min_mapq=args.min_map, dist_threshold=args.dist_threshold, ratio_of_seqsize=args.ratio_of_seqsize, verbose=args.verbose)
    elif args.shortread:
        anno = AnnSV_srs(sv_vcf=args.input, bam_list=args.bam, out_file=args.output, min_mapq=args.min_map, dist_threshold=args.dist_threshold, ratio_of_seqsize=args.ratio_of_seqsize, verbose=args.verbose)
    else:
        raise ValueError("Please specify either --longread or --shortread option.")
    anno.run()
