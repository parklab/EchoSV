#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/01/22
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import gzip, timeit, argparse, os, json
from pysam import AlignmentFile
from echosv.sv_cmp_liftover import sv_cmp_liftover
from echosv.sv_cmp_mbm import SVCombine

def match_main(params=None):
    import argparse
    parser = argparse.ArgumentParser(prog="match", description="Identify concordant and reference-exclusive SVs.")
    parser.add_argument("-i", "--input", help="Input config file, see example as ./test_data/test_colo829_config.json", default="./test_data/test_colo829_config.json")
    parser.add_argument("--multiplat", action="store_true", help="Use multi-platform comparison.")
    parser.add_argument("--merge", action="store_true", help="Merge concordant SVs across references and derive a single VCF.")
    parser.add_argument("--filter", action="store_true", help="Use genotyping-based filter for merging SVs.")
    parser.add_argument("-m", "--min_echo_score", type=float, help="Minimum echo score to consider an SV for matching (default: 0.5).", default=0.5)
    if params is not None:
        if isinstance(params, str):
            params = params.split()
        args = parser.parse_args(params)
    else:
        args = parser.parse_args()
    
    # load config.json file
    with open(args.input, 'r') as f:
        config = json.load(f)
    start_time = timeit.default_timer()
    df = sv_cmp_liftover(config=config)
    print(f"Match SVs through liftover: {timeit.default_timer() - start_time:.2f} seconds")

    start_time = timeit.default_timer()
    SVCombine(config=config, iffilter=args.filter, ifmerge=args.merge, ifmultiplat=args.multiplat, threshold=args.min_echo_score).run()
    print(f"Match SVs through graph-based matching: {timeit.default_timer() - start_time:.2f} seconds")
