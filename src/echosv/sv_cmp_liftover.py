#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/02/11
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu

Compare two breakpoints of any two SVs, For INS, only compare the inserted breakpoint
'''
import os, timeit, argparse, gzip, glob
from pysam import VariantFile
from echosv.vcf_utils import get_sv_end, get_sv_type, get_sv_len
from echosv.bed_utils import loadBed
from collections import defaultdict 
from intervaltree import IntervalTree
import pandas as pd    

def load_chain(chain_file, chain=None):
    if chain is None:
        chain = defaultdict(IntervalTree)
    # search for chain files that match the pattern as chain_file
    files = glob.glob(chain_file)
    if len(files) == 0:
        raise FileNotFoundError(f"No chain files found for pattern: {chain_file}")
    for chain_file in files:
        with gzip.open(chain_file, 'rt') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                contig, strand, start, end, ref, ref_start, ref_end, size = line.strip().split("\t")
                chain[contig].addi(int(start), int(end)+1, (int(strand), ref, int(ref_start), int(ref_end), int(size)))
    return chain

def get_pos(pos, interval):
    offset = pos - interval.begin
    if interval.data[0] == 1:   # reverse
        return interval.data[1], interval.data[3] - offset
    else:
        return interval.data[1], interval.data[2] + offset

def liftover(contig, pos, window_size=500, chain_dict=None, bed_tree=None):
    for n in range(window_size+1):
        chain = chain_dict[contig][pos-n]
        # if len(chain) == 1:
        for interval in chain:
            lift_chrom, lift_pos = get_pos(pos-n, interval)
            if bed_tree is None or bool(bed_tree[lift_chrom][lift_pos]):
                return lift_chrom, lift_pos
        # if len(chain) > 1:
        #     raise ValueError(f"Multiple chains found at {contig}:{pos-n}")
        chain = chain_dict[contig][pos+n]
        for interval in chain:
            lift_chrom, lift_pos = get_pos(pos+n, interval)
            if bed_tree is None or bool(bed_tree[lift_chrom][lift_pos]):
                return lift_chrom, lift_pos
        # if len(chain) > 1:
            # raise ValueError(f"Multiple chains found at {contig}:{pos+n}")
    return None, None

def sv_cmp_liftover(config, WINDOW_SIZE=500):
    # use the first ref1 as the baseline and collect ref1-based breakpoints
    sv_pos_dict = {}
    ref1 = config["refs"]["1"]
    for ref_index, ref in config["refs"].items():
        if ref != ref1:
            chain = load_chain(config["chains"][f"{ref_index}_to_1"])
            bedfiles = glob.glob(config["chains"][f"{ref_index}_to_1"].replace(".chain.gz", ".bed"))
            goodTree = []
            for bedfile in bedfiles:
                goodTree.append(loadBed(bedfile, goodonly=True))
        vcf = VariantFile(config["vcfs"][ref_index])
        sv_pos_dict[ref] = defaultdict(IntervalTree)
        for record in vcf.fetch():
            chrom1, pos1 = record.chrom, record.pos
            chrom2, pos2 = get_sv_end(record, ins_pseudoPos=False)
            end_pos = pos2
            if ref != ref1:
                if "h2tg" in chrom1:
                    chrom1, pos1 = liftover(chrom1, pos1, window_size=50, chain_dict=chain, bed_tree=goodTree[1])
                else:
                    chrom1, pos1 = liftover(chrom1, pos1, window_size=50, chain_dict=chain, bed_tree=goodTree[0])
                if chrom1 == None or pos1 == None:
                    continue
                if "h2tg" in chrom2:
                    chrom2, pos2 = liftover(chrom2, pos2, window_size=50, chain_dict=chain, bed_tree=goodTree[1])
                else:
                    chrom2, pos2 = liftover(chrom2, pos2, window_size=50, chain_dict=chain, bed_tree=goodTree[0])
                if chrom2 == None or pos2 == None:
                    continue
            sv_id = "%".join([record.info["SVTYPE"], chrom1, str(pos1), chrom2, str(pos2), str(record.pos), str(end_pos), str(get_sv_len(record)), record.id])
            sv_pos_dict[ref][chrom1].addi(pos1, pos1+1, (chrom2, pos2, sv_id))
            sv_pos_dict[ref][chrom2].addi(pos2, pos2+1, (chrom1, pos1, sv_id))
        vcf.close()
    matched_sv = set()
    match_items_list = []
    n_matches = 0
    for i in range(len(config["refs"])):
        for ichrom1, svs in sv_pos_dict[config["refs"][str(i+1)]].items():
            for sv in svs:  # each sv is an interval
                ipos1 = sv.begin
                ichrom2, ipos2, isv_id = sv.data
                if isv_id in matched_sv:
                    continue
                # check if the sv is matched
                match_item = []
                for index in range(0, i):
                    match_item.append(0)
                match_item.append(isv_id)
                matched_sv.add(isv_id)
                for j in range(i+1, len(config["refs"])):
                    ismatch = False
                    jsvs = sv_pos_dict[config["refs"][str(j+1)]][ichrom1][ipos1-WINDOW_SIZE:ipos1+WINDOW_SIZE]
                    if len(jsvs):
                        for jsv in jsvs:
                            jpos1 = jsv.begin
                            jchrom2, jpos2, jsv_id = jsv.data
                            if jsv_id in matched_sv:
                                continue
                            if abs(jpos1 - ipos1) < WINDOW_SIZE and jchrom2 == ichrom2 and abs(jpos2 - ipos2) < WINDOW_SIZE:
                                match_item.append(jsv_id)
                                matched_sv.add(jsv_id)
                                ismatch = True
                                break
                    if not ismatch:
                        match_item.append(0)
                    else:
                        n_matches += 1

                if len(match_item) - match_item.count(0) == len(config["refs"]):
                    n_matches += 1
                match_items_list.append(match_item)
    # output the match items as csv
    # print(f"Find {n_matches} pairwise matches within {WINDOW_SIZE}bp")
    df = pd.DataFrame(match_items_list, columns=config["refs"].values())  
    df.to_csv(config["output"].replace(".txt", "_liftover.txt"), index=False)       
    return match_items_list
