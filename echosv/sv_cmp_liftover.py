#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/02/11
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu

Compare two breakpoints of any two SVs, For INS, only compare the inserted breakpoint
'''
import os, timeit, argparse, gzip
from pysam import VariantFile
from echosv.vcf_utils import get_sv_end, get_sv_type, get_sv_len
from echosv.collections import defaultdict 
from intervaltree import IntervalTree
import pandas as pd    

def load_chain(chain_file):
    chain = defaultdict(IntervalTree)
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

def liftover(contig, pos, window_size=500, chain_dict=None, bed_file=None):
    easy_tree = defaultdict(IntervalTree)
    if bed_file is not None:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith("chrom"):
                    continue
                info = line.strip().split()
                if info[3] == "1":
                    easy_tree[info[0]].addi(int(info[1]), int(info[2]) + 1)

    for n in range(window_size+1):
        chain = chain_dict[contig][pos-n]
        # if len(chain) == 1:
        for interval in chain:
            lift_chrom, lift_pos = get_pos(pos-n, interval)
            if bed_file is None or bool(easy_tree[lift_chrom][lift_pos]):
                print(interval)
                return lift_chrom, lift_pos
        # if len(chain) > 1:
        #     raise ValueError(f"Multiple chains found at {contig}:{pos-n}")
        chain = chain_dict[contig][pos+n]
        for interval in chain:
            lift_chrom, lift_pos = get_pos(pos+n, interval)
            if bed_file is None or bool(easy_tree[lift_chrom][lift_pos]):
                print(interval)
                return lift_chrom, lift_pos
        # if len(chain) > 1:
            # raise ValueError(f"Multiple chains found at {contig}:{pos+n}")
    return None, None

def sv_cmp_liftover(refs, fileformat, corr_txt, WINDOW_SIZE=500):
    ref1, sample, platform = os.path.basename(fileformat).split("_")[:3]
    # use the first ref as the baseline and collect ref1-based breakpoints
    sv_pos_dict = {}
    for ref in refs:
        if ref != refs[0]:  # load chain file
            if ref == "grch38" and refs[0] == "chm13":
                chain = load_chain(f"/home/yuz006/disk/t2t_assembly/cancer_public/hg38_to_chm13.chain.gz")
                good_bed = f"/home/yuz006/disk/t2t_assembly/cancer_public/hg38_to_chm13.bed" 
            elif ref == "chm13" and refs[0] == "grch38":
                chain = load_chain(f"/home/yuz006/disk/t2t_assembly/cancer_public/chm13_to_grch38.chain.gz")
                good_bed = f"/home/yuz006/disk/t2t_assembly/cancer_public/chm13_to_grch38.bed"
            elif ref in ["hap1", "hap2"]:
                chain = load_chain(f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_{ref}_{refs[0]}.chain.gz")
                good_bed = f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_{ref}_to_{refs[0]}.bed"
                if not os.path.exists(good_bed):
                    good_bed = f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_{ref}_{refs[0]}.bed"
            elif ref == "dsa" and refs[0] in ["grch38", "chm13"]:
                chain = load_chain(f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_hap1_{refs[0]}.chain.gz")
                good_bed = f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_hap1_{refs[0]}.bed"
                chain_2 = load_chain(f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_hap2_{refs[0]}.chain.gz")
                good_bed_2 = f"/home/yuz006/disk/t2t_assembly/cancer_public/{sample}bl_hap2_{refs[0]}.bed"
        vcf = VariantFile(fileformat.replace(ref1, ref))
        sv_pos_dict[ref] = defaultdict(IntervalTree)
        for record in vcf.fetch():
            chrom1, pos1 = record.chrom, record.pos
            chrom2, pos2 = get_sv_end(record, ins_pseudoPos=False)
            end_pos = pos2
            if ref != refs[0]:
                if ref == "dsa" and chrom1.startswith("h2tg"):
                    chrom1, pos1 = liftover(chrom1, pos1, window_size=50, chain_dict=chain_2, bed_file=good_bed_2)
                    if chrom1 == None or pos1 == None:
                        continue
                    chrom2, pos2 = liftover(chrom2, pos2, window_size=50, chain_dict=chain_2, bed_file=good_bed_2)
                    if chrom2 == None or pos2 == None:
                        continue
                else:
                    chrom1, pos1 = liftover(chrom1, pos1, window_size=50, chain_dict=chain, bed_file=good_bed)
                    if chrom1 == None or pos1 == None:
                        continue
                    chrom2, pos2 = liftover(chrom2, pos2, window_size=50, chain_dict=chain, bed_file=good_bed)
                    if chrom2 == None or pos2 == None:
                        continue
            sv_id = "%".join([record.info["SVTYPE"], chrom1, str(pos1), chrom2, str(pos2), str(record.pos), str(end_pos), str(get_sv_len(record)), record.id])
            sv_pos_dict[ref][chrom1].addi(pos1, pos1+1, (chrom2, pos2, sv_id))
            sv_pos_dict[ref][chrom2].addi(pos2, pos2+1, (chrom1, pos1, sv_id))
        vcf.close()
    matched_sv = set()
    match_items_list = []
    n_matches = 0
    for i in range(len(refs)):
        for ichrom1, svs in sv_pos_dict[refs[i]].items():
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
                for j in range(i+1, len(refs)):
                    ismatch = False
                    jsvs = sv_pos_dict[refs[j]][ichrom1][ipos1-WINDOW_SIZE:ipos1+WINDOW_SIZE]
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
                if len(match_item) - match_item.count(0) == len(refs):
                    n_matches += 1
                match_items_list.append(match_item)
    # output the match items as csv
    print(f"Find {n_matches} pairwise matches for {sample} within {WINDOW_SIZE}bp")
    df = pd.DataFrame(match_items_list, columns=refs)  
    df.to_csv(corr_txt, index=False)       
    return match_items_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--refs", type=str, nargs='+', default=None)
    parser.add_argument("-f", "--fileList", type=str, nargs='+', default=None, required=True)
    parser.add_argument("-o", "--output", type=str, default="./corrMerge.txt")
    parser.add_argument("-w", "--window", type=int, default=500)
    args = parser.parse_args()
    start = timeit.default_timer()

    sv_cmp_liftover(args.refs, args.fileformat, args.output, args.window)
    print("Time used: ", timeit.default_timer()-start)
