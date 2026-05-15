#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2023/11/09
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import os, re, glob, gzip, argparse
from collections import defaultdict 
from intervaltree import IntervalTree, Interval
from pysam import VariantFile, tabix_index
from echosv.vcf_utils import _ispass, _checkVCF, _rowcount, get_sv_end, get_sv_len, get_sv_type, _getVAF
from echosv.bed_utils import loadBed
from time import perf_counter   # runtime
import pandas as pd
import numpy as np

    
class Merge:
    def __init__(self,
                 compFileList=None, 
                 outFile=None,
                 atol=500, 
                 sizetol=0.5,
                 includebed=None, 
                 checksvtype=False,
                 passonly=True,
                 verbose=False,
                 dowrite=False,
                 checkdup=False,
                 ins_pseudoPos=True,
                 ) -> None:
        if compFileList is None:
            raise ValueError("No vcf files input")
        self.comFileList = []
        if outFile is None:
            self.outFile = outFile
        else:
            self.outFile = outFile.split(".")[0] + ".txt"
        self.atol = atol
        self.sizetol = sizetol
        self.includebed = includebed
        self.checksvtype = checksvtype
        self.passonly = passonly
        self.verbose = verbose
        self.ins_pseudoPos = ins_pseudoPos
        if checkdup:
            for file in compFileList:
                cmpFile, n_cmp = _checkVCF(file, passonly=self.passonly)
                if not os.path.exists(cmpFile + ".tbi"):
                    tabix_index(cmpFile, preset="vcf", force=True)
                self.comFileList.append(cmpFile)
        else:
            self.comFileList = compFileList
        # vars
        self.sv_tree = {}
        self.cluster = defaultdict(list)
        self.dowrite = dowrite

    @staticmethod
    def _generate_id(record):
        chrom1, pos1 = record.chrom, record.pos
        chrom2, pos2 = get_sv_end(record, ins_pseudoPos=True)
        sv_id = record.id if record.id is not None else "."
        vaf, alt_supp = _getVAF(record, 0)
        sv_id = "%".join([chrom1, str(pos1), get_sv_type(record), str(pos2), str(get_sv_len(record)), sv_id, str(vaf), str(alt_supp)])
        return sv_id
    
    @staticmethod
    def _parse_sv_id(sv_id: str):
        """Parse the sv_id format produced by generate_index(record)."""
        parts = sv_id.split("%")
        # svtype, chrom1, pos1, chrom2, pos2, svlen, record_id, vaf
        if len(parts) < 8:
            raise ValueError(f"Unexpected sv_id format: {sv_id}")
        chrom1, pos1, svtype, pos2, svlen, rec_id, vaf, alt_supp = parts[:8]
        return chrom1, int(pos1), svtype, int(pos2), int(svlen), rec_id, float(vaf), int(alt_supp)

    def run(self):
        # load in regions
        if self.includebed is not None and os.path.isfile(self.includebed):
            self.include_tree = loadBed(self.includebed)

        for file in self.comFileList:
            if not os.path.isfile(file):
                raise FileNotFoundError(f"vcf file {file} not found or invalid")
            vcf = VariantFile(file)
            if len(file.split("/")) < 2:
                sample = file
            elif file.split("/")[-1].startswith("SMHT"):
                sample = file
            else:
                sample = file.split("/")[-2] + "-" + os.path.basename(file).split("_")[2] 
            self.sv_tree[sample] = defaultdict(IntervalTree)
            for record in vcf.fetch():
                if self.passonly and not _ispass(record):
                    continue
                chrom1, pos1 = record.chrom, record.pos
                chrom2, pos2 = get_sv_end(record, ins_pseudoPos=self.ins_pseudoPos)
                sv_id = self._generate_id(record)
                self.sv_tree[sample][chrom1].addi(pos1, pos1+1, (chrom2, pos2, sv_id))
                self.sv_tree[sample][chrom2].addi(pos2, pos2+1, (chrom1, pos1, sv_id))
            vcf.close()

        samples = list(self.sv_tree.keys())
        used_sv = {sample: set() for sample in samples}
        matched_sv = {sample: set() for sample in samples}
        n_matches = 0
        for i, file in enumerate(self.comFileList):
            vcf = VariantFile(file)
            sample = samples[i]
            for sv in vcf.fetch():
                if self.passonly and not _ispass(sv):
                    continue
                chrom2, pos2 = get_sv_end(sv, ins_pseudoPos=self.ins_pseudoPos)
                sv_id = self._generate_id(sv)
                if sv_id in used_sv[sample]:    # check if the sv is matched
                    continue
                support_list = [sample]
                backup_list = [sv_id]
                used_sv[sample].add(sv_id)
                for j in range(i, len(samples)):
                    jsvs = self.sv_tree[samples[j]][sv.chrom][sv.pos-self.atol:sv.pos+self.atol]
                    for jsv in jsvs:
                        jpos1 = jsv.begin
                        jchrom2, jpos2, jsv_id = jsv.data
                        if jsv_id in used_sv[samples[j]]:
                            continue
                        is_match = self.check_match((chrom2, pos2, sv_id.split("%")[2], int(sv_id.split("%")[4])), jsv)
                        if is_match:
                            support_list.append(samples[j])
                            backup_list.append(jsv_id)
                            used_sv[samples[j]].add(jsv_id)
                            matched_sv[samples[j]].add(jsv_id)
                            matched_sv[sample].add(sv_id)
                self.cluster[sv.chrom].append([sv.pos, sv, support_list, backup_list])
                n_matches += 1
        print(f"Total {n_matches} matched SVs found among {len(samples)} samples.")
        if self.dowrite:    # for pairwise comparison only
            n_cmp = 0
            n_fp = 0
            vcf_cmp = VariantFile(self.comFileList[1])
            tp_outh = VariantFile(self.comFileList[1].replace(".vcf", "_FP.vcf"), "w", header=vcf_cmp.header)
            for record in vcf_cmp.fetch():
                if self.passonly and not _ispass(record):
                    continue
                n_cmp += 1
                sv_id = self._generate_id(record)
                if sv_id not in matched_sv[samples[1]]:
                    tp_outh.write(record)
                    n_fp += 1
            tp_outh.close()
            vcf_cmp.close()

            n_base = 0
            n_fn = 0
            vcf_base = VariantFile(self.comFileList[0])
            fn_outh = VariantFile(self.comFileList[1].replace(".vcf", "_FN.vcf"), "w", header=vcf_base.header)
            for record in vcf_base.fetch():
                if self.passonly and not _ispass(record):
                    continue
                n_base += 1
                sv_id = self._generate_id(record)
                if sv_id not in matched_sv[samples[0]]:
                    fn_outh.write(record)
                    n_fn += 1
            fn_outh.close()
            vcf_base.close()
            n_cmp_tp = n_cmp - n_fp
            n_base_tp = n_base - n_fn
            print(f"{self.comFileList[1]} has {n_cmp_tp} TPs from {n_base} base SVs and {n_cmp} cmp SVs, precision {n_cmp_tp/n_cmp} recall {n_base_tp/n_base} f1 {2*(n_cmp_tp/n_cmp)*(n_base_tp/n_base)/((n_cmp_tp/n_cmp)+(n_base_tp/n_base))}")

    def check_match(self, i_info, j_interval):  # return True/False, distance
        chrom_i, pos_i, svtype_i, svlen_i = i_info
        chrom_j, pos_j, id_j = j_interval.data
        svlen_j = int(id_j.split("%")[4])
        svtype_j = id_j.split("%")[2]
        if chrom_i != chrom_j:
            return False
        if self.checksvtype and svtype_i != svtype_j:
            return False
        else:
            if {svtype_i, svtype_j} == {"INS", "DEL"}:
                return False
        if abs(pos_i - pos_j) < self.atol and self.size_cmp(svlen_i, svlen_j):
            return True
        return False

    def size_cmp(self, svlen_i, svlen_j):
        if max(svlen_i, svlen_j) == 0:
            return True
        if min(svlen_i, svlen_j)/max(svlen_i, svlen_j) > self.sizetol:
            return True
        return False
    
def quick_cmp(vcf1, vcf2, withid=False):
    sv_list = []
    for vcf in [vcf1, vcf2]:
        sv_set = set()
        for entry in VariantFile(vcf):
            chrom2, pos2 = get_sv_end(entry, ins_pseudoPos=False)
            if withid:
                sv_set.add((entry.chrom, entry.pos, chrom2, pos2, entry.id, entry.info["SVTYPE"], entry.alts[0]))
            else:
                sv_set.add((entry.chrom, entry.pos, chrom2, pos2, entry.info["SVTYPE"], entry.alts[0]))
        sv_list.append(sv_set)
    print("sv1 - sv2: ", len(sv_list[0] - sv_list[1]))
    a = list(sv_list[0] - sv_list[1])
    print("sv1 - sv2", sorted(a, key=lambda x: (x[0], x[1])))

def merge(out_vcf, cmp_object, ifnew=False, command_line=None):
    base_vcf = VariantFile(cmp_object.comFileList[0])

    if ifnew:
        f = VariantFile(out_vcf.replace(".txt", ".vcf.gz"), "w", header=None)
        f.header.add_line('##fileformat=VCFv4.2')
        f.header.add_line('##reference={}'.format(os.path.basename(cmp_object.comFileList[0]).split("_")[0]))
        # write merging params
        f.header.add_line(f'##Merging tolerance distance={cmp_object.atol}, Size similarity ratio>{cmp_object.sizetol}')
        f.header.add_line('##Command line={}'.format(command_line))
        # write contigs
        for contig in base_vcf.header.contigs:
            length = base_vcf.header.contigs[contig].length
            f.header.add_line('##contig=<ID={},length={}>'.format(contig, length))
        # write sv types
        f.header.add_line('##ALT=<ID=DEL,Description="Deletion">')
        f.header.add_line('##ALT=<ID=DUP,Description="Duplication">')
        f.header.add_line('##ALT=<ID=INV,Description="Inversion">')
        f.header.add_line('##ALT=<ID=INS,Description="Insertion">')
        f.header.add_line('##ALT=<ID=BND,Description="Translocation">')
        # write filter
        f.header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
        # write format
        f.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        # write info
        f.header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
        f.header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
        f.header.add_line('##INFO=<ID=CHROM2,Number=1,Type=String,Description="End chromosome of the variant described in this record">')
        f.header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')
        f.header.add_line('##INFO=<ID=AF,Number=1,Type=Float,Description="Variant Frequency">')
        f.header.add_line('##INFO=<ID=Support,Number=.,Type=String,Description="Caller and Technology supporting the structural variation">')
        f.header.add_line('##INFO=<ID=Backup,Number=.,Type=String,Description="Backup information when merging SVs, chrom-pos-svtype-stop-alt">')
        f.header.add_sample('SAMPLE')
        # f.header.add_line('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE')
        for chrom in base_vcf.header.contigs:
            sort_svs = sorted(cmp_object.cluster[chrom], key=lambda x: x[0])
            for sv in sort_svs:
                pos, sv_record, support_list, backup_list = sv
                chrom2, pos2 = get_sv_end(sv_record, ins_pseudoPos=False)
                svtype = get_sv_type(sv_record)
                svlen = get_sv_len(sv_record)
                af, _ = _getVAF(sv_record, 0)
                new_record = f.new_record()
                new_record.chrom = sv_record.chrom
                new_record.pos = sv_record.pos
                new_record.id = sv_record.id 
                new_record.ref = sv_record.ref if sv_record.ref is not None else "N"
                new_record.alts = sv_record.alts if sv_record.alts is not None else ("<{}>".format(svtype),)
                new_record.filter.add("PASS")
                new_record.info["SVTYPE"] = svtype
                new_record.info["CHROM2"] = chrom2
                new_record.stop = pos2
                new_record.info["SVLEN"] = svlen
                new_record.info["AF"] = af
                new_record.info["Support"] = ",".join(support_list)
                new_record.info["Backup"] = ",".join(backup_list)
                f.write(new_record)
        f.close()
        base_vcf.close()
    else:
        out_vcf = VariantFile(out_vcf.replace(".txt", ".vcf.gz"), "w", header=base_vcf.header)
        
        for chrom in base_vcf.header.contigs:
            sort_svs = sorted(cmp_object.cluster[chrom], key=lambda x: x[0])
            for sv in sort_svs:
                pos, sv_record, support_list, backup_list = sv
                out_vcf.write(sv_record)
        base_vcf.close()
        out_vcf.close()
    tabix_index(out_vcf.replace(".txt", ".vcf.gz"), preset="vcf", force=True)


def extract_truthset(raw_merge_vcf,
                     output_truthset,
                     n_supp=4,
                     n_platform=2,
                     removeGap=True,
                     gaps_bed=None,
                     command_line=None):
    vcf = VariantFile(raw_merge_vcf)
    out_vcf = VariantFile(output_truthset, 'w', header=vcf.header)
    out_vcf.header.add_line('##Command line={}'.format(command_line))
    n_total = 0
    n_extract = 0

    ref = os.path.basename(raw_merge_vcf).split("_")[0]
    print(f"Extracting truth set from {raw_merge_vcf} with at least {n_supp} supporting callers and {n_platform} platforms, removeGap={removeGap}")
    if removeGap:
        sample = os.path.basename(raw_merge_vcf).split("_")[1]
        if ref == "grch38" or "hg38":
            gaps_bed = "../../beds/hg38_no_alt_Ns.bed"
            gap_tree = loadBed(gaps_bed)
        elif os.path.exists(f"../../beds/{sample}bl_{ref}_Ns.bed"):
            gaps_bed = f"../../beds/{sample}bl_{ref}_Ns.bed"
            gap_tree = loadBed(gaps_bed)
        elif gaps_bed is None:
            print("Warning: --gaps-bed not provided; skipping gap-proximity filter.")
            removeGap = False
        elif not os.path.isfile(gaps_bed):
            print(f"Warning: gaps BED file not found at {gaps_bed}; skipping gap-proximity filter.")
            removeGap = False
        else:
            gap_tree = loadBed(gaps_bed)
    gap_dict = 1000
    # select SV with supporting >= n_supp and n_platform
    for record in vcf:
        n_total += 1
        # check if around reference gaps
        if removeGap:
            if ref != "chm13":
                overlap_gap = gap_tree[record.chrom][record.pos - gap_dict: record.pos + gap_dict]
                if len(overlap_gap):
                    if len(record.info["Support"]) >= n_supp:
                        # print(overlap_gap, record.chrom, record.pos, record.info["Support"])
                        pass
                    continue
                chrom2, pos2 = get_sv_end(record, ins_pseudoPos=False)
                overlap_gap = gap_tree[chrom2][pos2 - gap_dict: pos2 + gap_dict]
                if len(overlap_gap):
                    if len(record.info["Support"]) >= n_supp:
                        # print(overlap_gap, chrom2, pos2, record.info["Support"])
                        pass
                    continue        
        supp_platforms = set()
        for supp in record.info["Support"]:
            supp_platforms.add(supp.split("-")[1])
        pos1_list = []
        vaf_list = []
        svtype_list = []
        svlen_list = []
        for sv in record.info["Backup"]:
            chrom1, pos1, svtype, pos2, svlen, rec_id, vaf, alt_supp = sv.split("%")
            pos1_list.append(int(pos1))
            vaf_list.append(float(vaf))
            svtype_list.append(svtype)
            svlen = int(svlen)
            if svlen >= 50:
                svlen_list.append(svlen)
        if len(set(record.info["Support"])) >= n_supp and len(set(supp_platforms)) >= n_platform:
            record.info["AF"] = np.nanmedian(vaf_list)
            if len(set(svtype_list)) > 1:
                if "INS" in svtype_list and "DUP" in svtype_list:
                    record.info["SVTYPE"] = "DUP"
                elif "BND" in svtype_list:
                    svtype_list = [x for x in svtype_list if x != "BND"]
                    # find the most common type
                    most_svtype = max(set(svtype_list), key=svtype_list.count)
                    record.info["SVTYPE"] = most_svtype
                else:
                    most_svtype = max(set(svtype_list), key=svtype_list.count)
                    record.info["SVTYPE"] = most_svtype
            if record.info["SVTYPE"] != "BND":
                record.info["SVLEN"] = int(np.nanmedian(svlen_list))
                if record.info["SVTYPE"] != "INS":
                    record.stop = record.pos + record.info["SVLEN"] 
            else:
                record.info["SVLEN"] = 0
            out_vcf.write(record)
            n_extract += 1
    vcf.close()
    out_vcf.close()
    # tabix_index(output_truthset, preset="vcf", force=True)
    print(f"Extracted {n_extract} out of {n_total} records with at least {n_supp} supporting callers", raw_merge_vcf)


def merge_main(params=None):
    parser = argparse.ArgumentParser(prog="merge", description="Compare/merge SV call sets from a same reference genome.")
    parser.add_argument("-i", "--files", type=str, nargs='+', required=True)
    parser.add_argument("-o", "--output", type=str, default="./sv_cmp_output.txt")
    parser.add_argument("-a", "--atol", type=int, default=500)
    parser.add_argument("-s", "--sizetol", type=float, default=0.5)
    parser.add_argument("-b", "--includebed", type=str, default=None)
    parser.add_argument("-c", "--checksvtype", action="store_true", default=False)
    parser.add_argument("-d", "--checkdup", action="store_true", default=False)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-w", "--dowrite", action="store_true", default=False)
    # for merging
    parser.add_argument("--merge", action="store_true", default=False)  # for merging multiple VCFs
    parser.add_argument("--new", action="store_true", default=False, help="If given, create variant object from scratch") 
    # for extracting high-confidence SVs
    parser.add_argument("--extract", action="store_true", default=False, help="If given, extract high-confidence SVs with at least 4 supporting callers and 2 platforms")
    parser.add_argument("--gapbed", type=str, default=None, help="Path to a BED file of reference gap (N) regions for gap-proximity filtering (used with --extract)")
    if params is not None:
        if isinstance(params, str):
            params = params.split()
        args = parser.parse_args(params)
    else:
        args = parser.parse_args()

    if args.extract:
        if args.files is None or len(args.files) != 1:
            raise ValueError("Exactly one merged VCF file is required for extracting truth set.")
        raw_merge_vcf = args.files[0]
        command_line = "echosv merge -i {} -o {} --extract".format(raw_merge_vcf, args.output)
        extract_truthset(raw_merge_vcf, args.output, n_supp=4, n_platform=2, removeGap=True, gaps_bed=args.gapbed, command_line=command_line)
        exit(0)

    if args.files is None or len(args.files) < 2:
        raise ValueError("At least two VCF files are required for comparison.")

    start_time = perf_counter()
    model = Merge(compFileList=args.files,  
            outFile=args.output,
            atol=args.atol,
            sizetol=args.sizetol,
            includebed=args.includebed,
            checksvtype=args.checksvtype,
            verbose=args.verbose,
            dowrite=args.dowrite,
            checkdup=args.checkdup,
            )
    model.run()
    print("Comparison finished in {:.2f} seconds.".format(perf_counter() - start_time))
    if args.merge:
        print("Merging vcfs into {}".format(args.output))
        command_line = "echosv merge " + " ".join(["--files"] + args.files) + "--output {} --atol {} --sizetol {} --merge --new".format(args.output, args.atol, args.sizetol)
        merge(args.output, model, args.new, command_line=command_line)
        
