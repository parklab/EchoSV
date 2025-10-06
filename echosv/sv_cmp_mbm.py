#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2024/08/04
@Author: Yu Wei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

from pysam import VariantFile, VariantHeader, tabix_index
import argparse, os, timeit
from scipy.sparse import lil_matrix, csr_matrix
import pandas as pd
import numpy as np
from echosv.vcf_utils import get_sv_len, get_sv_end

def get_suppRead(vcf_file, sampleidx=0, remove0suppread=False):
    '''
    Get the supporting reads for each SV in the vcf file
    '''
    svs_dict = {}   # id: [supp_reads]
    suppread_dict = {}  # read: [sv_id]
    vcf = VariantFile(vcf_file)
    sampleName = list(vcf.header.samples)[sampleidx]
    for record in vcf:
        chrom2, pos2 = get_sv_end(record)
        sv_id = str(record.pos) + "%" + str(pos2) + "%" + str(get_sv_len(record)) + "%" + record.id
        supp_reads = record.samples[sampleName]['RNAMES']
        if not isinstance(supp_reads, tuple):
            if supp_reads == "." or supp_reads == None:   # no supporting reads
                supp_reads = []
            else:
                supp_reads = [supp_reads]
        if remove0suppread and len(supp_reads) == 0:
            continue
        svs_dict[sv_id] = supp_reads
        for read in supp_reads:
            if read not in suppread_dict:
                suppread_dict[read] = []
            suppread_dict[read].append(sv_id)
    return svs_dict, suppread_dict

def get_match_score(suppRead1, suppRead2, threshold=0.5):   # Jaccard similarity
    supp_set1 = set(suppRead1)
    supp_set2 = set(suppRead2)
    if len(supp_set1) == 0 or len(supp_set2) == 0:
        return 0
    score = len(supp_set1 & supp_set2) / len(supp_set1 | supp_set2)
    return score if score >= threshold else 0

def get_max_match_dense(match_matrix):  # ndarray input
    row_ind, col_ind = linear_sum_assignment(match_matrix, maximize=True)
    matching = list(zip(row_ind, col_ind))
    total_cost = match_matrix[row_ind, col_ind].sum()
    # only keep the valid matches
    matching = [(row, col) for row, col in matching if match_matrix[row, col] > 0]
    return matching, total_cost

class SVCombine:
    '''
    create the sv match table for the input vcf files, and output the table to a csv file.
    rows: SVs from the base vcf file
    columns: SVs from the comparison vcf files
    '''
    def __init__(self, corr_txt, 
                 fileformat, 
                 refs, 
                 output,
                 verbose=False, 
                 threshold=0.3, 
                 ifhighconfi=True):
        self.corr_txt = corr_txt
        self.fileformat=fileformat
        self.reflist = refs
        self.fileList = []
        self.out_file = output
        self.verbose = verbose
        self.threshold = threshold
        self.ifhighconfi = ifhighconfi
        self.svs_items = {ref: set() for ref in self.reflist}

    def _loadMatches(self):
        pre_matches = set()
        pre_svs = {ref: set() for ref in self.reflist}
        corr_matrix = pd.read_csv(self.corr_txt, sep=",", header=0)
        for index, row in corr_matrix.iterrows():
            match_ids = []
            for ref in self.reflist:
                if row[ref] != "0":
                    # pos1%pos2%svlen%sv_id
                    sv_id = row[ref].split("%")[-4] + "%" + row[ref].split("%")[-3] + "%" + row[ref].split("%")[-2] + "%" + row[ref].split("%")[-1]
                    if sv_id in pre_svs[ref]:
                        raise ValueError("Find duplicate svs in the given corr matrix, please check!", sv_id, ref)    
                    pre_svs[ref].add(sv_id)
                    match_ids.append(sv_id)
            if len(match_ids) > 1:  
                for i in range(len(match_ids)):
                    for j in range(i+1, len(match_ids)):
                        pre_matches.add(match_ids[i] + "-" + match_ids[j])
        return pre_matches, pre_svs

    def run(self):
        # collect all the svs
        if len(self.fileformat) == 1:
            self.fileformat = self.fileformat[0]
            refo, sample = os.path.basename(self.fileformat).split("_")[:2]
            for ref in self.reflist:
                self.fileList.append(self.fileformat.replace(refo, ref))
                vcf = VariantFile(self.fileformat.replace(refo, ref))
                for record in vcf:
                    _, pos2 = get_sv_end(record)
                    self.svs_items[ref].add(str(record.pos) + "%" + str(pos2) + "%" + str(get_sv_len(record)) + "%" + record.id)
                vcf.close()
        else:
            for file_i in range(len(self.fileformat)):
                self.fileList.append(self.fileformat[file_i])
                vcf = VariantFile(self.fileformat[file_i])
                for record in vcf:
                    _, pos2 = get_sv_end(record)
                    self.svs_items[self.reflist[file_i]].add(str(record.pos) + "%" + str(pos2) + "%" + str(get_sv_len(record)) + "%" + record.id)
                vcf.close()
        if not self.corr_txt is None:
            pre_matches, pre_svs = self._loadMatches()
        # create a table with 12 columns 
        result_df = pd.DataFrame(columns=self.reflist)
        n_matches = 0
        n_novel_match_easy = 0
        n_novel_match_difficult = 0
        # base is the first vcf file
        for base_index in range(0, len(self.reflist)):
            base_ref = self.reflist[base_index]
            base_sv_dict, base_suppRead_dict = get_suppRead(self.fileList[base_index])
            base_svs = list(base_sv_dict.keys())
            base_table = []
            for i, sv in enumerate(base_sv_dict):
                base_table.append([sv] + [0]*(len(self.reflist)-base_index-1))
            for cmp_index in range(base_index+1, len(self.reflist)):
                cmp_ref = self.reflist[cmp_index]
                cmp_sv_dict, cmp_suppRead_dict = get_suppRead(self.fileList[cmp_index])
                cmp_svs = list(cmp_sv_dict.keys())
                match_matrix = self.get_matrix(base_sv_dict, cmp_sv_dict, base_suppRead_dict, cmp_suppRead_dict)
                if self.ifhighconfi:
                    # ont similarity matrix
                    sv_dict1, sv_suppRead_dict1 = get_suppRead(self.fileList[base_index], sampleidx=1)
                    sv_dict2, sv_suppRead_dict2 = get_suppRead(self.fileList[cmp_index], sampleidx=1)
                    match_matrix = self.get_matrix(sv_dict1, sv_dict2, sv_suppRead_dict1, sv_suppRead_dict2, match_matrix)
                    # ill similarity matrix
                    # if os.path.getsize(self.fileList[base_index].replace("lr", "sr")) and os.path.getsize(self.fileList[cmp_index].replace("lr", "sr")):
                    sv_dict1, sv_suppRead_dict1 = get_suppRead(self.fileList[base_index].replace("lr", "sr"), sampleidx=0)
                    sv_dict2, sv_suppRead_dict2 = get_suppRead(self.fileList[cmp_index].replace("lr", "sr"), sampleidx=0)
                    match_matrix = self.get_matrix(sv_dict1, sv_dict2, sv_suppRead_dict1, sv_suppRead_dict2, match_matrix)
                if self.corr_txt is not None:   # add pre_matches info
                    # remove the previous matches from the lil matrix
                    base_to_del = set()    # index: count
                    cmp_to_del = set()
                    for i in range(len(base_svs)):
                        for j in range(len(cmp_svs)):
                            if base_svs[i] + "-" + cmp_svs[j] in pre_matches and match_matrix[i, j] > 0:
                                base_table[i][cmp_index-base_index] = cmp_svs[j]
                                n_matches += 1
                                base_to_del.add(i)
                                cmp_to_del.add(j)
                    print("Remove previous matches from the match matrix", base_ref, cmp_ref, len(base_to_del), len(cmp_to_del))
                    # print([base_svs[i] for i in base_to_del], [cmp_svs[i] for i in cmp_to_del])
                    rows_to_keep = np.delete(np.arange(len(base_svs)), list(base_to_del))
                    cols_to_keep = np.delete(np.arange(len(cmp_svs)), list(cmp_to_del))
                    match_matrix = match_matrix[rows_to_keep, :][:, cols_to_keep]
                    base_new_order = [base_svs[i] for i in range(len(base_svs)) if i not in base_to_del]
                    cmp_new_order = [cmp_svs[i] for i in range(len(cmp_svs)) if i not in cmp_to_del]

                matches = self.get_matches(match_matrix)
                if self.verbose:
                    print(f"Matched {len(matches)} out of {len(base_sv_dict)} base svs and {len(cmp_sv_dict)} cmp svs from {base_ref} and {cmp_ref}")
                for base_sv_index, cmp_sv_index in matches:
                    n_matches += 1
                    if self.corr_txt is not None:
                        base_sv = base_new_order[base_sv_index]
                        cmp_sv = cmp_new_order[cmp_sv_index]
                        if base_sv not in pre_svs[base_ref] or cmp_sv not in pre_svs[cmp_ref]:
                            # print("Find novel match for SV in difficult region", base_ref, base_sv, base_sv not in pre_svs[base_ref], cmp_ref, cmp_sv, cmp_sv not in pre_svs[cmp_ref], "score", match_matrix[base_sv_index, cmp_sv_index], len(base_sv_dict[base_sv]), len(cmp_sv_dict[cmp_sv]))
                            n_novel_match_difficult += 1
                        else:
                            # print("Find novel match for SV in easy region", base_ref, base_sv, cmp_ref, cmp_sv, "score", match_matrix[base_sv_index, cmp_sv_index], len(base_sv_dict[base_sv]), len(cmp_sv_dict[cmp_sv]))
                            n_novel_match_easy += 1
                        base_sv_index = base_svs.index(base_sv)
                        cmp_sv_index = cmp_svs.index(cmp_sv)
                    base_table[base_sv_index][cmp_index-base_index] = cmp_svs[cmp_sv_index]
            # remove the existing rows
            existing_svs = result_df[base_ref].drop_duplicates().tolist()
            # print(f"Existing SVs: {len(existing_svs)}, {existing_svs}")
            base_table = [row for row in base_table if row[0] not in existing_svs]
            # convert the base_table to a dataframe
            new_df = pd.DataFrame(base_table, columns=self.reflist[base_index:])
            # fill the previous columns with NA
            for i in range(base_index):
                new_df[self.reflist[i]] = 0
            new_df = new_df[self.reflist]
            result_df = pd.concat([result_df, new_df], axis=0, ignore_index=True)
        if self.verbose:
            print(f"Find {n_matches} matches, {n_novel_match_easy} novel matches in easy region and {n_novel_match_difficult} novel matches in difficult region.")
        result_df.to_csv(self.out_file, sep=",", index=False)
        if self.reflist == ["hap1", "hap2"]:    # dsa merging
            # write into a new file
            if self.ifhighconfi:
                self.output_sv(result_df)
            else:
                self.output_sv_quick(result_df)
        
    
    def get_matrix(self, svDict1, svDict2, supp_read_dict1, supp_read_dict2, match_matrix=None):
        len1 = len(svDict1)
        len2 = len(svDict2)
        supp_read_set = set(supp_read_dict1.keys()) & set(supp_read_dict2.keys())
        start = timeit.default_timer()
        # print("Creating the match matrix...", len1, len2, len(supp_read_set))
        if match_matrix is None:
            match_matrix = lil_matrix((len1, len2), dtype=float)
        cal_state_matrix = lil_matrix((len1, len2), dtype=int)   # 0: not calculated, 1: calculated
        sv1_ids = {}
        index = 0
        for id in svDict1.keys():
            sv1_ids[id] = index
            index += 1
        sv2_ids = {}
        index = 0
        for id in svDict2.keys():
            sv2_ids[id] = index
            index += 1  
        for i in supp_read_set:
            for sv1_id in supp_read_dict1[i]:
                for sv2_id in supp_read_dict2[i]:
                    sv1_ind = sv1_ids[sv1_id]
                    sv2_ind = sv2_ids[sv2_id]
                    if cal_state_matrix[sv1_ind, sv2_ind] == 0:
                        score = get_match_score(svDict1[sv1_id], svDict2[sv2_id], threshold=0)
                        cal_state_matrix[sv1_ind, sv2_ind] = 1
                        if score > match_matrix[sv1_ind, sv2_ind]:
                            match_matrix[sv1_ind, sv2_ind] = score
        return match_matrix

    def get_matches(self, match_matrix):
        start = timeit.default_timer()
        # set values of match_matrix (lil_matrix) lower than the threshold to 0
        match_matrix[match_matrix < self.threshold] = 0
        matches, max_score = get_max_match_sparse(match_matrix.tocsr())
        # print(max_score, "Time taken to get the max match: ", timeit.default_timer() - start)
        return matches
    
    def output_sv(self, result_df):
        # out_vcf = self.out_file.replace(".txt", ".vcf.gz")
        out_vcf = self.fileList[0].replace("hap1", "dsa")
        chosen_svs = set()
        for mode in ["lr", "sr"]:
            n_hap_spec = 0
            vcf1 = VariantFile(self.fileList[0].replace("lr", mode))
            # if "MatchScore" not in vcf1.header.info:
            #     vcf1.header.info.add("MatchScore", number=1, type="Float", description="Match score between two SVs")
            if "MatchSV" not in vcf1.header.info:
                vcf1.header.info.add("MatchSV", number=1, type="String", description="Matched SV id")
            n_out = 0
            vcf2 = VariantFile(self.fileList[1].replace("lr", mode))
            header = VariantHeader()
            for record in vcf1.header.records:
                if record.key != "contig":
                    header.add_record(record)
            for contig_name, contig_info in vcf1.header.contigs.items():
                header.contigs.add(contig_name, length=contig_info.length)
            for contig_name, contig_info in vcf2.header.contigs.items():
                header.contigs.add(contig_name, length=contig_info.length)
            if mode == "lr":
                for i in range(4):
                    header.add_sample(f"SAMPLE{i}") 
            else:   # sr
                for i in range(2):
                    header.add_sample(f"SAMPLE{i}")
            outVCF = VariantFile(out_vcf.replace("lr", mode), "w", header=header)
            for svItem in vcf1:
                chrom2, pos2 = get_sv_end(svItem)
                sv_id = str(svItem.pos) + "%" + str(pos2) + "%" + str(get_sv_len(svItem)) + "%" + svItem.id
                # check if the sv is in the result_df["hap1"]
                # print(sv_id, result_df[result_df["hap1"] == sv_id]["hap2"].values[0])
                match_sv = result_df[result_df["hap1"] == sv_id]["hap2"].values[0]
                if match_sv != 0:
                    svItem.info["MatchSV"] = match_sv
                    outVCF.write(svItem)
                    chosen_svs.add(sv_id)
                    n_out += 1
                elif mode == "lr" and self.check_normal_supp(svItem):
                    outVCF.write(svItem)
                    chosen_svs.add(sv_id)
                    n_out += 1
                    n_hap_spec += 1
                elif mode == "sr" and sv_id in chosen_svs:
                    outVCF.write(svItem)
                    n_out += 1
                    n_hap_spec += 1
                else:
                    # print("potential FP", sv_id, match_sv) 
                    pass
            for svItem in vcf2.fetch():
                chrom2, pos2 = get_sv_end(svItem)
                sv_id = str(svItem.pos) + "%" + str(pos2) + "%" + str(get_sv_len(svItem)) + "%" + svItem.id
                match_sv = result_df[result_df["hap2"] == sv_id]["hap1"].values[0]
                if match_sv == 0:
                    if mode == "lr" and not self.check_normal_supp(svItem):
                        continue
                    if mode == "sr" and sv_id not in chosen_svs:
                        continue
                    new_record = header.new_record(contig=svItem.contig, start=svItem.pos, stop=svItem.stop, id=svItem.id, alleles=svItem.alleles, qual=svItem.qual, filter=svItem.filter)
                    for key in svItem.info.keys():
                        new_record.info[key] = svItem.info[key]
                    for index, key in enumerate(svItem.samples.keys()):
                        new_record.samples[f"SAMPLE{index}"]["AF"] = svItem.samples[key]["AF"]
                        if svItem.samples[key]["AF"] == 0:
                            new_record.samples[f"SAMPLE{index}"]["RNAMES"] = ""
                            if "RHAP" in svItem.samples[key]:
                                new_record.samples[f"SAMPLE{index}"]["RHAP"] = ""
                        else: 
                            new_record.samples[f"SAMPLE{index}"]["RNAMES"] = ",".join(svItem.samples[key]["RNAMES"])
                            if "RHAP" in svItem.samples[key]:
                                new_record.samples[f"SAMPLE{index}"]["RHAP"] = svItem.samples[key]["RHAP"]
                    new_record.pos = svItem.pos
                    new_record.stop = svItem.stop
                    new_record.info["SVTYPE"] = svItem.info["SVTYPE"]
                    if not (new_record.chrom == svItem.chrom and new_record.pos == svItem.pos and new_record.stop == svItem.stop):
                        raise ValueError("The new record is not the same as the original record", str(svItem), str(new_record))
                    outVCF.write(new_record)
                    n_hap_spec += 1
                    if mode == "lr":
                        chosen_svs.add(sv_id)
                    n_out += 1
            outVCF.close()
            vcf1.close()
            vcf2.close()
            print("Output {} from {} SVs to {}, filtered haplotype-specific {}".format(n_out, len(result_df), out_vcf, n_hap_spec))
            tabix_index(out_vcf.replace("lr", mode), force=True, preset="vcf")
    
    def output_sv_quick(self, result_df):
        out_vcf = self.fileList[0].replace("hap1", "dsa")
        vcf1 = VariantFile(self.fileList[0])
        # if "MatchScore" not in vcf1.header.info:
        #     vcf1.header.info.add("MatchScore", number=1, type="Float", description="Match score between two SVs")
        if "MatchSV" not in vcf1.header.info:
            vcf1.header.info.add("MatchSV", number=1, type="String", description="Matched SV id")
        n_out = 0
        vcf2 = VariantFile(self.fileList[1])
        header = VariantHeader()
        for record in vcf1.header.records:
            if record.key != "contig":
                header.add_record(record)
        for contig_name, contig_info in vcf1.header.contigs.items():
            header.contigs.add(contig_name, length=contig_info.length)
        for contig_name, contig_info in vcf2.header.contigs.items():
            header.contigs.add(contig_name, length=contig_info.length)
        for i in range(2):
            header.add_sample(f"SAMPLE{i}") 
        outVCF = VariantFile(out_vcf, "w", header=header)
        for svItem in vcf1:
            chrom2, pos2 = get_sv_end(svItem)
            sv_id = str(svItem.pos) + "%" + str(pos2) + "%" + str(get_sv_len(svItem)) + "%" + svItem.id
            # check if the sv is in the result_df["hap1"]
            match_sv = result_df[result_df["hap1"] == sv_id]["hap2"].values[0]
            if match_sv != 0:
                svItem.info["MatchSV"] = match_sv
                outVCF.write(svItem)
                n_out += 1
            elif self.check_normal_supp(svItem):
                outVCF.write(svItem)
                n_out += 1
        for svItem in vcf2.fetch():
            chrom2, pos2 = get_sv_end(svItem)
            sv_id = str(svItem.pos) + "%" + str(pos2) + "%" + str(get_sv_len(svItem)) + "%" + svItem.id
            match_sv = result_df[result_df["hap2"] == sv_id]["hap1"].values[0]
            if match_sv == 0:
                if not self.check_normal_supp(svItem):
                    continue
                new_record = header.new_record(contig=svItem.contig, start=svItem.pos, stop=svItem.stop, id=svItem.id, alleles=svItem.alleles, qual=svItem.qual, filter=svItem.filter)
                for key in svItem.info.keys():
                    new_record.info[key] = svItem.info[key]
                for index, key in enumerate(svItem.samples.keys()):
                    new_record.samples[f"SAMPLE{index}"]["AF"] = svItem.samples[key]["AF"]
                    af = svItem.samples[key]["AF"]
                    if isinstance(af, tuple):
                        af = af[0]
                    if af == 0:  
                        new_record.samples[f"SAMPLE{index}"]["RNAMES"] = ""
                        if "RHAP" in svItem.samples[key]:
                            new_record.samples[f"SAMPLE{index}"]["RHAP"] = ""
                    else:
                        new_record.samples[f"SAMPLE{index}"]["RNAMES"] = ",".join(svItem.samples[key]["RNAMES"])
                        if "RHAP" in svItem.samples[key]:
                            new_record.samples[f"SAMPLE{index}"]["RHAP"] = svItem.samples[key]["RHAP"]
                new_record.pos = svItem.pos
                new_record.stop = svItem.stop
                new_record.info["SVTYPE"] = svItem.info["SVTYPE"]
                if not (new_record.chrom == svItem.chrom and new_record.pos == svItem.pos and new_record.stop == svItem.stop):
                    raise ValueError("The new record is not the same as the original record", str(svItem), str(new_record))
                outVCF.write(new_record)
                n_out += 1
        outVCF.close()
        vcf1.close()
        vcf2.close()
        print("Output {} from {} SVs to {}".format(n_out, len(result_df), out_vcf))
        tabix_index(out_vcf, force=True, preset="vcf")

    def check_normal_supp(self, svItem):
        if self.ifhighconfi:
            normal_hifi = svItem.samples[list(svItem.samples)[2]]['RHAP']
            normal_ont = svItem.samples[list(svItem.samples)[3]]['RHAP']
            if normal_hifi is None or normal_hifi == ".":
                normal_hifi = 0
            else:
                normal_hifi = len(normal_hifi)
            if normal_ont is None or normal_ont == ".":
                normal_ont = 0
            else:
                normal_ont = len(normal_ont)
            return normal_hifi <= 3 or normal_ont <= 3
        else:
            norm = svItem.samples[list(svItem.samples)[1]]['RNAMES']
            if norm is None or norm == ".":
                norm = 0
            elif isinstance(norm, tuple):
                norm = len(norm)
            else:
                norm = 1
            return norm <= 3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge SVs between two SV vcf files from different ref genomes (using Jaccard similarity and supporting reads.)")
    parser.add_argument("-i", "--corr_txt", type=str, default=None)
    parser.add_argument("-f", "--fileformat", type=str, nargs='+', default=["/home/yuz006/disk/sv_calls/cancer_public/grch38_h2009_merge_15callset_highconfi_lr.vcf.gz"])
    parser.add_argument("-r", "--refs", type=str, nargs='+', default=["grch38", "hap1", "hap2"])
    parser.add_argument("-o", "--output", type=str, default="./corrMerge_mbm.txt")
    parser.add_argument("-t", "--threshold", type=float, default=0.5, help="Threshold for min match score, default is 0.3 (excluded)")
    parser.add_argument("--highconfi", action="store_true", help="If given, use high confidence genotyping across ill, hifi and ont.")
    parser.add_argument("-v", "--verbose", action="store_true", help="If given, show debug messages.")
    args = parser.parse_args()

    # SVCombine(args.corr_txt, args.fileformat, args.refs, args.output, verbose=args.verbose, threshold=args.threshold, ifhighconfi=args.highconfi).run()
    # load corrMerge.txt
    # for c in ["h1437", "h2009", "hcc1937"]:
    #     corr_txt = f"/home/yuz006/disk/t2t_assembly/cancer_public/verkko/hap_only_vali/{c}_corrMax.txt"
    #     corr_matrix = pd.read_csv(corr_txt, sep=",", header=0)
    #     n = 0
    #     yes = 0
    #     matched_svs =  {}
    #     yes_svs = set()
    #     for index, row in corr_matrix.iterrows():
    #         if row["dsa"] != "0":
    #             n += 1
    #             matched_svs[row["dsa"]] = row["hap1"] + "-" + row["hap2"]
    #             if row["hap1"] != "0" or row["hap2"] != "0":
    #                 yes += 1
    #                 yes_svs.add(row["dsa"])
        
    #     from utils import loadDict
    #     region_dict = loadDict(f"/home/yuz006/disk/t2t_assembly/cancer_public/a_{c}_for_igv/{c}bl_hap1.fa.satellite.bed")
    #     region_dict = loadDict(f"/home/yuz006/disk/t2t_assembly/cancer_public/a_{c}_for_igv/{c}bl_hap2.fa.satellite.bed", region_dict)
    

    #     vcf = VariantFile(f"/home/yuz006/disk/t2t_assembly/cancer_public/verkko/hap_only_vali/dsa_{c}_merge_15callset_highconfi_lr.vcf.gz")
    #     severus = 0
    #     yes_severus = 0
    #     for record in vcf:
    #         chrom2, pos2 = get_sv_end(record)
    #         sv_id = str(record.pos) + "%" + str(pos2) + "%" + str(get_sv_len(record)) + "%" + record.id
    #         ifcentr = False
    #         ifsatellite = False
    #         overlap_region = region_dict[record.chrom][record.pos]
    #         overlap_centr = [reg for reg in list(overlap_region) if reg.data == "Satellite/centr"]
    #         overlap_region2 = region_dict[chrom2][pos2]
    #         overlap_centr += [reg for reg in list(overlap_region2) if reg.data == "Satellite/centr"]
    #         if len(overlap_centr) > 0:
    #             ifcentr = True
    #         if len(overlap_region) > 0 or len(overlap_region2) > 0:
    #             ifsatellite = True
    #         if sv_id in matched_svs:
    #             if ifsatellite:
    #                 severus += 1
    #             if ifcentr:
    #                 yes_severus += 1
    #             # if "nanomonsv-hifi" in record.info["Support"] or "pbsv-hifi" in record.info["Support"] or \
    #             #     "severus-hifi" in record.info["Support"] :
    #             #     severus += 1
    #             #     if sv_id in yes_svs:
    #             #         yes_severus += 1
    #             #     print(c, sv_id, record.chrom, record.pos, record.info["SVTYPE"], record.info["SVLEN"] if "SVLEN" in record.info else 0, matched_svs[sv_id], record.info["Support"], ifsatellite, ifcentr)
    #     vcf.close()
    #     print(c, n, yes, len(matched_svs), severus, yes_severus, yes_severus/severus)

    df = pd.read_csv("/home/yuz006/ha.csv", sep=",", header=0)
    
    # traverse each line
    for i, row in df.iterrows():
        vcf = VariantFile(f"/home/yuz006/disk/sv_calls/cancer_public/dsa_{row[0]}_merge_15callset_highconfi_lr.vcf.gz")
        for record in vcf:
            chrom2, pos2 = get_sv_end(record)
            sv_id = str(record.pos) + "%" + str(pos2) + "%" + str(get_sv_len(record)) + "%" + record.id
            if sv_id == row[1]:
                print(record.info["AF"])
        vcf.close()
