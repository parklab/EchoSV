#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2024/08/18
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu

Given a SV set, extract the reads that support the SVs from given dataset.

Update:
2024/09/01: correct the search range, correct the check of proper_pair; regard INS as DUP
'''

from pysam import VariantFile, tabix_index, AlignmentFile, VariantHeader
import os
from echosv.vcf_utils import _checkVCF, _ispass, _rowcount
from echosv.bed_utils import loadBed
import re, timeit

class AnnSV:
    def __init__(self, sv_vcf=None, 
                 bam_list=[], 
                 bed_file=None,
                 out_file=None,
                 min_mapq=1,
                 dist_threshold=500, 
                 ifWriteReadNames=True,
                 ratio_of_seqsize=0.5,
                 verbose=False):
        self.sv_vcf = sv_vcf
        self.bam_list = bam_list
        self.bam = None
        self.sampleList = []
        self.bed_file = bed_file
        self.interested_regions = None
        self.out_file = out_file
        self.min_mapq = min_mapq
        self.dist_threshold = dist_threshold
        self.ifWriteReadNames = ifWriteReadNames
        self.ratio_of_seqsize = ratio_of_seqsize
        self.verbose = verbose

    def _checkFile(self):
        if not os.path.exists(self.sv_vcf):
            raise FileNotFoundError("SV vcf file not found.")
        else:
            self.sv_vcf, n_sv = _checkVCF(self.sv_vcf)
        for bam in self.bam_list:
            if not os.path.exists(bam):
                raise FileNotFoundError("Bam file not found.")
            # Extract the file name of self.bam_file
            self.sampleList.append(os.path.basename(bam).replace(".bam", ""))
        
        if (not self.bed_file is None) and os.path.exists(self.bed_file):
            self.interested_regions = loadBed(self.bed_file)

        if self.out_file is None:
            if "nodup" in self.sv_vcf:
                self.out_file = self.sv_vcf.replace("_nodup.vcf", "_reads.vcf")
            else:
                self.out_file = self.sv_vcf.replace(".vcf", "_reads.vcf")
        return n_sv

    def run(self):
        n_sv = self._checkFile()
        vcf = VariantFile(self.sv_vcf)
        # copy all header records except for the sample information
        header = VariantHeader()
        header.add_line("##params=MinMapQ={},DistThreshold={},WriteReadNames={},RatioOfSeqSize={}".format(
            self.min_mapq, self.dist_threshold, self.ifWriteReadNames, self.ratio_of_seqsize))
        for record in vcf.header.records:
            header.add_record(record)
        for sample in self.sampleList:
            header.add_sample(sample)  # add sample for each BAM input
        if not "AF" in vcf.header.formats:
            header.formats.add("AF", number=1, type="Float", description="Allele frequency of the variant")
        if not "RNAMES" in vcf.header.formats:
            header.formats.add("RNAMES", number=1, type="String", description="Read names supporting the variant")
        if self.bed_file is not None:
            header.info.add("Region", number=0, type="Flag", description="Region annotation of the variant")
        annoVcf = VariantFile(self.out_file, "w", header=header)
        n_region_sv = 0
        n_sv = 0
        for svItem in vcf.fetch():
            if not _ispass(svItem):
                continue
            record = annoVcf.new_record(contig=svItem.contig, start=svItem.start, alleles=svItem.alleles,
                                        qual=svItem.qual, filter=svItem.filter.keys(), id=svItem.id, info=svItem.info)
            record.stop = svItem.stop
            if self.bed_file is not None:   
                if self.regionAnno(svItem):
                    n_region_sv += 1
                    record.info["Region"] = True
            #for sample in svItem.samples.keys():
                #record.samples[sample]['AF'] = svItem.samples[sample]['AF']
                #record.samples[sample]['RNAMES'] = svItem.samples[sample]['RNAMES']
            # add sample information from each BAM
            for sample_index in range(len(self.sampleList)):
                sampleName = self.sampleList[sample_index]
                self.bam = AlignmentFile(self.bam_list[sample_index], "rb")
                supp_reads, af = self.get_sv_reads(svItem)
                self.bam.close()
            
                record.samples[sampleName]['AF'] = af
                if self.ifWriteReadNames:
                    supp_reads = [read.replace(":", "_") for read in supp_reads]
                    record.samples[sampleName]['RNAMES'] = ",".join(supp_reads)
                else:
                    record.samples[sampleName]['RNAMES'] = str(len(supp_reads))

            annoVcf.write(record)
            if self.verbose:
                print(record)
            n_sv += 1
            if n_sv % 1000 == 0:
                print("Processed {} SVs.".format(n_sv))
            #if abs(svItem.info['AF'] - record.samples[sampleName]['AF']) > 0.1:
            #    print(svItem, svItem.info['AF'], af, len(supp_reads), cov)
            #    n_region_sv += 1
        print("Number of SVs in interested regions: {} out of {}".format(n_region_sv, n_sv))
        annoVcf.close()
        tabix_index(self.out_file, preset="vcf", force=True)

    def get_sv_type(self, svItem):
        svtype = None
        if "SVTYPE" in svItem.info:
            svtype = svItem.info["SVTYPE"]
        elif "REPTYPE" in svItem.info:
            svtype = svItem.info["REPTYPE"]
        else:
            raise Exception("No SVTYPE or REPTYPE found in vcf file.")
        return svtype

    def regionAnno(self, svItem):
        if "BND" in self.get_sv_type(svItem):
            chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
            if bool(self.interested_regions[svItem.chrom][svItem.pos]) and \
                bool(self.interested_regions[chrom2][int(pos2)]):
                return True
            else:
                return False
        if bool(self.interested_regions[svItem.chrom][svItem.pos: svItem.stop]):
            return True
        else:
            return False

    def get_sv_reads(self, svItem): # collect variant supporting reads: clipped and discordant reads
        svtype = self.get_sv_type(svItem)
        chrom1 = svItem.chrom
        pos1 = svItem.pos
        pos1_end = pos1 + 1 + self.dist_threshold
        if "BND" in svtype:
            chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
            pos2 = int(pos2)
        else:
            chrom2 = svItem.chrom
            pos2 = svItem.stop
            if "INS" in svtype:
                pos2 = svItem.pos + self._getSVlen(svItem)
            mid_pos = (pos1 + pos2) // 2
        supp_reads = set()  # hash table
        unsupp_reads = set()

        left_aln_pos = {}
        # if self.verbose:
        #     print(chrom1, pos1, chrom2, pos2, svtype)
        n_left_read = 0
        n_right_read = 0
        improper_aln = 0
        clipped_aln = 0
        cigar_reads = 0
        n_unsupp_cigar = 0
        n_unsupp_pair = 0
        n_unsupp_right = 0
        n_unsupp_left = 0
        if "BND" in svtype:
            end_pos = pos1 + self.dist_threshold
        else:
            end_pos = min(mid_pos, pos1 + self.dist_threshold)
        if self.verbose:
            print("left region", chrom1, max(0, pos1-self.dist_threshold), end_pos)
        for read in self.bam.fetch(chrom1, max(0, pos1-self.dist_threshold), end_pos):
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality >= self.min_mapq:
                if read.query_name not in left_aln_pos:
                    left_aln_pos[read.query_name] = {}
                # record each read alignment located near the left breakpoint of SV
                left_aln_pos[read.query_name][f"{read.reference_start}-{read.reference_end}"] = read
        # print("left aln pos", len(left_aln_pos), left_aln_pos)
        if "BND" in svtype:
            start_pos = max(0, pos2 - self.dist_threshold)
        else:
            start_pos = max(mid_pos, pos2 - self.dist_threshold)
        if self.verbose:
            print("right region", chrom2, start_pos, pos2+1+self.dist_threshold)
        for read in self.bam.fetch(chrom2, start_pos, pos2+1+self.dist_threshold):
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality >= self.min_mapq and read.query_name not in supp_reads:
                isSupport = None
                mapper = read.get_aligned_pairs(matches_only=True)  # (query_pos, ref_pos)
                if read.query_name in left_aln_pos: # a same read
                    stored_aln = left_aln_pos[read.query_name]
                    if f"{read.reference_start}-{read.reference_end}" in stored_aln:   # a same segment of a read
                        if self._getSVlen(svItem) < 150:
                            isSupport = self._checkSVinAln(read, svtype, pos1, pos2, self._getSVlen(svItem), mapper)
                            if isSupport:
                                supp_reads.add(read.query_name)
                                cigar_reads += 1
                                if self.verbose:
                                    print("within", read.query_name, isSupport, read.reference_start, [item for item in left_aln_pos[read.query_name]])

                            elif isSupport == False:
                                unsupp_reads.add(read.query_name)
                                n_unsupp_cigar += 1
                            if isSupport != None:
                                left_aln_pos[read.query_name].pop(f"{read.reference_start}-{read.reference_end}")
                                if len(left_aln_pos[read.query_name]) == 0:
                                    left_aln_pos.pop(read.query_name)
                    else:
                        for posStr, aln in stored_aln.items():
                            if read.is_read1 == aln.is_read1:
                                isSupport = True
                                clipped_aln += 1
                                if self._getSVlen(svItem) > 300:    ######
                                    supp_reads.add(read.query_name)
                                break
                        if not isSupport:
                            if not read.is_proper_pair:
                                improper_aln += 1
                                supp_reads.add(read.query_name)   #####
                            else:
                                isSupport = self._checkBND(read, read.reference_start, read.reference_end-1, [pos2])
                                if isSupport:
                                    supp_reads.add(read.query_name)
                                    n_right_read += 1
                                elif isSupport == False and read.reference_start < pos2 < read.reference_end:
                                    unsupp_reads.add(read.query_name)
                                    n_unsupp_pair += 1
                                # check the left breakpoint
                                if read.query_name not in supp_reads:
                                    for posStr, aln in stored_aln.items():
                                        isSupport = self._checkBND(aln, aln.reference_start, aln.reference_end-1, [pos1])
                                        if isSupport:
                                            supp_reads.add(read.query_name)
                                            n_left_read += 1
                                        elif isSupport == False and aln.reference_start < pos1 < aln.reference_end:
                                            unsupp_reads.add(read.query_name)
                                            n_unsupp_left += 1
                                    # print("right unsupp", read.query_name, read.reference_start, read.reference_end, pos2)

                else:
                    isSupport = self._checkBND(read, read.reference_start, read.reference_end-1, [pos2])
                    # if isSupport:   # the read has matched breakpoint with soft-clip, target variant supporting
                    #     supp_reads.add(read.query_name)
                    #     n_right_read += 1
                    if isSupport == False and read.reference_start < pos2 < read.reference_end and read.is_proper_pair:    # ref-supporting read
                        unsupp_reads.add(read.query_name)
                        n_unsupp_right += 1
                #     else:
                #         isSupport = None
        for read_name, pos_dict in left_aln_pos.items():
            if read_name not in supp_reads:
                for posStr, aln in pos_dict.items():
            #         # print("left soft-clipped check")
                    isSupport = self._checkBND(aln, aln.reference_start, aln.reference_end-1, [pos1])
                    # if isSupport:   # the read has matched breakpoint with soft-clip
                    #     supp_reads.add(read_name)
                    #     n_left_read += 1
                    if isSupport == False and aln.reference_start < pos1 < aln.reference_end and read.is_proper_pair: # ref-supporting read
                        unsupp_reads.add(read_name)
                        n_unsupp_left += 1
        #         else:
        #             isSupport = None
        #         if self.verbose:
        #             print("left", read_name, posStr, isSupport)
        # if supp_reads.intersection(unsupp_reads):
        #     print(svItem)
        #     print(supp_reads.intersection(unsupp_reads))
            # raise Exception("Reads in both supported and unsupported sets.")
        unsupp_reads = unsupp_reads - supp_reads
        if self.verbose:
            print(svItem)
            print("total supp: {}, cigar reads: {}, cross improper reads: {}, cross clipped reads: {}, left reads {} right reads {}".format(
                len(supp_reads), cigar_reads, improper_aln, clipped_aln, n_left_read, n_right_read))
            print("total unsupp: {}, unsupp cigar reads: {}, unsupp pair reads: {}, unsupp left reads {}, unsupp right reads {}".format(
                len(unsupp_reads), n_unsupp_cigar, n_unsupp_pair, n_unsupp_left, n_unsupp_right))
        # nominator = len(supp_reads) - (n_left_read + n_right_read)*0.5
        nominator = len(supp_reads)
        cov = nominator + len(unsupp_reads) - (n_unsupp_left + n_unsupp_right)*0.5
        af = nominator / cov if cov != 0 else 0
        return supp_reads, af
        # return supp_reads, len(supp_reads) + len(unsupp_reads) - (n_left_unsupp_reads + n_right_unsupp_reads)*0.5
    
    def _getSVlen(self, svItem):
        if "SVLEN" in svItem.info:
            svlen = svItem.info["SVLEN"]
            if isinstance(svlen, int):
                return abs(svlen)
            if isinstance(svlen, tuple):
                return abs(svlen[0])
            if isinstance(svlen, float):
                return abs(int(svlen))
            # print(svItem)
            raise Exception("No valid SVLEN found in vcf file.")
        else:
            if len(svItem.alts) > 1:
                # print(svItem)
                raise Exception("Can't compute valid SVLEN in vcf file.")
            return abs(len(svItem.alts[0]) - len(svItem.ref))
    
    def _checkSVinAln(self, aln, svtype, ref_pos1, ref_pos2, svlen=None, mapper=None):
        q_pos1 = q_pos2 = r_pos1 = r_pos2 = None
        for pair in mapper:
            if q_pos1 is None and ref_pos1 - self.dist_threshold <= pair[1] <= ref_pos1 + self.dist_threshold:
                q_pos1 = pair[0]
                r_pos1 = pair[1]
                break
        for pair in mapper[::-1]:
            if q_pos2 is None and ref_pos2 - self.dist_threshold <= pair[1] <= ref_pos2 + self.dist_threshold:
                q_pos2 = pair[0]
                r_pos2 = pair[1]
                break
        if r_pos1 is None or r_pos2 is None:
            return None
        q_dist = abs(q_pos2 - q_pos1)   # query distance
        r_dist = abs(r_pos2 - r_pos1)   # read distance
        if self.verbose:
            print(aln.query_name, svlen, q_pos1, q_pos2, r_pos1, r_pos2, svtype, q_dist, r_dist)
        if q_dist < r_dist - 30 and svtype == "DEL" and abs(abs(q_dist - r_dist) - svlen) <= min(self.dist_threshold, svlen * self.ratio_of_seqsize):
            # if self.verbose:
            #     print("DEL", aln.query_name, q_dist, r_dist, svlen, q_pos1, q_pos2, r_pos1, r_pos2)
        
            return True
        if q_dist > r_dist + 30 and svtype in ["INS", "DUP"] and abs(abs(q_dist - r_dist) - svlen) <= min(self.dist_threshold, svlen * self.ratio_of_seqsize):
            # if self.verbose:
            #     print("INS", aln.query_name, q_dist, r_dist, svlen, q_pos1, q_pos2, r_pos1, r_pos2)
        
            return True
        # print("within one segment of one read")
        isSupport = self._checkBND(aln, mapper[0][1], mapper[-1][1], [ref_pos1, ref_pos2])
        if isSupport == False and (r_pos1 > ref_pos2 or r_pos2 < ref_pos1): # ref-supporting read
            return None
        return isSupport
    
    def _checkBND(self, read, left_end, right_end, posList, mismatch_thre=50):
        '''
        Check if the read alignment support a breakpoint of SV.
        If its end is solf-clipped around the breakpoint and has no other supp alignment, added into supp_reads,
        elif the aln support the reference sequence of the breakpoint, added into unsupp_reads.
        else, discard the aln.

        TO DO: soft-clipped part may align to different genomic part rather than SV infered region. 
        Overestimate the supp reads.
        '''
        cigarstring = read.cigarstring
        soft_clip_start_pattern = re.compile(r'^[0-9]+[SH]')
        soft_clip_end_pattern = re.compile(r'[0-9]+[SH]$')
        clp_pos_list = []
        start_match = soft_clip_start_pattern.search(cigarstring)
        if bool(start_match) and int(start_match.group(0)[:-1]) >= mismatch_thre:
            clp_pos_list.append(left_end)
        end_match = soft_clip_end_pattern.search(cigarstring)
        if bool(end_match) and int(end_match.group(0)[:-1]) >= mismatch_thre:
            clp_pos_list.append(right_end)
        if len(clp_pos_list) == 0:
            return False
        for pos in clp_pos_list:
            for ref_pos in posList:
                # if abs(pos - ref_pos) <= self.dist_threshold:
                if abs(pos - ref_pos) <= mismatch_thre:
                    return True
        
        