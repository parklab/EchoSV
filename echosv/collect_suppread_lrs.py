#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2024/01/25
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu

Given a SV set, extract the reads that support the SVs, AF and region annotation from given dataset.
Update:
2024/08/30: Regard any INS as DUP, to correct the potential mistake of SVTYPE in vcf file. check the SuppReads across both breakpoints.
2024/10/12: Speed up the script by multi-processing of each contig.
2024/10/14: Correct the INS bug: reads that only span the left breakpoint of INS (not DUP) will be omited.
2025/01/25: Add pbsv-VNTR regions for extensive searching of supporting reads.
TO DO:
Add orientation for SVs, consider the orientation of SuppReads with blunt end (>=50bp soft-clipped). 
'''

from pysam import VariantFile, tabix_index, AlignmentFile, VariantHeader
import os
from echosv.vcf_utils import _checkVCF, _ispass, _rowcount
from echosv.bed_utils import loadBed
import re, warnings, timeit
warnings.filterwarnings("ignore")
import multiprocessing, timeit


def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def get_hp_tag(read):
    if read.has_tag("HP"):
        return str(read.get_tag("HP"))
    else:
        return "0"

class AnnSV:
    def __init__(self, sv_vcf=None, 
                 bam_list=[], 
                 bed_file=None,
                 out_file=None,
                 min_mapq=1,
                 dist_threshold=500, 
                 ifWriteReadNames=True,
                 ifINSasDUP=True,
                 svs_per_job=100,
                 ratio_of_seqsize=0.5,
                 vntr=None,
                 verbose=False):
        self.sv_vcf = sv_vcf
        self.bam_list = bam_list
        self.sampleList = []
        self.bed_file = bed_file
        self.interested_regions = None
        self.out_file = out_file
        self.min_mapq = min_mapq
        self.dist_threshold = dist_threshold
        self.ifWriteReadNames = ifWriteReadNames
        self.ifINSasDUP = ifINSasDUP
        self.svs_per_job = svs_per_job
        self.ratio_of_seqsize = ratio_of_seqsize
        self.vntr = loadBed(vntr, max_size=10000) if vntr else None
        self.verbose = verbose

    def _checkFile(self):
        if not os.path.exists(self.sv_vcf):
            raise FileNotFoundError("SV vcf file not found.")
        elif "merge" in self.sv_vcf or "smaht" in self.sv_vcf or "truthset" in self.sv_vcf:
            n_sv = _rowcount(VariantFile(self.sv_vcf))
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
            self.out_file = self.out_file.replace("_reads.vcf", "_reads_lr.vcf")
        return n_sv

    def run(self):
        n_sv = self._checkFile()
        vcf = VariantFile(self.sv_vcf)
        svs_info_per_run = []
        n_jobs = n_sv // self.svs_per_job + 1
        for i in range(n_jobs):
            svs_info_per_run.append({})
        n_index = 0
        for svItem in vcf.fetch():
            if not _ispass(svItem):
                continue
            n_job_index = n_index // self.svs_per_job
            sv_id = "{}%{}%{}%{}%{}".format(svItem.chrom, svItem.pos, svItem.stop, svItem.id, svItem.alts[0])
            svs_info_per_run[n_job_index][sv_id] = [self.get_sv_type(svItem), self._getSVlen(svItem), svItem.alts]  # svtype, svlen, alts
            n_index += 1
        vcf.close()

        # process each bin of SVs
        pool = multiprocessing.Pool(processes=16)
        svs_anno = pool.map(self._processBin, svs_info_per_run)
        pool.close()
        pool.join()
        # merge a list of dict into one dict
        svs_anno = {k: v for d in svs_anno for k, v in d.items()}

        vcf = VariantFile(self.sv_vcf)
        # copy all header records except for the sample information
        header = VariantHeader()
        # add values of params into header, i.e. the number of min_mapping_quality, dist_threshold
        header.add_line("##params=MinMapQ={},DistThreshold={},WriteReadNames={},INSasDUP={},RatioOfSeqSize={}".format(
            self.min_mapq, self.dist_threshold, self.ifWriteReadNames, self.ifINSasDUP, self.ratio_of_seqsize))
        for record in vcf.header.records:
            header.add_record(record)
        for sample in self.sampleList:
            header.add_sample(sample)  # add sample for each BAM input
        if not "AF" in vcf.header.formats:
            header.formats.add("AF", number=1, type="Float", description="Allele frequency of the variant")
        if not "RNAMES" in vcf.header.formats:
            header.formats.add("RNAMES", number=1, type="String", description="Read names supporting the variant")
        if not "RHAP" in vcf.header.formats:
            header.formats.add("RHAP", number=1, type="String", description="Haplotag of reads supporting the variant")
        if self.bed_file is not None:
            header.info.add("Region", number=0, type="Flag", description="Region annotation of the variant")
        
        n_region_sv = 0
        n_sv = 0
        annoVcf = VariantFile(self.out_file, "w", header=header)
        for svItem in vcf.fetch():
            record = annoVcf.new_record(contig=svItem.contig, start=svItem.start, alleles=svItem.alleles,
                                            qual=svItem.qual, filter=svItem.filter.keys(), id=svItem.id, info=svItem.info)
            record.stop = svItem.stop
            sv_id = "{}%{}%{}%{}%{}".format(svItem.chrom, svItem.pos, svItem.stop, svItem.id, svItem.alts[0])
            for sample_index in range(len(self.sampleList)):
                sampleName = self.sampleList[sample_index]
                af, supp_info, haps_string = svs_anno[sv_id][sample_index]
                
                record.samples[sampleName]['AF'] = af
                record.samples[sampleName]['RNAMES'] = supp_info
                record.samples[sampleName]['RHAP'] = haps_string
            if self.bed_file is not None:   
                if self.regionAnno(svItem):
                    record.info["Region"] = True
                    n_region_sv += 1
            annoVcf.write(record)
            n_sv += 1
        annoVcf.close()
        tabix_index(self.out_file, preset="vcf", force=True)
        vcf.close()
        print("Number of SVs in interested regions: {} out of {}".format(n_region_sv, n_sv))

    def _processBin(self, svBin):    # dict
        start = timeit.default_timer()
        svs_anno = {}
        first_sv_id = None
        for sv_id, info in svBin.items():
            svs_anno[sv_id] = []    # list of (af, supp_reads, ifRegion)
            if first_sv_id is None:
                first_sv_id = sv_id
        for sample_index in range(len(self.sampleList)):    
            bam = AlignmentFile(self.bam_list[sample_index], "rb")
            for sv_id, info in svBin.items():
                supp_reads, cov, supp_read_haps = self.get_sv_reads((sv_id, info), bam)
                af = len(supp_reads) / cov if cov != 0 else 0
                if self.ifWriteReadNames:
                    supp_info = ",".join(supp_reads)
                else:
                    supp_info = str(len(supp_reads))
                svs_anno[sv_id].append((af, supp_info, supp_read_haps))
            bam.close()
        if self.verbose:
            print("Processed SVs from {} in {}s.".format(first_sv_id, timeit.default_timer()-start))
        return svs_anno

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

    def get_sv_reads(self, svItem, bam=None):
        sv_chrom, sv_pos, sv_end, sv_id, alt = svItem[0].split("%")
        sv_type, sv_len, sv_alts = svItem[1]
        sv_pos = int(sv_pos)
        sv_end = int(sv_end)
        chrom1 = sv_chrom
        pos1 = sv_pos
        if "BND" in sv_type:
            chrom2, pos2 = sv_alts[0].replace('[', ']').split(']')[1].split(":")
            pos2 = int(pos2)
        elif "INS" in sv_type:
            chrom2 = sv_chrom
            pos2 = sv_pos + 1
            if self.ifINSasDUP:
                pos2 = sv_pos + sv_len
        elif sv_end == sv_pos:
            if "[" in sv_alts[0] or "]" in sv_alts[0]:
                chrom2, pos2 = sv_alts[0].replace('[', ']').split(']')[1].split(":")
                pos2 = int(pos2)
            else:
                pos2 = sv_pos + sv_len
                chrom2 = chrom1
        else:
            chrom2 = sv_chrom
            pos2 = sv_end
        supp_reads = set()  # hash table
        supp_read_haps = ""
        unsupp_reads = set()

        left_aln_pos = {}
        n_left_unsupp_reads = 0
        n_right_unsupp_reads = 0
        # if self.verbose:
        #     print(chrom1, pos1, chrom2, pos2, svtype)
        n_cross_read = 0
        n_left_read = 0
        n_right_read = 0
        # introduce vntr checking
        left_start = max(0, pos1-self.dist_threshold)
        left_end = pos1 + 1 + self.dist_threshold
        if self.vntr is not None:
            vntr_region = self.vntr[sv_chrom][sv_pos]
            if len(vntr_region) == 1:
                # print(vntr_region, list(vntr_region)[0].begin, list(vntr_region)[0].end)
                left_start = min(list(vntr_region)[0].begin, left_start)
                left_end = max(list(vntr_region)[0].end, left_end)
            
        for read in bam.fetch(chrom1, left_start, left_end):
            if read.mapping_quality >= self.min_mapq:
                if read.query_name not in left_aln_pos:
                    left_aln_pos[read.query_name] = {}
                # record each read alignment located near the left breakpoint of SV
                left_aln_pos[read.query_name][f"{read.reference_start}-{read.reference_end}"] = read
        # introduce vntr checking
        right_start = max(0, pos2-self.dist_threshold)
        right_end = pos2 + 1 + self.dist_threshold
        if self.vntr is not None:
            vntr_region = self.vntr[chrom2][pos2]
            if len(vntr_region) == 1:
                # print(vntr_region, list(vntr_region)[0].begin, list(vntr_region)[0].end)
                right_start = min(list(vntr_region)[0].begin, right_start)
                right_end = max(list(vntr_region)[0].end, right_end)
        for read in bam.fetch(chrom2, right_start, right_end):
            if read.mapping_quality >= self.min_mapq and read.query_name not in supp_reads:
                isSupport = None
                mapper = read.get_aligned_pairs(matches_only=True)  # (query_pos, ref_pos)
                if read.query_name in left_aln_pos: # a same read
                    if f"{read.reference_start}-{read.reference_end}" in left_aln_pos[read.query_name] and len(left_aln_pos[read.query_name].keys()) == 1:   # a same segment of a read
                        isSupport = self._checkSVinAln(read, sv_type, pos1, pos2, sv_len, mapper, (left_start, left_end, right_start, right_end))
                        if isSupport:
                            supp_reads.add(read.query_name)
                            supp_read_haps += get_hp_tag(read)
                        elif isSupport == False:
                            unsupp_reads.add(read.query_name)
                        if self.verbose:
                            print("span within", read.query_name, isSupport, read.reference_start, [item for item in left_aln_pos[read.query_name]])
                        if isSupport != None:
                            left_aln_pos[read.query_name].pop(f"{read.reference_start}-{read.reference_end}")
                            if len(left_aln_pos[read.query_name]) == 0:
                                left_aln_pos.pop(read.query_name)
                    else:
                        isSupport = True
                        supp_reads.add(read.query_name)
                        supp_read_haps += get_hp_tag(read)
                        if self.verbose:
                            print("span", read.query_name, isSupport, read.reference_start, [item for item in left_aln_pos[read.query_name]], len(supp_reads))
                        left_aln_pos.pop(read.query_name)
                        n_cross_read += 1
                else:
                    # print("right soft-clipped check")
                    isSupport = self._checkBND(read, read.reference_start, read.reference_end-1, [pos2])
                    if isSupport:   # the read has matched breakpoint with soft-clip, target variant supporting
                        supp_reads.add(read.query_name)
                        supp_read_haps += get_hp_tag(read)
                        n_right_read += 1
                    elif isSupport == False and read.reference_start < pos2 < read.reference_end:    # ref-supporting read
                        unsupp_reads.add(read.query_name)
                        n_right_unsupp_reads += 1
                    else:
                        isSupport = None
                    # print("right", read.query_name, isSupport)
        for read_name, pos_dict in left_aln_pos.items():
            for posStr, aln in pos_dict.items():
                # print("left soft-clipped check")
                isSupport = None
                if "INS" in sv_type:
                    isSupport = self._checkSVinAln(aln, sv_type, pos1, pos2-sv_len, sv_len, aln.get_aligned_pairs(matches_only=True), (left_start, left_end, left_start, left_end))
                if not isSupport:
                    isSupport = self._checkBND(aln, aln.reference_start, aln.reference_end-1, [pos1])
                if isSupport:   # the read has matched breakpoint with soft-clip
                    supp_reads.add(read_name)
                    supp_read_haps += get_hp_tag(aln)
                    n_left_read += 1
                elif isSupport == False and aln.reference_start < pos1 < aln.reference_end: # ref-supporting read
                    unsupp_reads.add(read_name)
                    n_left_unsupp_reads += 1
                else:
                    isSupport = None
                # print("left", read_name, isSupport)
        if supp_reads.intersection(unsupp_reads):
            print(svItem)
            print(supp_reads.intersection(unsupp_reads))
            # raise Exception("Reads in both supported and unsupported sets.")
        if self.verbose:
            print(svItem, len(supp_reads), supp_reads)
            print("cross read: {}, left reads {} right reads {}".format(n_cross_read, n_left_read, n_right_read))
        return supp_reads, len(supp_reads)+len(unsupp_reads)-int(0.5*(n_left_unsupp_reads+n_right_unsupp_reads)), supp_read_haps
        # return supp_reads, len(supp_reads) + len(unsupp_reads)
    
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
    
    def _checkSVinAln(self, aln, svtype, ref_pos1, ref_pos2, svlen=None, mapper=None, search_region=(0, 0, 0, 0)):
        q_pos1 = q_pos2 = r_pos1 = r_pos2 = None
        left_start, left_end, right_start, right_end = search_region
        for pair in mapper:
            if q_pos1 is None and left_start <= pair[1] <= left_end:
                q_pos1 = pair[0]
                r_pos1 = pair[1]
                break
        for pair in mapper[::-1]:
            if q_pos2 is None and right_start <= pair[1] <= right_end:
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
            if self.verbose:
                print("DEL", aln.query_name, q_dist, r_dist, svlen, q_pos1, q_pos2, r_pos1, r_pos2)
        
            return True
        if q_dist > r_dist + 30 and svtype in ["INS", "DUP"] and abs(abs(q_dist - r_dist) - svlen) <= min(self.dist_threshold, svlen * self.ratio_of_seqsize):
            if self.verbose:
                print("INS", aln.query_name, q_dist, r_dist, svlen, q_pos1, q_pos2, r_pos1, r_pos2)
        
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
        mismatch_thre = min(mismatch_thre, self.dist_threshold)
        isSupport = False
        cigarstring = read.cigarstring
        soft_clip_start_pattern = re.compile(r'^[0-9]+S')
        soft_clip_end_pattern = re.compile(r'[0-9]+S$')
        clp_pos_list = []
        start_match = soft_clip_start_pattern.search(cigarstring)
        if bool(start_match) and int(start_match.group(0)[:-1]) >= mismatch_thre:
            clp_pos_list.append(left_end)
        end_match = soft_clip_end_pattern.search(cigarstring)
        if bool(end_match) and int(end_match.group(0)[:-1]) >= mismatch_thre:
            clp_pos_list.append(right_end)
        for pos in clp_pos_list:
            for ref_pos in posList:
                # if abs(pos - ref_pos) <= self.dist_threshold:
                if abs(pos - ref_pos) <= mismatch_thre:
                    if not read.has_tag('SA'):  # supp alignment could also have SA tag
                        return True
                    supp_segments = read.get_tag('SA').split(";")
                    valid_supp_segments = [x for x in supp_segments if len(x.split(","))>1 and int(x.split(",")[4]) > self.min_mapq]
                    if len(valid_supp_segments) == 0:
                        return True
                    # if self.verbose:
                    #     print("different alternative allele", read.reference_name, clp_pos_list, posList, supp_segments, valid_supp_segments, read.query_name)
                    
                    isSupport = None
        return isSupport
    
def mergeVcf(vcf_files, out_file):
    '''
    Merge multiple vcf files into one.
    '''
    vcf_list = [VariantFile(vcf_file) for vcf_file in vcf_files]
    samples = [vcf.header.samples[0] for vcf in vcf_list]
    header = VariantHeader()
    for record in vcf_list[0].header.records:
        header.add_record(record)
    for sample in samples:
        header.add_sample(sample)
    outvcf = VariantFile(out_file, "w", header=header)
    print("Merging {} vcf files into one.".format(len(vcf_files)))
    for svItem in vcf_list[0]:
        record = outvcf.new_record(contig=svItem.contig, start=svItem.start, alleles=svItem.alleles,
                                        qual=svItem.qual, filter=svItem.filter.keys(), id=svItem.id, info=svItem.info)
        record.stop = svItem.stop
        record.samples[samples[0]]['AF'] = svItem.samples[samples[0]]['AF']
        record.samples[samples[0]]['RNAMES'] = svItem.samples[samples[0]]['RNAMES']
        for i in range(1, len(vcf_list)):
            vcf = vcf_list[i]
            sample = samples[i]
            # print(sample, vcf)
            svs = vcf.fetch(svItem.contig, svItem.start, svItem.stop)
            for svv in svs:
                if svv.id == svItem.id:
                    # print(svv)
                    record.samples[sample]['AF'] = svv.samples[sample]['AF']
                    record.samples[sample]['RNAMES'] = svv.samples[sample]['RNAMES']
        outvcf.write(record)
    outvcf.close()
    tabix_index(out_file, preset="vcf", force=True)

def run(params):
    svs, contig = params
    result = []
    for sv_id, info in svs.items():
        if sv_id.split("%")[0] == contig:
            result.append(info)
    return result
    
