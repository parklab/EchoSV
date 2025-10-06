#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2023/11/09
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''
import os
from pysam import VariantFile, tabix_index, FastaFile

def svtype_extract(vcfFile, 
                   svtype=None, 
                   passonly=True, 
                   svlenThre=50, 
                   verbose=True,
                   ):    
    vcf = VariantFile(vcfFile)
    fileName = os.path.splitext(vcfFile)
    if svtype is None:
        outFileName = os.path.splitext(fileName[0])[0]+"_nodup"+ ".vcf.gz"
    else:
        outFileName = fileName[0]+"_"+svtype+fileName[1]
    extractVcf = VariantFile(outFileName, "w", header=vcf.header)
    n_total = 0
    n_select = 0
    for svItem in vcf.fetch():
        n_total += 1
        if passonly and not _ispass(svItem):  # check filter
            continue
        # if not _ishom(svItem):
        #     continue
            
        if svtype is None or ("SVTYPE" in svItem.info and svtype in svItem.info["SVTYPE"].upper()) \
            or ("REPTYPE" in svItem.info and svtype in svItem.info["REPTYPE"].upper()): # check sv type 
            if "SVLEN" in svItem.info:  # check svlen
                svlen = svItem.info["SVLEN"]
                if not isinstance(svlen, int):
                    svlen = svlen[0]
                if abs(svlen) >= svlenThre:
                    extractVcf.write(svItem)    
                    n_select += 1
            elif "gridss" in vcfFile:
                infolist = svItem.alts[0].replace('[', ']').split(']')
                if len(infolist) == 1:
                    svlen = len(infolist[0])
                    if svlen >= svlenThre:
                        svItem.info["SVTYPE"] = "INS"
                        extractVcf.write(svItem)
                        n_select += 1
                else:
                    chrom, region = infolist[1].split(":")
                    if chrom != svItem.chrom or (chrom == svItem.chrom and abs(int(region) - svItem.pos) >= svlenThre):
                        extractVcf.write(svItem)
                        n_select += 1 
            else:            
                extractVcf.write(svItem)
                n_select += 1
    extractVcf.close()
    tabix_index(outFileName, preset="vcf", force=True)
    return outFileName  # saved fileName

def check_svtype(vcfFile, passonly=True):
    svtype = {}
    vcf = VariantFile(vcfFile)
    for svItem in vcf.fetch():
        if passonly and not _ispass(svItem):
            continue
        if not "SVTYPE" in svItem.info:
            raise ValueError("No SVTYPE in VCF file!")
        if svItem.info["SVTYPE"] in svtype:
            svtype[svItem.info["SVTYPE"]] += 1
        else:
            svtype[svItem.info["SVTYPE"]] = 1
    print(f"SVTYPE in {vcfFile}: ", svtype)

def check_svfilter(vcfFile):
    svfilter = {}
    vcf = VariantFile(vcfFile)
    for svItem in vcf.fetch():
        label = list(svItem.filter)[0]
        if label not in svfilter:
            svfilter[label] = 1
        else:
            svfilter[label] += 1
    print(f"SV filter in {vcfFile}: ", svfilter)

def chroms(chr):
    if chr.startswith("chr"):
        chr = chr[3:]
    elif chr.startswith("haplotype"):
        chr = chr.split("-")[1]
    chr = chr.split("_")[0]
    if chr == "X":
        return 23
    if chr == "Y":
        return 24
    if chr.isdigit():
        return int(chr)
    return 25

def removeDup(vcfFile, dupSet, passonly=True):
    vcf = VariantFile(vcfFile)
    contigs_order = list(vcf.header.contigs)
    fileName = os.path.splitext(vcfFile)[0]
    fileName = os.path.splitext(fileName)
    outFileName = fileName[0]+"_nodup"+ ".vcf.gz"
    extractVcf = VariantFile(outFileName, "w", header=vcf.header)
    _n_total = 0

    for svItem in vcf.fetch():
        if (passonly and _ispass(svItem)) or not passonly: 
            if "pbsv" in svItem.id or "sawfish" in svItem.id:
                if svItem.info["SVTYPE"] == "BND":
                    if svItem.id in dupSet:
                        extractVcf.write(svItem)
                        _n_total += 1
                else:
                    extractVcf.write(svItem)
                    _n_total += 1
            elif "MATE_ID" in svItem.info and (svItem.id[-2:] == "_1" or svItem.id[-2:] == "_2") and svItem.id[:-2] in dupSet:  # severus
                chrom1 = svItem.chrom
                pos1 = svItem.pos
                chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
                if chroms(chrom1) < chroms(chrom2) or (chroms(chrom1) == chroms(chrom2) and pos1 <= int(pos2)):
                    # if "SVLEN" in svItem.info and svItem.info["SVLEN"] < 50:
                    #     continue
                    extractVcf.write(svItem)
                    _n_total += 1
            elif "MATEID" in svItem.info and svItem.id[:-2] in dupSet:   # nanomonsv
                if svItem.id[-2:] == "_0":
                    extractVcf.write(svItem)
                    _n_total += 1

            elif svItem.id in dupSet: # duplicated variant
                chrom1 = svItem.chrom
                pos1 = svItem.pos
                chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
                if chroms(chrom1) < chroms(chrom2) or (chroms(chrom1) == chroms(chrom2) and pos1 <= int(pos2)):
                    # if "SVLEN" in svItem.info and svItem.info["SVLEN"] < 50:
                    #     continue
                    extractVcf.write(svItem)
                    _n_total += 1
            else:
                # if "SVLEN" in svItem.info and svItem.info["SVLEN"] < 50:
                #     continue
                extractVcf.write(svItem)
                _n_total += 1
    extractVcf.close()
    tabix_index(outFileName, preset="vcf", force=True)
    return outFileName, _n_total  # saved fileName
    
def _ispass(svItem):
    if ('PASS' in list(svItem.filter) and list(svItem.filter) != []) \
        or list(svItem.filter) == []:
        return True
    else:
        return False
    
def _ishom(svItem):
    gt = svItem.samples[0]["GT"]
    gt_count = (gt[0] if gt[0] is not None else 0) + (gt[1] if gt[1] is not None else 0)    # Calculate Genotype allele count
    return gt_count == 2

def _count(VariantObj, passonly=True, verbose=False):
    n_total = 0
    variant = {}
    for item in VariantObj.fetch():
        _exist = False
        if ((passonly and _ispass(item)) or not passonly) and item.id is not None:
            if "pbsv" in item.id and item.info["SVTYPE"] == "BND":
                pos_list = item.id.split(".")[2].split("-")
                if len(pos_list) == 2:
                    pos1, pos2 = pos_list
                elif len(pos_list) == 4:    # for DSA-result
                    pos1 = pos_list[0] + "-" + pos_list[1]
                    pos2 = pos_list[2] + "-" + pos_list[3]
                
                name1 = "pbsv.BND."+pos1+"-"+pos2
                name2 = "pbsv.BND."+pos2+"-"+pos1
                if name1 not in variant and name2 not in variant:
                    variant[item.id] = 1
                elif name1 in variant:
                    variant[name1] += 1
                    _exist = True
                elif name2 in variant:
                    variant[name2] += 1
                    _exist = True
            elif "sawfish" in item.id and item.info["SVTYPE"] == "BND": # DEL may be reported as BND and show only once
                name1 = item.id
                name2 = item.info["MATEID"][0]
                if name1 not in variant and name2 not in variant:
                    if name1[-1] == "0":
                        variant[name1] = 1
                    else:
                        variant[name2] = 1
                elif name1 in variant:
                    variant[name1] += 1
                    _exist = True
                elif name2 in variant:
                    variant[name2] += 1
                    _exist = True
            elif "MATE_ID" in item.info and (item.id[-2:] == "_1" or item.id[-2:] == "_2"):    # severus, two breakpoint are reported for BND
                # check if item.id ends with "_1" or "_2", remove the suffix
                name = item.id[:-2]
                if name not in variant:
                    variant[name] = 1
                else:
                    variant[name] += 1
                    _exist = True
            elif "MATEID" in item.info and (item.id[-2:] == "_0" or item.id[-2:] == "_1"):    # nanomonsv, two breakpoint are reported for BND
                # check if item.id ends with "_0" or "_1", remove the suffix
                name = item.id[:-2]
                if name not in variant:
                    variant[name] = 1
                else:
                    variant[name] += 1
                    _exist = True
            elif item.id not in variant:
                variant[item.id] = 1
            else:
                variant[item.id] += 1
                _exist = True
            if not _exist:
                n_total += 1
        elif item.id is None:
            n_total += 1
    if verbose:
        return n_total, dict((k, v) for k, v in variant.items() if v > 1)
    return n_total

def regionAnno(svItem, interested_regions):
    if svItem.info["SVTYPE"] == "BND":
        chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
        if bool(interested_regions[svItem.chrom][svItem.pos]) or\
            bool(interested_regions[chrom2][int(pos2)]):
            return True
        else:
            return False
    # if bool(interested_regions[svItem.chrom][svItem.pos: svItem.stop+1]):
    if bool(interested_regions[svItem.chrom][svItem.pos]) or bool(interested_regions[svItem.chrom][svItem.stop]):
        # if not bool(interested_regions[svItem.chrom][svItem.pos: svItem.stop+1]):
            # print(svItem, svItem.pos, svItem.stop, interested_regions[svItem.chrom][svItem.pos: svItem.stop])
        return True
    else:
        return False

def _rowcount(VariantObj, passonly=True, include_tree=None):
    n_total = 0
    for item in VariantObj.fetch():
        if (passonly and _ispass(item)) or not passonly:
            if include_tree is None:
                n_total += 1
            elif regionAnno(item, include_tree):
                n_total += 1
    return n_total

def _checkVCF(vcfFile, passonly=True):
    if not os.path.exists(vcfFile):
        return None, None
    if os.path.getsize(vcfFile) < 30:
        return None, None

    variantObj = VariantFile(vcfFile)
    rawcount = _rowcount(variantObj, passonly=passonly)
    count, dup_variant = _count(variantObj, verbose=True, passonly=passonly)
    variantObj.close()
    # if rawcount == 0:
        # raise ValueError("VCF file is empty!")
    # elif count == 0:
        # return vcfFile, rawcount
    if rawcount != count:
        if "svdss" not in vcfFile:
            # print("VCF file have multiple rows for the same variant! ", rawcount, count, dup_variant)
            newVcfFile, newcount = removeDup(vcfFile, dup_variant, passonly=passonly)
            if newcount != count:
                print("Duplication removal failed!", rawcount, count, newcount, len(dup_variant)) 
                # raise ValueError("Duplication removal failed!", rawcount, count, newcount, len(dup_variant))
            return newVcfFile, newcount
    return vcfFile, rawcount

def _addHeader(vcfFile, new="SVLEN", outFileName=None):
    vcf = VariantFile(vcfFile)
    vcf.header.info.add("SVLEN", ".", "Integer", "Length of the structural variant")
    outvcf = VariantFile(outFileName, "w", header=vcf.header)
    for svItem in vcf.fetch():
        if _ispass(svItem):
            outvcf.write(svItem)
    outvcf.close()
    # tabix_index(outFileName, preset="vcf", force=True)
    return outFileName  # saved fileName

def get_sv_end(svItem, ins_pseudoPos=False):
    if "CHR2" in svItem.info.keys() and "POS2" in svItem.info.keys():
        chrom2 = svItem.info["CHR2"]
        pos2 = svItem.info["POS2"]
    elif "SVTYPE" in svItem.info and (svItem.info['SVTYPE'] == "BND" or svItem.info['SVTYPE'] == "TRANS"):
        if len(svItem.alts) == 1:
            chrom2, pos2 = svItem.alts[0].replace('[', ']').split(']')[1].split(":")
        else:
            raise ValueError("BND format error")
    else:
        chrom2 = svItem.chrom
        pos2 = svItem.stop
        if svItem.info["SVTYPE"] == "INS" and ins_pseudoPos:
            svlen = get_sv_len(svItem)
            pos2 += svlen
    return chrom2, int(pos2)

def get_sv_type(svItem):
    svtype = None
    if "SVTYPE" in svItem.info:
        svtype = svItem.info["SVTYPE"]
    elif "REPTYPE" in svItem.info:
        svtype = svItem.info["REPTYPE"]
    else:
        raise Exception("No SVTYPE or REPTYPE found in vcf file.")
    return svtype

def get_sv_len(svItem):
    svlen = 0
    if "SVLEN" in svItem.info:
        if isinstance(svItem.info["SVLEN"], int) or isinstance(svItem.info["SVLEN"], float):
            svlen = svItem.info["SVLEN"]
        else:
            svlen = svItem.info["SVLEN"][0]
    elif "SVTYPE" in svItem.info and svItem.info['SVTYPE'] == "INS" or "I-Mosaic" in svItem.id:
        if svItem.alts[0] == "<INS>":
            svlen = svItem.info["SVINSLEN"]
        else:
            svlen = len(svItem.alts[0]) - len(svItem.ref)
    elif "SVTYPE" in svItem.info and svItem.info['SVTYPE'] == "BND":
        pass
    else:
        raise ValueError("wrong sv length extraction")
    return int(svlen)