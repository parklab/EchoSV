#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/01/22
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import gzip, timeit, argparse, os
from Bio import SeqIO
from pysam import AlignmentFile

def extract_contig_lengths(fasta_file, write_dict=False, verbose=True):
    '''
    Extract contig lengths from a fasta file.
    '''
    contig_lengths = {}
    if os.path.exists(fasta_file.replace(".fa", ".fa.dict")):
        dict_file = fasta_file.replace(".fa", ".fa.dict")
        with open(dict_file, "r") as f:
            for line in f:
                info = line.strip().split("\t")
                if info[0] == "@SQ":
                    contig_lengths[info[1].split(":")[1]] = int(info[2].split(":")[1])
    elif os.path.exists(fasta_file.replace(".fa", ".dict")):
        dict_file = fasta_file.replace(".fa", ".dict")
        with open(dict_file, "r") as f:
            for line in f:
                info = line.strip().split("\t")
                contig_lengths[info[0]] = int(info[1])
    else:
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contig_lengths[record.id] = len(record.seq)

    import json
    if write_dict:
        # Save the dictionary to a JSON file
        dict_file = fasta_file.replace(".fa", ".dict")
        with open(dict_file, "w") as f:
            for contig, length in contig_lengths.items():
                f.write("{}\t{}\n".format(contig, length))
        print(f"Contig lengths saved to {dict_file}")
    else:
        total_length = sum(contig_lengths.values())
        if verbose:
            print(f"{len(contig_lengths)} contigs and {total_length} bps assemblied in fasta {fasta_file}")
    return contig_lengths

def get_coverage_bam(bam_file, output_file, mapq_threshold=1):
    '''
    Given a bam file, find the regions with uniquely mappable coverage.
    '''
    alignment_file = AlignmentFile(bam_file, "rb")
    chroms = {} # chroms of reference genome
    for chrom in alignment_file.references:
        chroms[chrom] = alignment_file.get_reference_length(chrom)
    def cal_coverage(chrom):
        dsa_map_region = {}
        coverage = [0] * chroms[chrom]
        for read in alignment_file.fetch(chrom):
            if read.mapping_quality < mapq_threshold:
                continue
            start = read.reference_start
            end = read.reference_end
            if end > chroms[chrom]:
                print(f"Warning: {read.query_name} has alignment {start}-{end} out of chrom {chrom} range {chroms[chrom]}.")
            # for i in range(start, min(end, chroms[chrom])):
            #     coverage[i] += 1
            if read.query_name not in dsa_map_region:
                dsa_map_region[read.query_name] = []
            prefix_len = 0
            if read.is_reverse and read.cigartuples[-1][0] == 5:    # hard clip, 5: H
                prefix_len = read.cigartuples[-1][1]
            elif not read.is_reverse and read.cigartuples[0][0] == 5:
                prefix_len = read.cigartuples[0][1]
            # UPDATE: refine the region with CIGAR tuples
            query_start = 0
            for op, length in read.cigartuples:  # M=0, I=1, D=2, N=3(ref_skip), S=4, H=5, Equal=7, Diff=8
                if op == 1:  # I=1 (Insertion to the reference, consumes read but not reference)
                    query_start += length
                elif op == 2:  # D=2 (Deletion from the reference, consumes reference but not read)
                    if length <= 20:
                        for i in range(start, start+length):
                            coverage[i] += 1
                    start += length
                elif op == 3:  # N=3 (Skipped region from the reference, like intron, consumes reference)
                    start += length
                elif op == 4 or op == 5:  # S=4 (Soft clip on the read, consumes read but not used in alignment)
                    # query_start += length
                    pass
                elif op in [0, 7, 8]:   # Match/mismatch, consumes both
                    if read.is_reverse:
                        dsa_map_region[read.query_name].append([prefix_len + read.query_length - query_start - length, prefix_len+read.query_length - query_start])
                    else:
                        dsa_map_region[read.query_name].append([prefix_len + query_start, prefix_len + query_start + length])
                    query_start += length
                    for i in range(start, start+length):
                        coverage[i] += 1
                    start += length
            # dsa_map_region[read.query_name].append([prefix_len+read.query_alignment_start, prefix_len+read.query_alignment_end])
        # track regions, find all the regions >= 1kb [no restriction on size]
        chrom_regions = []
        current_cov = None
        for pos, cov in enumerate(coverage):
            if current_cov is None:
                current_cov = cov
                start = pos
            elif current_cov != cov:
                end = pos
                # if end - start >= 1000:
                chrom_regions.append([chrom, start+1, end, current_cov])    # 1-based
                current_cov = cov
                start = pos
            else:
                continue
        if len(chrom_regions) == 0 or chrom_regions[-1][2] != chroms[chrom]:
            # if chroms[chrom] - start >= 1000:
            chrom_regions.append([chrom, start+1, chroms[chrom], current_cov])
        # print(f"Chrom {chrom} has {len(chrom_regions)} regions >= 1kb")
        return chrom_regions, dsa_map_region
    bed_file = output_file.replace('.chain.gz', '.bed')
    with open(bed_file, "w") as f:
        for chrom in chroms.keys():
            regions, dsa_mapped = cal_coverage(chrom)
            for region in regions:
                f.write("{}\t{}\t{}\t{}\n".format(region[0], region[1], region[2], region[3]))

def generate_chain(aln_bam, output_file=None, mapq_threshold=1, fa_file=None):
    '''
    Generate chain file from hap_to_ref alignment bam file.
    '''
    start_time = timeit.default_timer()
    if output_file is None:
        output_file = aln_bam.replace('.bam', '.chain.gz')
    elif not output_file.endswith('.gz'):
        output_file += '.gz'
    hap_to_ref_bam = AlignmentFile(aln_bam, 'rb')
    if fa_file is None:
        hap_contig = extract_contig_lengths(aln_bam.replace('_hg38.bam', '.fa').replace('_chm13.bam', '.fa'))
    else:
        hap_contig = extract_contig_lengths(fa_file)
    with gzip.open(output_file, 'wt') as f:
        f.write("#contig\tifReverse\tstart\tend\tref\tref_start\tref_end\tsize\n")
        for read in hap_to_ref_bam.fetch():
            if read.mapping_quality < mapq_threshold:
                continue
            contig_length = hap_contig[read.query_name]
            prefix_len = 0
            if read.query_length == contig_length:
                if read.is_reverse:
                    start = contig_length - read.query_alignment_end
                    end = contig_length - read.query_alignment_start
                else:
                    start = read.query_alignment_start
                    end = read.query_alignment_end
            else:
                if read.is_reverse:
                    if read.cigartuples[-1][0] == 5:
                        prefix_len = read.cigartuples[-1][1]
                else:
                    if read.cigartuples[0][0] == 5:
                        prefix_len = read.cigartuples[0][1]
                start = read.query_alignment_start + prefix_len
                end = read.query_alignment_end + prefix_len
            # contig strand start end ref start end size
            if read.is_reverse:
                ref_end = read.reference_end
                for op, length in read.cigartuples[::-1]:
                    if op == 1: # I=1 (Insertion to the reference, consumes read but not reference)
                        start += length
                    elif op == 2: # D=2 (Deletion from the reference, consumes reference but not read)
                        ref_end -= length
                    elif op == 3: # N=3 (Skipped region from the reference, like intron, consumes reference)
                        ref_end -= length
                    elif op == 4 or op == 5: # S=4 (Soft clip on the read, not affecting alignment) H=5 (Hard clip on the read, affecting alignment)
                        pass
                    elif op in [0, 7, 8]:   # Match/mismatch, consumes both
                        f.write(f"{read.query_name}\t{int(read.is_reverse)}\t{start+1}\t{start+length+1}\t{read.reference_name}\t{ref_end-length+1}\t{ref_end+1}\t{length}\n")
                        ref_end -= length
                        start += length
            else:
                ref_start = read.reference_start
                for op, length in read.cigartuples:
                    if op == 1:
                        start += length
                    elif op == 2:
                        ref_start += length
                    elif op == 3:
                        ref_start += length
                    elif op == 4 or op == 5:
                        pass
                    elif op in [0, 7, 8]:
                        f.write(f"{read.query_name}\t{int(read.is_reverse)}\t{start+1}\t{start+length+1}\t{read.reference_name}\t{ref_start+1}\t{ref_start+length+1}\t{length}\n")
                        ref_start += length
                        start += length
    print(f"Chain file has been generated: {output_file} in {timeit.default_timer()-start_time:.2f} seconds.")
    start_time = timeit.default_timer()
    get_coverage_bam(aln_bam, output_file, mapq_threshold)
    print("Bed file has been generated in {:.2f} seconds.".format(timeit.default_timer()-start_time))
    

def chain_main(params=None):
    parser = argparse.ArgumentParser(
        prog="chain", 
        description="Generate chain file from ref1_to_ref2 alignment bam file."
    )
    parser.add_argument(
        '-b', '--input_bam',
        type=str,
        required=True,
        help='Input ref1_to_ref2 alignment bam file.'
    )
    parser.add_argument(
        '-o', '--output_chain',
        type=str,
        default=None,
        help='Output chain file. Default: input.chain.gz'
    )
    parser.add_argument(
        '-q', '--mapq_threshold',
        type=int,
        default=1,
        help='Mapping quality threshold. Default: 1'
    )
    parser.add_argument(
        '-f', '--fasta',
        type=str,
        default=None,
        help='Fasta file of the ref1. Default: None'
    )

    # Allow passing a string or list of args programmatically
    if params is not None:
        if isinstance(params, str):
            params = params.split()
        args = parser.parse_args(params)
    else:
        args = parser.parse_args()

    generate_chain(
        args.input_bam,
        args.output_chain,
        args.mapq_threshold,
        args.fasta
    )

if __name__ == "__main__":
    chain_main()