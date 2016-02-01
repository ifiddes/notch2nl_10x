"""
10x genomics data haplotype phaser for notch locus
"""
import sys
import os
import pysam
import vcf
import itertools
import argparse
import cPickle as pickle
from collections import defaultdict
from pybedtools import BedTool
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from lib.seq_lib import ChromosomeInterval
from lib.general_lib import format_ratio


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bamfile", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pickleFile", required=True)
    return parser.parse_args()


regions = [['chr1', 119989614, 120087517, 'Notch2'], 
           ['chr1', 146149503, 146248224, 'Notch2NL-A'], 
           ['chr1', 148600445, 148698970, 'Notch2NL-B'], 
           ['chr1', 149374498, 149469354, 'Notch2NL-C'], 
           ['chr1', 120707777, 120799019, 'Notch2NL-D']]


def bin_reads(chrom, start, stop, bam_handle, offset=10000):
    """
    bins reads based on their BX tag
    """
    aln_iter = bam_handle.fetch(region="{}:{}-{}".format(chrom, start - offset, stop + offset), multiple_iterators=True)
    binned_reads = defaultdict(list)
    for rec in aln_iter:
        tags = dict(rec.tags)  # pysam is weird - why is this not a dict?
        if 'BX' in tags:
            t = tags['BX']
            binned_reads[t].append(rec)
    return binned_reads


def build_vcf_intervals(reads, vcf_recs, bam_handle):
    """
    Find if any of these reads match a known SUN/indel by simple bedtools intersections
    """
    vcf_bed_recs = [ChromosomeInterval(x.CHROM, x.start, x.end, None) for x in vcf_recs]
    vcf_bed = BedTool(vcf_bed_recs)
    reads_bed_recs = [(bam_handle.getrname(x.tid), x.positions[0], x.positions[-1]) for x in reads if len(x.positions) > 2]
    reads_bed = BedTool(reads_bed_recs)
    return list(vcf_bed.intersect(reads_bed))


def find_read_positions(chrom_positions, read):
    read_positions = []
    for read_pos, ref_pos in read.aligned_pairs:
        if ref_pos in chrom_positions:
            read_positions.append(read_pos)
    assert len(read_positions) > 0
    return read_positions


def determine_if_bin_is_right_paralog(reads, vcf_recs, bam_handle, min_aln_size=37):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    vcf_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in vcf_recs]
    num_support = 0
    total_intersections = 0
    supporting_reads = []
    for read in reads:
        # ignore un/poorly aligned reads
        if len(read.aligned_pairs) < min_aln_size:
            continue
        read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
        intersections = [[vcf_recs[i], x] for i, x in enumerate(vcf_intervals) if x.proper_subset(read_interval)]
        for vcf_rec, vcf_interval in intersections:
            chrom_positions = range(vcf_rec.start, vcf_rec.end)
            read_positions = find_read_positions(chrom_positions, read)
            if len(chrom_positions) != len(read_positions) or None in read_positions:
                total_intersections += 1
                continue
            read_seq = "".join([read.seq[i] for i in read_positions])
            if read_seq == vcf_rec.REF:
                num_support += 1
                supporting_reads.append(read)
            total_intersections += 1
    return num_support, total_intersections, supporting_reads


def main():
    args = parse_args()
    bam_handle = pysam.Samfile(args.bamfile)
    v_h = vcf.Reader(open(args.vcf))
    counts = {}
    interesting_read_holder = {}
    for chrom, start, stop, name in regions:
        binned_reads = bin_reads(chrom, start, stop, bam_handle)
        vcf_recs = list(v_h.fetch(chrom, start, stop))
        # todo: allow insertions
        vcf_recs = [x for x in vcf_recs if not (x.is_indel is False and x.is_deletion is True)]
        c = {}
        reads_holder = {}
        for tag, reads in binned_reads.iteritems():
            num_support, total_intersections, supporting_reads = determine_if_bin_is_right_paralog(reads, vcf_recs, bam_handle)
            c[tag] = [num_support, total_intersections, len(reads)]
            if len(supporting_reads) > 0:
                nonsupport_reads = [x for x in reads if x not in supporting_reads]
                reads_holder[tag] = [supporting_reads, nonsupport_reads]
        counts[name] = c
        interesting_read_holder[name] = reads_holder
    pickle.dump(interesting_read_holder, open(args.pickleFile, "w"))


if __name__ == "__main__":
    main()