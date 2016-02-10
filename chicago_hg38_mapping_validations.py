"""
Mapping chicago reads to hg38 with bwa and validating against pre-phased SNPs
"""

import sys
import os
import pysam
import vcf
import itertools
import argparse
import cPickle as pickle
import math
from pyfasta import Fasta
import pandas as pd
from collections import defaultdict, Counter
from pybedtools import BedTool
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from lib.seq_lib import ChromosomeInterval
from lib.general_lib import format_ratio, mkdir_p, DefaultOrderedDict
from sonLib.bioio import fastaRead
from collections import OrderedDict
# i need to sort the BAM by read tag
bam_handle = pysam.Samfile("CP499_chicago_500kb_notch2nl.bam.pairs.fq.hg38.read_sorted.bam")
read_dict = defaultdict(list)
for read in bam_handle:
    read_dict[read.qname].append(read)


def find_read_positions(chrom_positions, read):
    read_positions = []
    for read_pos, ref_pos in read.aligned_pairs:
        if ref_pos in chrom_positions:
            read_positions.append(read_pos)
    assert len(read_positions) > 0
    return read_positions


def build_sun_gt_map(vcf_rec):
    gt_map = {}
    for s in vcf_rec.samples:
        if s.gt_bases not in gt_map:
            gt_map[s.gt_bases] = []  # avoid defaultdict to prevent bugs
        gt_map[s.gt_bases].append(s.sample)
    return gt_map


def compare_sun_to_reference(intersections, read):
    for vcf_rec, vcf_interval in intersections:
        chrom_positions = range(vcf_rec.start, vcf_rec.end)
        read_positions = find_read_positions(chrom_positions, read)
        if len(chrom_positions) != len(read_positions) or None in read_positions:
            continue  # can't handle insertions right now
        read_seq = "".join([read.seq[i] for i in read_positions])
        #read_seq = "".join([seqs[read.qname][i] for i in read_positions])
        gt_map = build_sun_gt_map(vcf_rec)
        if read_seq not in gt_map:
            yield "NoMatches", vcf_rec, read
        elif len(gt_map[read_seq]) == 1:
            yield gt_map[read_seq][0], vcf_rec, read


def build_snp_gt_map(vcf_rec):
    gt_map = {}
    ps = vcf_rec.samples[0]['PS']
    bases = vcf_rec.samples[0].gt_bases.split("|")
    for i, b in enumerate(bases):
        gt_map[b] = "{}_{}".format(ps, i)
    return gt_map


def compare_snp_to_reference(intersections, read):
    for vcf_rec, vcf_interval in intersections:
        chrom_positions = range(vcf_rec.start, vcf_rec.end)
        read_positions = find_read_positions(chrom_positions, read)
        if len(chrom_positions) != len(read_positions) or None in read_positions:
            continue  # can't handle insertions right now
        read_seq = "".join([read.seq[i] for i in read_positions])
        #read_seq = "".join([seqs[read.qname][i] for i in read_positions])
        gt_map = build_snp_gt_map(vcf_rec)
        if read_seq not in gt_map:
            yield "NoHaplotype", "SnpNoMatches", vcf_rec, read
        else:
            yield gt_map[read_seq], vcf_rec, read


def find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
    sun_intersections = [[sun_recs[i], x] for i, x in enumerate(sun_intervals) if x.proper_subset(read_interval)]
    snp_intersections = [[snp_recs[i], x] for i, x in enumerate(snp_intervals) if x.proper_subset(read_interval)]
    sun_results = list(compare_sun_to_reference(sun_intersections, read))
    snp_results = list(compare_snp_to_reference(snp_intersections, read))
    return snp_results, sun_results


snp_h = vcf.Reader(open("/hive/users/ifiddes/longranger-1.2.0/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz"))
snp_recs = [x for x in snp_h.fetch('chr1', 119989614 - 550000, 149469354 + 550000) if not (x.is_indel is False and x.is_deletion is True)]
v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels.vcf.gz"))
sun_recs = [x for x in v_h if not (x.is_indel is False and x.is_deletion is True)]
sun_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in sun_recs]
snp_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in snp_recs]


analyzed_reads = defaultdict(lambda: defaultdict(list))
for qname, reads in read_dict.iteritems():
    if len(reads) != 2 or any([x.is_unmapped for x in reads]):
        continue
    for read in reads:
        snp_results, sun_results = find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle)
        if len(snp_results) > 0:
            analyzed_reads[qname]['snp'].append(snp_results)
        if len(sun_results) > 0:
            analyzed_reads[qname]['sun'].append(sun_results)


# TODO: try using this information to bridge together 10x-derived phase blocks.