"""
10x genomics haplotype analysis making use of PS tag

Extract reads that are NOT phased and remap per-barcode to the collapsed consensus
"""
import sys
import os
import pysam
import vcf
import itertools
import argparse
import cPickle as pickle
import math
import pandas as pd
from collections import defaultdict, Counter
from pybedtools import BedTool
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from lib.seq_lib import ChromosomeInterval
from lib.general_lib import format_ratio, mkdir_p, DefaultOrderedDict

genome = 'NA12878'
regions = [['chr1', 119989614, 120087517, 'Notch2'], 
           ['chr1', 146149503, 146248224, 'Notch2NL-A'], 
           ['chr1', 148600445, 148698970, 'Notch2NL-B'], 
           ['chr1', 149374498, 149469354, 'Notch2NL-C'], 
           ['chr1', 120707777, 120799019, 'Notch2NL-D']]


bam_path = "/hive/users/ifiddes/longranger-1.2.0/{}/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam".format(genome)
bam_handle = pysam.Samfile(bam_path)


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



read_holder = defaultdict(set)
for chrom, start, stop, name in regions:
    binned_reads = bin_reads(chrom, start, stop, bam_handle)
    for tag, reads in binned_reads.iteritems():
        read_holder[tag].update(reads)


# this is ugly, but whatever. I will just extract all reads in this to fastq, remap, then build a mapping of read name
# to barcode to rescue the BX tags.

with pysam.Samfile("remap_reads/notch2nl_reads.bam", "wb", template=bam_handle) as outf:
    for r in read_holder.itervalues():
        for x in r:
            outf.write(x)



tag_holder = defaultdict(set)
for tag, reads in read_holder.iteritems():
    tag_holder[tag] = [x.qname for x in reads]


# lets map these
#samtools sort -O bam -T tmp notch2nl_reads.bam > notch2nl_reads.sorted.bam
# samtools index notch2nl_reads.sorted.bam
#samtools bamshuf -Ou notch2nl_reads.sorted.bam tmp | samtools bam2fq -s singles.fq - > pairs.fq
# lets try blat...
# fastqToFa singles.fq singles.fa
# fastqToFa pairs.fq pairs.fa
# cat pairs.fa singles.fa > reads.fa
# blat /hive/users/ifiddes/notch2nl_suns/test_index/notch2_aligned_consensus.fasta reads.fa blat.psl
# psl2sam.pl blat.psl > blat.sam
# samtools view -b blat.sam | samtools sort -T tmp -O bam - > blat.sorted.bam
# bwa mem /hive/users/ifiddes/notch2nl_suns/test_index/notch2_aligned_consensus.fasta pairs.fq -p > pairs_bwa.sam
# bwa mem /hive/users/ifiddes/notch2nl_suns/test_index/notch2_aligned_consensus.fasta singles.fq > singles_bwa.sam
# head -n 1 pairs_bwa.sam > blat.header.sam
# cat blat.sam >> blat.header.sam
# samtools view -b blat.header.sam | samtools sort -T tmp -O bam - > blat.sorted.bam
# samtools view -b pairs_bwa.sam | samtools sort -T tmp -O bam - > pairs_bwa.sorted.bam
# samtools view -b singles_bwa.sam | samtools sort -T tmp -O bam - > singles_bwa.sorted.bam
# samtools merge bwa.bam singles_bwa.sorted.bam  pairs_bwa.sorted.bam
# samtools index bwa.bam


inverted_tag_holder = {}
for tag, r in tag_holder.iteritems():
    for x in r:
        inverted_tag_holder[x] = tag


def fix_tags(path):
    h = pysam.Samfile(path)
    new_recs = defaultdict(list)
    for rec in h:
        if len(rec.aligned_pairs) < 50:
            continue
        if rec.qname[-2:] == "/2" or rec.qname[-2:] == "/1":
            qname = rec.qname[:-2]
        else:
            qname = rec.qname
        tag = inverted_tag_holder[qname]
        rec.tags.append(('BX', tag))
        new_recs[tag].append(rec)
    return new_recs


blat_recs = fix_tags("remap_reads/blat.sorted.bam")
bwa_recs = fix_tags("remap_reads/bwa.bam")


#>>> sum([len(x) for x in blat_recs.itervalues()])
#220907
#>>> sum([len(x) for x in bwa_recs.itervalues()])
#198029

# blat recs lack sequence
from pyfasta import Fasta
seqs = Fasta("remap_reads/reads.fa")


# now that we have fixed the reads, we need to somehow determine what paralog they fit to.
# I feel like this should be some fancy statistical model, but I am not sure how.
# I also think I should try and intersect this with the SNP calls produced by longranger
# I can extract only phased heterozygous calls and use my coordinate remapping to map to the consensus.
# Algorithm:

# for each aligned barcode:
#    if the aligned barcode matches two or more (phased) heterozygous SNPs:
#      count it as part of that phase set

# problems with this approach:
# I may not have enough heterozygous SNPs. The heterozygous SNP calls may be wrong. I don't think so, however -
# they are very conservative, which is why we have so few reads.

# lifted the heterozygous SNPs from the 10x results to my consensus

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
            yield gt_map[read_seq], vcf_rec.INFO['PARA'][0], vcf_rec, read


def overall_metrics(results):
    results["sun_matches"] = Counter([x for x in zip(*results["sun"])[0]]) if len([x for x in results["sun"] if x[0] != "NoMatches"]) > 0 else None
    results["snp_hit"] = True if len([x for x in results["snp"] if x[0] != "SnpNoMatches"]) > 0 else False


def find_read_snp_sun_intersections(reads, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    results = {"sun": [], "snp": []}
    for read in reads:
        read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
        sun_intersections = [[sun_recs[i], x] for i, x in enumerate(sun_intervals) if x.proper_subset(read_interval)]
        snp_intersections = [[snp_recs[i], x] for i, x in enumerate(snp_intervals) if x.proper_subset(read_interval)]
        if len(sun_intersections) == 0 and len(snp_intersections) == 0:
            continue
        sun_results = list(compare_sun_to_reference(sun_intersections, read))
        if len(sun_results) > 0:
            results["sun"].extend(sun_results)
        snp_results = list(compare_snp_to_reference(snp_intersections, read))
        if len(snp_results) > 0:
            results["snp"].extend(snp_results)
    overall_metrics(results)
    return results


snp_h = vcf.Reader(open("remapped_heterozygous_snps_consensus.vcf"))
snp_recs = [x for x in snp_h if not (x.is_indel is False and x.is_deletion is True)]
v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels_ConsensusRef.vcf.gz"))
sun_recs = [x for x in v_h if not (x.is_indel is False and x.is_deletion is True)]
bam_handle = pysam.Samfile("remap_reads/bwa.bam")  # using bwa because psl2sam.pl didn't put sequences in
sun_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in sun_recs]
snp_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in snp_recs]


analyzed_tags = {}
for tag, reads in bwa_recs.iteritems():
    analyzed_tags[tag] = find_read_snp_sun_intersections(reads, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle)


# now lets get some metrics - can we phase these?
len([x for x in analyzed_tags.itervalues() if x['snp_hit'] is not False and x['sun_matches'] is not None])
# 2073 potentially informative barcodes


# remap to haplotypes, look for entries that are incongruous
remapped_tags = defaultdict(list)
bad_tags = {}
for tag, results in analyzed_tags.iteritems():
    if results["snp_hit"] is True and results["sun_matches"] is not None:
        haplotypes, snp_paras, _, __ = zip(*results['snp'])
        r = zip(*results['snp'])
        if len(snp_paras) > 1:
            bad_tags[tag] = results
        else:
            remapped_tags[haplotypes[0]] = [tag, results]


# tags seen in entire 10x data:
tenx_tags = set(['40120000161_1', '40120000161_0', '30117450697_1', '30117450697_0', '40149418773_0', '40149418773_1', '40146126759_1', '40146126759_0', '40149366446_1', '40149366446_0'])
my_tags = set(remapped_tags.keys())
# associated with paralogs:
tenx_para_tags = {'40120000161_1': 'Notch2', '40120000161_0': 'Notch2', '30117450697_1': 'Notch2', '30117450697_0': 'Notch2', '40149418773_0': 'Notch2NL-C', '40149418773_1': 'Notch2NL-C', '40146126759_1': 'Notch2NL-B', '40146126759_0': 'Notch2NL-B', '40149366446_1': 'Notch2NL-C', '40149366446_0': 'Notch2NL-C'}

for tag, para in tenx_para_tags.iteritems():
    if tag not in my_tags:
        print tag, para


# 40149366446_1 Notch2NL-C
# 40149366446_0 Notch2NL-C
# 30117450697_0 Notch2

# lets try visualizing some bad ones
for tag, results in bad_tags.iteritems():
    haplotypes, snp_paras, _, __ = zip(*results['snp'])
    assert not('Notch2NL-A' in snp_paras and len(snp_paras) > 1)


with pysam.Samfile("test_reads.bam", "wb", template=bam_handle) as outf:
    for q in ['snp', 'sun']:
        reads = zip(*results[q])[-1]
        for read in reads:
            outf.write(read)