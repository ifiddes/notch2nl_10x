"""
10x genomics haplotype analysis making use of PS tag

Version 2 uses a improved VCF that has all SUNs marked on all paralogs.
However, this does not seem to work that well.

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
from lib.general_lib import format_ratio, mkdir_p


genome = 'NA12878'


regions = [['chr1', 119989614, 120087517, 'Notch2'], 
           ['chr1', 146149503, 146248224, 'Notch2NL-A'], 
           ['chr1', 148600445, 148698970, 'Notch2NL-B'], 
           ['chr1', 149374498, 149469354, 'Notch2NL-C'], 
           ['chr1', 120707777, 120799019, 'Notch2NL-D']]


def bin_phased_reads(chrom, start, stop, bam_handle, offset=10000):
    """
    bins reads based on their BX tag
    """
    aln_iter = bam_handle.fetch(region="{}:{}-{}".format(chrom, start - offset, stop + offset), multiple_iterators=True)
    phased_reads = defaultdict(list)
    for rec in aln_iter:
        tags = dict(rec.tags)  # pysam is weird - why is this not a dict?
        if 'PS' in tags:
            t = (tags['PS'], tags['HP'])
            phased_reads[t].append(rec)
    return phased_reads


bam_path = "/hive/users/ifiddes/longranger-1.2.0/{}/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam".format(genome)
bam_handle = pysam.Samfile(bam_path)
v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels_AllPerspectives.vcf.gz"))


phased_read_holder = {}
bam_handle = pysam.Samfile(bamfile)
for chrom, start, stop, name in regions:
    phased_read_holder[name] = bin_phased_reads(chrom, start, stop, bam_handle)


# validate these phasing results according to SUN positions


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


def build_gt_map(vcf_rec):
    gt_map = {}
    for s in vcf_rec.samples:
        if s.gt_bases not in gt_map:
            gt_map[s.gt_bases] = []  # avoid defaultdict to prevent bugs
        gt_map[s.gt_bases].append(s.sample)
    return gt_map


def find_read_sun_intersections(reads, vcf_recs, bam_handle, min_aln_size=37):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    vcf_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in vcf_recs]
    evaled_read_list = []
    for read in reads:
        # ignore un/poorly aligned reads
        if len(read.aligned_pairs) < min_aln_size:
            continue
        read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
        intersections = [[vcf_recs[i], x] for i, x in enumerate(vcf_intervals) if x.proper_subset(read_interval)]
        if len(intersections) == 0:
            evaled_read_list.append([None, read])
            continue
        for vcf_rec, vcf_interval in intersections:
            chrom_positions = range(vcf_rec.start, vcf_rec.end)
            read_positions = find_read_positions(chrom_positions, read)
            if len(chrom_positions) != len(read_positions) or None in read_positions:
                evaled_read_list.append([None, read])
                continue  # can't handle insertions right now
            read_seq = "".join([read.seq[i] for i in read_positions])
            gt_map = build_gt_map(vcf_rec)
            if read_seq not in gt_map:
                evaled_read_list.append(["NoMatches", read])
            elif len(gt_map[read_seq]) == 1:
                evaled_read_list.append([gt_map[read_seq][0], read])
            else:
                evaled_read_list.append(["MultipleParalogs", read])
    return evaled_read_list


read_map_holder = {}
region_map = {x[-1]: x for x in regions}
for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C"]:
    chrom, start, stop, para = region_map[para]
    read_map_holder[para] = {}
    vcf_recs = list(v_h.fetch(chrom, start, stop))
    # todo: allow insertions
    vcf_recs = [x for x in vcf_recs if not (x.is_indel is False and x.is_deletion is True)]
    for tag, reads in phased_read_holder[para].iteritems():
        read_map_holder[para][tag] = find_read_sun_intersections(reads, vcf_recs, bam_handle)


from sonLib.bioio import fastaRead
from collections import OrderedDict
# We need the start positions for each in the actual genome, used a browser BLAT
start_pos = {"Notch2": 120087516, "Notch2NL-A": 146248223, "Notch2NL-B": 148698969, "Notch2NL-C": 149374496, "Notch2NL-D": 120707775}
# which of these are backwards?
backwards = {"Notch2NL-C", "Notch2NL-D"}

names = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D'] # same as in VCF


# build a map of alignment positions to sequence positions
r = {name: seq for name, seq in fastaRead("/hive/users/ifiddes/notch2nl_suns/notch2_aligned.fasta")}
r_sort = sorted(r.iteritems(),key=lambda x: x[0])
names, seqs = zip(*r_sort)
tgt_is = {n: 0 for n in names}


pos_map = defaultdict(dict)
for ref_i, cs in enumerate(zip(*seqs)):
    for name, tgt_i in tgt_is.iteritems():
        pos_map[name][ref_i] = tgt_i
    for name, c in zip(*[names, cs]):
        if c != "-":
            tgt_is[name] += 1


# invert pos_map
pos_map_inverted = defaultdict(dict)
for para, vals in pos_map.iteritems():
    if para not in backwards:
        vals = {start_pos[para] - y: x for x, y in vals.iteritems()}
    else:
        vals = {start_pos[para] + y: x for x, y in vals.iteritems()}
    pos_map_inverted[para] = vals



# map interesting reads
read_holder = defaultdict(dict)
for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C"]:
    for phase_set in read_map_holder[para]:
        read_holder[para][phase_set] = []
        for snp_para, read in read_map_holder[para][phase_set]:
            read_pos, ref_pos = zip(*read.aligned_pairs)
            ref_pos = [x for x in ref_pos if x is not None]
            try:
                start = pos_map_inverted[para][min(ref_pos)]
                stop = pos_map_inverted[para][max(ref_pos)]
                left = min([start, stop])
                right = max([start, stop])
            except:
                continue
            read_pos = [snp_para, left, right]
            read_holder[para][phase_set].append(read_pos)


sorted_read_holder = OrderedDict()
for para in read_holder:
    sorted_read_holder[para] = OrderedDict()
    for phase_set, reads in read_holder[para].iteritems():
        reads = sorted(reads, key=lambda x: x[1])
        sorted_read_holder[para][phase_set] = reads


import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
sns.set_style("white")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter

colors = {"NoMatches": sns.color_palette()[-1], "MultipleParalogs": "grey", None: "grey"}
for i, n in enumerate(names):
    colors[n] = sns.color_palette()[i]


with PdfPages('linked_bam_analysis/{}_phased_haplotypes_v2.pdf'.format(genome)) as pdf:
    for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C"]:
        fig, plots = plt.subplots(len(sorted_read_holder[para]), sharey=True, sharex=True)
        for i, tag in enumerate(sorted(sorted_read_holder[para].keys())):
            reads = sorted_read_holder[para][tag]
            p = plots[i]
            p.set_ylabel("Haplotype {}".format(i))
            plt.yticks([])
            plt.ylim((0, 1))
            plt.xticks((0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))
            plt.xlim((0, 100143))
            p.hlines(color="black", y=0.5, xmin=0, xmax=100143, linewidth=0.5)
            for read_para, left, right in reads:
                if read_para is None:
                    continue
                elif read_para == "MultipleParalogs":
                    box = {"xy": (left, 0.4), "width": right - left, "height": 0.2, "color": colors[read_para], "alpha": 0.5}
                else:
                    box = {"xy": (left, 0.0), "width": right - left, "height": 1.0, "color": colors[read_para], "alpha": 0.8}
                _ = p.add_patch(patches.Rectangle(**box))
        plt.suptitle("Visualization of {} phased reads in genome {}".format(para, genome))
        plt.xlabel("Notch2NL Consensus Position")
        allnames = list(names) + ["NoMatches", "MultipleParalogs"]
        my_patches = [patches.Patch(color=colors[n], label=n) for n in allnames]
        fig.legend(handles=my_patches, labels=allnames, title="SUN Matches", fontsize=9, loc="center right")
        plt.tight_layout(rect=[0, 0.07, 0.86, 0.94])
        pdf.savefig()
        plt.close('all')
