"""
10x genomics haplotype analysis making use of PS tag

Version 5 plan:

Haynes told me that when a barcode is phased, only reads between stretches of MapQ are included. We want to relax this.
I also did not get the barcode sorting code to work. I don't think it matters. I will simply extract all reads mapping
to all notch loci, like before. For each barcode, if any of the reads are phased, phase all reads with MapQ>mapq_cutoff.
Discard all others.
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

mapq_cutoff = int(sys.argv[1])

genome = 'NA12878'
regions = [['chr1', 119989614, 120087517, 'Notch2'], 
           ['chr1', 146149503, 146248224, 'Notch2NL-A'], 
           ['chr1', 148600445, 148698970, 'Notch2NL-B'], 
           ['chr1', 149374498, 149469354, 'Notch2NL-C'], 
           ['chr1', 120707777, 120799019, 'Notch2NL-D']]


bam_path = "/hive/users/ifiddes/longranger-1.2.0/{}/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam".format(genome)
bam_handle = pysam.Samfile(bam_path)


def find_phased_barcodes(chrom, start, stop, bam_handle, offset=10000, mapq_cutoff=10):
    aln_iter = bam_handle.fetch(chrom, start - offset, stop + offset, multiple_iterators=True)
    phased_barcodes = {}
    recs = []
    for rec in aln_iter:
        if rec.mapq < mapq_cutoff:
            continue
        tags = dict(rec.tags)  # pysam is weird - why is this not a dict?
        if 'PS' in tags:
            t = (tags['PS'], tags['HP'])
            phased_barcodes[tags['BX']] = t
        recs.append(rec)
    return phased_barcodes, recs



phased_barcodes_holder = {}
read_holder = {}
for chrom, start, stop, name in regions:
    phased_barcodes_holder[name], recs = find_phased_barcodes(chrom, start, stop, bam_handle, mapq_cutoff=mapq_cutoff)
    read_holder[name] = recs


phased_read_holder = {}
rescued_counts = defaultdict(lambda: Counter())
tot_counts = defaultdict(lambda: Counter())
for name in read_holder:
    phased_read_holder[name] = defaultdict(list)
    for read in read_holder[name]:
        tags = dict(read.tags)
        if 'BX' in tags and tags['BX'] in phased_barcodes_holder[name]:
            if 'PS' not in tags:
                rescued_counts[name][tags['BX']] += 1
            ps = phased_barcodes_holder[name][tags['BX']]
            phased_read_holder[name][ps].append(read)
            tot_counts[name][tags['BX']] += 1


d = DefaultOrderedDict(list)
for para in sorted(tot_counts):
    # fix rescued_counts
    for tag in tot_counts[para]:
        if tag not in rescued_counts[para]:
            rescued_counts[para][tag] = 0
    tc = sorted(tot_counts[para].items(), key=lambda x: x[0])
    tc = zip(*tc)[1]
    c = sorted(rescued_counts[para].items(), key=lambda x: x[0])
    c = zip(*c)[1]
    d["Total Reads"].extend(tc)
    d["Rescued Reads"].extend(c)
    d["Paralog"].extend([para] * len(c))


from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.DataFrame.from_dict(d)
p = sns.pairplot(df, hue="Paralog", plot_kws=dict(s=8))
p.savefig("linked_bam_analysis/rescued_reads_{}.pdf".format(mapq_cutoff), format="pdf")


# now use code from previous versions to see how this changes the plots.


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
    evaled_vcfs = []
    for read in reads:
        # ignore un/poorly aligned reads
        if len(read.aligned_pairs) < min_aln_size:
            continue
        read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
        intersections = [[vcf_recs[i], x] for i, x in enumerate(vcf_intervals) if x.proper_subset(read_interval)]
        if len(intersections) == 0:
            continue
        for vcf_rec, vcf_interval in intersections:
            chrom_positions = range(vcf_rec.start, vcf_rec.end)
            read_positions = find_read_positions(chrom_positions, read)
            if len(chrom_positions) != len(read_positions) or None in read_positions:
                continue  # can't handle insertions right now
            read_seq = "".join([read.seq[i] for i in read_positions])
            gt_map = build_gt_map(vcf_rec)
            if read_seq not in gt_map:
                evaled_vcfs.append(["NoMatches", vcf_interval.start])
            elif len(gt_map[read_seq]) == 1:
                evaled_vcfs.append([gt_map[read_seq][0], vcf_interval.start])
            else:
                evaled_vcfs.append(["MultipleParalogs", vcf_interval.start])
    return evaled_vcfs



v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels_AllPerspectives.vcf.gz"))
read_map_holder = {}
region_map = {x[-1]: x for x in regions}
for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D", "Notch2"]:
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


# build data structure
data_holder = defaultdict(dict)
for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D", "Notch2"]:
    for phase_set in read_map_holder[para]:
        data_holder[para][phase_set] = defaultdict(Counter)
        for snp_para, pos in read_map_holder[para][phase_set]:
            pos = pos_map_inverted[para][pos]
            data_holder[para][phase_set][pos][snp_para] += 1
        # filter out entries that don't matter
        to_keep = {}
        for pos, c in data_holder[para][phase_set].iteritems():
            if len(set(names) & c.viewkeys()) != 0:
                to_keep[pos] = sorted(c.items(), reverse=True)
        data_holder[para][phase_set] = to_keep


sorted_data_holder = OrderedDict()
for para in data_holder:
    sorted_data_holder[para] = OrderedDict()
    for phase_set, d in data_holder[para].iteritems():
        d = sorted(d.iteritems(), key=lambda x: x[0])
        sorted_data_holder[para][phase_set] = d


import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
sns.set_style("whitegrid")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter

colors = {"NoMatches": sns.color_palette()[-1], "MultipleParalogs": "grey", None: "grey"}
for i, n in enumerate(names):
    colors[n] = sns.color_palette()[i]


def find_largest(d):
    vals = []
    for x in d:
        s = 0
        for y in x[1]:
            s += y[1]
        vals.append(s)
    return log_transform(max(vals))


def log_transform(x):
    return math.log(x + math.e - 1)


with PdfPages('linked_bam_analysis/{}_phased_haplotypes_rescued_reads_{}.pdf'.format(genome, mapq_cutoff)) as pdf:
    for para in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D", "Notch2"]:
        fig, plots = plt.subplots(len(sorted_data_holder[para]), sharey=True, sharex=True)
        for i, tag in enumerate(sorted(sorted_data_holder[para].keys())):
            d = sorted_data_holder[para][tag]
            p = plots[i]
            p.set_ylabel("Haplotype {}".format(i), fontsize=9)
            p.set_yticklabels([])
            plt.ylim((0, find_largest(d)))
            plt.xticks((0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))
            plt.xlim((0, 100143))
            for pos, counts in d:
                offset = 0
                for name, c in counts:
                    v = log_transform(c)
                    _ = p.vlines(x=pos, ymin=offset, ymax=v, color=colors[name])
                    offset += v
        plt.suptitle("Visualization of {} phased reads in genome {} with MapQ cutoff {}".format(para, genome, mapq_cutoff))
        plt.figtext(0.01, 0.65, "log(# reads + e - 1)", rotation="vertical", fontsize=9)
        plt.xlabel("Notch2NL Consensus Position")
        allnames = list(names) + ["NoMatches", "MultipleParalogs"]
        my_patches = [patches.Patch(color=colors[n], label=n) for n in allnames]
        fig.legend(handles=my_patches, labels=allnames, title="SUN Matches", fontsize=9, loc="center right")
        plt.tight_layout(rect=[0.01, 0.07, 0.86, 0.94])
        pdf.savefig()
        plt.close('all')
