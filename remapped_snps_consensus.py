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
# We need the start positions for each in the actual genome, used a browser BLAT
start_pos = {"Notch2": 120087516, "Notch2NL-A": 146248223, "Notch2NL-B": 148698969, "Notch2NL-C": 149374496, "Notch2NL-D": 120707775}
# which of these are backwards?
backwards = {"Notch2NL-C", "Notch2NL-D"}

names = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D'] # same as in VCF

genome = 'NA12878'
regions = [['chr1', 119989614, 120087517, 'Notch2'], 
           ['chr1', 146149503, 146248224, 'Notch2NL-A'], 
           ['chr1', 148600445, 148698970, 'Notch2NL-B'], 
           ['chr1', 149374498, 149469354, 'Notch2NL-C'], 
           ['chr1', 120707777, 120799019, 'Notch2NL-D']]


# build a map of alignment positions to sequence positions
r = {name: seq for name, seq in fastaRead("/hive/users/ifiddes/notch2nl_suns/notch2_aligned.fasta")}
r_sort = sorted(r.iteritems(),key=lambda x: x[0])
names, seqs = zip(*r_sort)
tgt_is = {n: 0 for n in names}

# map sequence coordinates to alignment coordinates
pos_map = defaultdict(dict)
for ref_i, cs in enumerate(zip(*seqs)):
    for name, tgt_i in tgt_is.iteritems():
        pos_map[name][ref_i] = tgt_i
    for name, c in zip(*[names, cs]):
        if c != "-":
            tgt_is[name] += 1


# map genome coordinates to alignment coordinates
pos_map_inverted = defaultdict(dict)
for para, vals in pos_map.iteritems():
    if para not in backwards:
        vals = {start_pos[para] - y: x for x, y in vals.iteritems()}
    else:
        vals = {start_pos[para] + y: x for x, y in vals.iteritems()}
    pos_map_inverted[para] = vals


# map alignment coordinates to genome coordinates
pos_map_genome = defaultdict(dict)
for para, vals in pos_map_inverted.iteritems():
    pos_map_genome[para] = {y: x for x, y in vals.iteritems()}


snp_vcf = vcf.Reader(open("/hive/users/ifiddes/longranger-1.2.0/platgene_NA12878-hg38-2_0_1.vcf.gz"))
#snp_vcf = vcf.Reader(open("/hive/users/ifiddes/longranger-1.2.0/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz"))


consensus_fasta = Fasta("/hive/users/ifiddes/notch2nl_suns/notch2_aligned_consensus.fasta")['Notch2NL_consensus']


import string
def rev_comp(s, rc=string.maketrans("ATGC", "TACG")):
    return string.translate(s, rc)[::-1]


# we need to build a set of any mismatches in the consensus - I.E., we want only to evaluate columns that all match.
bad_columns = []
for i, x in enumerate(zip(*seqs)):
    if len(Counter(x)) != 1:
        bad_columns.append(i)


# merge intervals
from operator import itemgetter
from itertools import groupby
f = []
for k, g in groupby(enumerate(bad_columns), lambda (i,x):i-x):
    f.append(map(itemgetter(1), g))

# this filters out deletions, maybe it won't matter
bad_intervals = [ChromosomeInterval('Notch2NL_consensus', x[0], x[-1] + 1, ".") for x in f]

import copy
remapped_snps = []
for chrom, start, stop, para in regions:
    vcf_recs = snp_vcf.fetch(chrom, start, stop)
    for rec in vcf_recs:
        if rec.heterozygosity == 0.5 and rec.samples[0].phased is True:
            orig_rec = copy.deepcopy(rec)
            rec.CHROM = 'Notch2NL_consensus'
            rec.add_info("PARA", para)
            if para not in backwards:
                rec.REF = rev_comp(rec.REF)
                rec.ALT[0].sequence = rev_comp(rec.ALT[0].sequence)
                rec.start = pos_map_inverted[para][orig_rec.end]
                rec.end = pos_map_inverted[para][orig_rec.start]
                rec.POS = rec.end
            else:
                rec.start = pos_map_inverted[para][orig_rec.start] - 1
                rec.end = pos_map_inverted[para][orig_rec.end] - 1
                rec.POS = rec.end
            rec_interval = ChromosomeInterval(rec.CHROM, rec.start, rec.end, None)
            intersections = [True for i, x in enumerate(bad_intervals) if x.overlap(rec_interval)]
            if len(intersections) > 0:
                continue
            assert consensus_fasta[rec.start: rec.end] == rec.REF
            remapped_snps.append(rec)


with open("remapped_heterozygous_snps_IlluminaGoldStandard_NA12878.vcf", "w") as outf:
    w = vcf.Writer(outf, snp_vcf)
    for rec in remapped_snps:
        w.write_record(rec)

# we have 159 heterozygous SNPs for NA12878 that are not in fact SUNs and are remappable
# Counter({'Notch2': 51, 'Notch2NL-D': 50, 'Notch2NL-A': 40, 'Notch2NL-B': 12, 'Notch2NL-C': 6})
# not a lot of B/C specific SNPs...
# I am going to use the exact same code for the NA12878 gold standard VCF I found from Illumina
# Counter({'Notch2': 47, 'Notch2NL-D': 12, 'Notch2NL-A': 9, 'Notch2NL-C': 4})
# Gold standard has far fewer SNPs, but the ratios seem more realistic accurate.
# how many are unique to each set?

# do some BED intersections

from pybedtools import BedTool

ill = BedTool("remapped_heterozygous_snps_IlluminaGoldStandard_NA12878.vcf")
tenx = BedTool("remapped_heterozygous_snps_NA12878.vcf")
len(list(ill.intersect(tenx)))
# 60 overlaps overall - breakdown by paralog
Counter([x.fields[7].split(";")[-1].split("=")[1] for x in ill.intersect(tenx)])
# Counter({u'Notch2': 39, u'Notch2NL-D': 9, u'Notch2NL-A': 8, u'Notch2NL-C': 4})
# generally it seems that both sets agree, except for Notch2NL-B, where Illumina has nothing.
# Therefore, I am going to stick with what 10x has, until the NA12878 run with the input VCF is done.
