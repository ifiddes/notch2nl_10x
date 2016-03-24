"""
Investigating the region of NUDT4P near Notch2NL-A that appears to have segregating
deletions.
"""

roi = ['chr1', 146231874, 146358436]

na12878_bam = '/hive/users/ifiddes/longranger-1.2.0/NA12878_combined/outs/phased_possorted_bam.bam'

import pysam
import numpy as np
from extract_linked_bams import bin_reads

bam_handle = pysam.Samfile(na12878_bam)
reads = bin_reads(roi[0], roi[1], roi[2], bam_handle, offset=10000)
# 4730 tags
# now, I want to grab any read with these tags, mapped or unmapped.


def eval_rec(rec, read_dict):
    tags = dict(rec.tags)
    if 'BX' in tags:
        t = tags['BX']
        if t in read_dict:
            read_dict[t].append(rec)



full_reads = reads.copy()
for rec in bam_handle.fetch(until_eof=True):
    eval_rec(rec, full_reads)

# while that is running, I am going to compare doing the same thing but with only notch regions

# adding in notch region reads
notch_region = ['chr1', 144950000, 150750000]
notch_reads = reads.copy()
for rec in bam_handle.fetch(notch_region[0], notch_region[1], notch_region[2]):
    eval_rec(rec, notch_reads)

# i am not sure how to get pysam to fetch ONLY unmapped, so using samtools
# samtools view -b phased_possorted_bam.bam '*' > unmapped.bam

reads_per_bin = [len(x) for x in notch_reads.itervalues()]
# average of 72 reads per bin

unmapped_bam = '/hive/users/ifiddes/longranger-1.2.0/NA12878_combined/outs/unmapped.bam'
notch_and_unmapped = notch_reads.copy()
unmapped_handle = pysam.Samfile(unmapped_bam)
for rec in unmapped_handle.fetch(until_eof=True):
    eval_rec(rec, notch_and_unmapped)


# how many did we add?
np.mean([len(x) for x in notch_and_unmapped.itervalues()])
# 293: a lot!

# I am not sure of the best way to extract the barcodes we are interested in
# first idea: I zoomed in on the browser until I reached the deletion window
deletion_window = ['chr1', 146281804, 146308501, '.']
sys.path.append("/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio")
from pycbio.bio.intervals import ChromosomeInterval
interval = ChromosomeInterval(*deletion_window)
overlaps = {}
for tag, reads in notch_and_unmapped.iteritems():
    num_overlap = 0
    tot_mapped = 0
    num_unmapped = 0
    for read in reads:
        if read.is_unmapped:
            num_unmapped += 1
            continue
        tot_mapped += 1
        r = zip(*read.aligned_pairs)[1]
        r = filter(None, r)
        l, r = min(r), max(r)
        c = ChromosomeInterval('chr1', l, r, '.')
        if c.proper_subset(interval):
            num_overlap += 1
    overlaps[tag] = [num_overlap, tot_mapped, num_unmapped]


# lets explore this dataset.
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
df = pd.DataFrame.from_dict(overlaps)
df = df.transpose()
df.columns = ['Number Overlaps', 'Total Mapped', 'Number Unmapped']
p = sns.pairplot(df)
p.savefig("pairplot.pdf", format='pdf')

# based on this pairplot, I think we can safely pull anything with exactly 0 unmapped
to_investigate = {}
for tag, reads in notch_and_unmapped.iteritems():
    num_overlap = overlaps[tag][0]
    if num_overlap == 0:
        to_investigate[tag] = reads

# that only dropped the number of tags from 4730 to 3714
# lets also filter for number mapped >= 10 and unmapped >= 50
# because we need to anchor things

to_investigate_v2 = {}
for tag, reads in notch_and_unmapped.iteritems():
    num_overlap, tot_mapped, num_unmapped = overlaps[tag]
    if num_overlap == 0 and tot_mapped >= 10 and num_unmapped >= 50:
        to_investigate_v2[tag] = reads


# that drops us to 1706 tags
# Ok... for each tag, extract all unmapped reads to a fasta, marking them.
with open("reads_to_map.fa", "w") as outf:
    for tag, reads in to_investigate_v2.iteritems():
        for read in reads:
            if read.is_unmapped:
                outf.write(">{}\n{}\n".format(tag, read.seq))


# BLAT these:
# blat ~/ifiddes_hive/chr1_fastas/hg38_chr1.fa reads_to_map.fa mapped.psl

# try chaining the results
# simpleChain mapped.psl -outPsl mapped.chained.psl
# num lines went from 26174 to 23270
# just convert to BED12 for each tag
from collections import defaultdict
recs = defaultdict(list)
from pycbio.bio.psl import PslRow
for l in open("mapped.psl"):
    p = PslRow(l.split())
    recs[p.q_name].append(p)


# 1360 tags have mappings - we need to evaluate the other 400!
filtered_recs = defaultdict(list)
for tag, psls in recs.iteritems():
    for p in psls:
        if p.block_count == 1 and p.coverage > 80:
            filtered_recs[tag].append(p)


np.mean([len(x) for x in filtered_recs.values()])

# this is dumb. lets just generate a coverage track and visualize.

# going to try manually investigating a few barcodes
for tag, reads in to_investigate.iteritems():
    if len(reads) > 10:
        for i, read in enumerate(reads):
            print ">{}_{}\n{}".format(tag, i, read.seq)
        assert False


# I took these 45 reads (4 of which are mapped) and BLAT'd against the human genome.
# I got hits for #0,1,2,3,7,8,10,11,12,14,15,19,20,23,24 and no others
# I stuck all into NCBI BLAST and looked around
# 19, 20, 21, 28, 37, 38 all had no hits
# the remainder all seem to hit the Carp genome...

# try one more bin for sanity's sake
for tag, reads in to_investigate.iteritems():
    if len(reads) > 10 and len(reads) != 45:
        for i, read in enumerate(reads):
            if read.is_unmapped:
                print ">{}_{}\n{}".format(tag, i, read.seq)
        assert False

# this bin has 184 reads!
# we can only BLAT 25 at a time though.



# this is stupid. I need a generalized script that converts mapped barcoded reads to BED12


from pycbio.bio.intervals import gap_merge_intervals

def get_read_interval(read):
    positions = zip(*read.aligned_pairs)[1]
    positions = filter(None, positions)
    l, r = min(positions), max(positions)
    return ChromosomeInterval(read.reference_id, l, r, '.')


def reads_to_merged_intervals(reads):
    intervals = [get_read_interval(x) for x in reads if not x.is_unmapped]
    intervals = sorted(intervals)
    intervals = gap_merge_intervals(intervals, 0)
    return intervals


def intervals_to_bed12(intervals, name, chrom=None):
    if chrom is None:
        chrom = intervals[0].chromosome
    chrom_start = intervals[0].start
    chrom_end = intervals[-1].stop
    strand = '.'
    score = 0
    thick_start = thick_end = 0
    rgb = '128,0,0'
    block_count = len(intervals)
    block_sizes = [len(x) for x in intervals]
    block_starts = [0]
    for i in intervals[1:]:
        block_starts.append(i.start - chrom_start)
    return map(str, [chrom, chrom_start, chrom_end, name, score, strand, thick_start,
                thick_end, rgb, block_count, ','.join(map(str, block_sizes)),
                ','.join(map(str, block_starts))])


# lets use this on the to_investigate_v2 set as well as the notch_and_unmapped set

with open("to_investigate.bed", "w") as outf:
    outf.write('track name="Possibly deletion haplotype"\n')
    for tag, reads in to_investigate_v2.iteritems():
        intervals = reads_to_merged_intervals(reads)
        rec = intervals_to_bed12(intervals, tag, chrom='chr1')
        outf.write('\t'.join(rec) + '\n')


with open("notch_and_unmapped.bed", "w") as outf:
    outf.write('track name="All tags in region"\n')
    for tag, reads in notch_and_unmapped.iteritems():
        intervals = reads_to_merged_intervals(reads)
        rec = intervals_to_bed12(intervals, tag, chrom='chr1')
        outf.write('\t'.join(rec) + '\n')
