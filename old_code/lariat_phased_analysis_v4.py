"""
10x genomics haplotype analysis making use of PS tag

Version 4 plan:

1. Remap all reads of each haplotype to my consensus reference. 
2. Build a per-haplotype VCF representing non-SUN variants found on this haplotype.
3. Start pulling out unphased GEMs that map to the Notch locus and see if I can phase them.
"""
import sys
import os
import pysam
import vcf
import itertools
import argparse
import cPickle as pickle
import math
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
phased_read_holder = {}
bam_handle = pysam.Samfile(bamfile)
for chrom, start, stop, name in regions:
    phased_read_holder[name] = bin_phased_reads(chrom, start, stop, bam_handle)


out_dir = "/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/linked_bam_analysis/split_bams/{}".format(genome)
mkdir_p(out_dir)
for para in phased_read_holder:
    for tag in phased_read_holder[para]:
        with pysam.Samfile(os.path.join(out_dir, "{}.{}.bam".format(para, "_".join(map(str, tag)))), "wb", template=bam_handle) as outf:
            for read in phased_read_holder[para][tag]:
                outf.write(read)

fastq_dir = "/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/linked_bam_analysis/split_fastqs/{}".format(genome)
mkdir_p(fastq_dir)
fastqs = {}
for bam in [x for x in os.listdir(out_dir) if x.endswith("bam")]:
    outbase = os.path.join(fastq_dir, bam.replace(".bam", ""))
    bampath = os.path.join(out_dir, bam)
    outpaired = outbase + ".paired.fq"
    outsingle = outbase + ".single.fq"
    fastqs[tuple(bam.replace(".bam", "").split("."))] = [outpaired, outsingle]
    subprocess.call("samtools bamshuf -Ou {} tmp | samtools bam2fq -s {} - > {}".format(bampath, outsingle, outpaired), shell=True)


index = "/hive/users/ifiddes/notch2nl_suns/notch2_aligned_consensus.fasta"

aligned_dir = "/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/linked_bam_analysis/split_alignments/{}".format(genome)
mkdir_p(aligned_dir)

paired_cmd = "bwa mem -t 4 -p {} {} | samtools view -b - | samtools sort -O bam -T tmp - > {}"
single_cmd = "bwa mem -t 4 {} {} | samtools view -b - | samtools sort -O bam -T tmp - > {}"


merged_bams = {}
for haplotype, (pairedfastq, singlefastq) in fastqs.iteritems():
    outbase = os.path.join(aligned_dir, ".".join(haplotype))
    outpair = outbase + ".paired.sorted.bam"
    outsingle = outbase + ".single.sorted.bam"
    outmerged = outbase + ".merged.bam"
    subprocess.call(paired_cmd.format(index, pairedfastq, outpair), shell=True)
    subprocess.call(single_cmd.format(index, singlefastq, outsingle), shell=True)
    subprocess.call("samtools merge {} {} {}".format(outmerged, outsingle, outpair), shell=True)
    subprocess.call("samtools index {}".format(outmerged), shell=True)
    merged_bams[haplotype] = outmerged


# produce basic variant calls
freebayes_cmd = "freebayes -f {} -p 1 -0 {} > {} &"
for haplotype, bam in merged_bams.iteritems():
    out_vcf = os.path.join(aligned_dir, ".".join(haplotype) + ".vcf")
    cmd = freebayes_cmd.format(index, bam, out_vcf)
    subprocess.call(cmd, shell=True)

# NEXT STEPS:
# reverse intersect these VCFs by finding variants unique to only a haplotype (non other, non-SUN) that pass mapQ.
# validate these calls against what 10x did - maybe you should convert 10x coordinates and use as input VCF?
# next, begin remapping procedure - extract GEMs that map to any notch locus and align one at a time.
# look for matching known SUN positions to determine paralog then phase based on variants.