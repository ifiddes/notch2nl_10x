"""
Run the traditional WGS-SUN based pipeline on 10x data to compare to the results
"""

import pysam
import sys
import vcf
import string
import itertools
import numpy as np
import argparse
import tempfile
import os
import subprocess
from pyfasta import Fasta
from operator import itemgetter
from itertools import groupby
from collections import Counter, defaultdict
sys.path.append("/hive/users/ifiddes/pycbio")
from pycbio.sys.procOps import runProc, callProc
from pycbio.sys.fileOps import tmpFileGet
from pycbio.sys.mathOps import format_ratio
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('inBam', help='(10x) bamfile to remap')
    parser.add_argument('outPdf', help='path to write plot to')
    parser.add_argument('--outBam', default=None, help='path to write consensus aligned bam to')
    parser.add_argument('--consensusVcf', default='/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels_ConsensusRef.vcf.gz')
    parser.add_argument('--consensusRef', default='/hive/users/ifiddes/notch2nl_suns/notch2_aligned_consensus.fasta')
    return parser.parse_args()


regions = [['chr1', 119990189, 120163923, 'Notch2'],
           ['chr1', 146149601, 146329308, 'Notch2NL-A'],
           ['chr1', 148597945, 148786127, 'Notch2NL-B'],
           ['chr1', 149349476, 149477855, 'Notch2NL-C'],
           ['chr1', 120706154, 120801963, 'Notch2NL-D']]


def extract_reads(bam, offset=25000):
    tmp_paired = tmpFileGet(suffix='paired.fq')
    tmp_single = tmpFileGet(suffix='single.fq')
    tmp_shuf = tmpFileGet()
    region_strs = ['{}:{}-{}'.format(chrom, start - offset, stop + offset) for chrom, start, stop, para in regions]
    view_cmd = ['samtools', 'view', '-b', bam]
    view_cmd.extend(region_strs)
    cmd = [view_cmd,
           ['samtools', 'bamshuf', '-Ou', '-', tmp_shuf],
           ['samtools', 'bam2fq', '-s', tmp_single, '-']]
    with open(tmp_paired, 'w') as tmp_paired_h:
        runProc(cmd, stdout=tmp_paired_h)
    return tmp_paired, tmp_single


def remap_reads(tmp_paired, tmp_single, index):
    sort_tmp = tmpFileGet()
    paired_bam = tmpFileGet(suffix='paired.sorted.bam')
    unpaired_bam = tmpFileGet(suffix='unpaired.sorted.bam')
    paired_cmd = [['bwa', 'mem', '-p', index, tmp_paired],
                  ['samtools', 'view', '-b', '-'],
                  ['samtools', 'sort', '-T', sort_tmp, '-O', 'bam', '-']]
    unpaired_cmd = [['bwa', 'mem', index, tmp_single],
                   ['samtools', 'view', '-b', '-'],
                   ['samtools', 'sort', '-T', sort_tmp, '-O', 'bam', '-']]
    for cmd, path in [[paired_cmd, paired_bam], [unpaired_cmd, unpaired_bam]]:
        with open(path, 'w') as f_h:
            runProc(cmd, stdout=f_h)
    return paired_bam, unpaired_bam


def merge_bams(paired_bam, unpaired_bam, out_bam):
    cmd = ['samtools', 'merge', out_bam, paired_bam, unpaired_bam]
    runProc(cmd)
    cmd = ['samtools', 'index', out_bam]
    runProc(cmd)


def build_remapped_bam(in_bam, consensus_ref, out_bam):
    tmp_paired, tmp_single = extract_reads(in_bam)
    paired_bam, unpaired_bam = remap_reads(tmp_paired, tmp_single, consensus_ref)
    merge_bams(paired_bam, unpaired_bam, out_bam)
    for p in [tmp_paired, tmp_single, paired_bam, unpaired_bam]:
        os.remove(p)


def pileup(out_bam, vcf_path):
    bases = {"A", "T", "G", "C", "a", "t", "g", "c"}
    vcf_handle = vcf.Reader(open(vcf_path))
    wgs_results = defaultdict(list)
    for vcf_rec in vcf_handle:
        if vcf_rec.is_indel:
            continue
        pos_str = "{0}:{1}-{1}".format(vcf_rec.CHROM, vcf_rec.POS)
        cmd = ['samtools', 'mpileup', '-q', '20', '-Q', '20', '-r', pos_str, out_bam]
        mpileup_rec = callProc(cmd).split()
        pile_up_result = Counter(x.upper() for x in mpileup_rec[4] if x in bases)
        sample_dict = {s.sample: s.gt_bases for s in vcf_rec.samples}
        for s in vcf_rec.samples:
            if len([x for x in sample_dict.itervalues() if x == s.gt_bases]) != 1:
                continue
            c = 1.0 * pile_up_result[s.gt_bases] / len(mpileup_rec[4])
            wgs_results[s.sample].append([vcf_rec.POS, c])
    return wgs_results


def plot_results(wgs_results, out_pdf, aln_size):
    paralogs = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D']
    fig, plots = plt.subplots(5, sharey=True, sharex=True)
    plt.yticks((0, 0.1, 0.2, 0.3, 0.4))
    plt.ylim((0, 0.4))
    xticks = range(0, int(round(aln_size / 10000.0) * 10000.0))
    plt.xticks(xticks)
    plt.xlim((0, aln_size))
    plt.xlabel("Alignment position")
    for i, (p, para) in enumerate(zip(plots, paralogs)):
        p.set_title(para)
        wgs = wgs_results[para]
        xvals, yvals = zip(*wgs)
        p.vlines(xvals, np.zeros(len(xvals)), yvals, color=sns.color_palette()[0], alpha=0.7, linewidth=0.8)
        # mark the zeros
        zero_wgs = [[x, y + 0.02] for x, y in wgs if y == 0]
        if len(zero_wgs) > 0:
            z_xvals, z_yvals = zip(*zero_wgs)
            p.vlines(z_xvals, np.zeros(len(z_xvals)), z_yvals, color=sns.color_palette()[2], alpha=0.7, linewidth=0.8)
    plt.tight_layout(pad=2.5, h_pad=0.25)
    zero_line = matplotlib.lines.Line2D([], [], color=sns.color_palette()[2])
    reg_line = matplotlib.lines.Line2D([], [], color=sns.color_palette()[0])
    fig.legend(handles=(reg_line, zero_line), labels=["WGS SUN Fraction", "WGS Missing SUN"], loc="upper right")
    fig.text(0.01, 0.5, 'SUN fraction of reads', va='center', rotation='vertical')
    plt.savefig(out_pdf, format="pdf")
    plt.close()


def get_aln_size(consensus_ref):
    f = Fasta(consensus_ref)
    assert len(f) == 1
    return len(f[f.keys()[0]])


def main():
    args = parse_args()
    if args.outBam is None:
        out_bam = tmpFileGet(suffix='merged.sorted.bam')
    else:
        out_bam = args.outBam
    build_remapped_bam(args.inBam, args.consensusRef, out_bam)
    wgs_results = pileup(out_bam, args.consensusVcf)
    aln_size = get_aln_size(args.consensusRef)
    plot_results(wgs_results, args.outPdf, aln_size)
    if args.outBam is None:
        os.remove(out_bam)


if __name__ == '__main__':
    main()
