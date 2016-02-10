"""
Extract phase blocks from VCF
"""

import vcf
import sys
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('invcf', type=argparse.FileType('r'))
    parser.add_argument('outbed', type=argparse.FileType('w'))
    parser.add_argument('--startPos', default=118989614)
    parser.add_argument('--endPos', default=150469354)
    parser.add_argument('--chrom', default='chr1')
    parser.add_argument('--genome', default=None, help='add ID trackline')
    return parser.parse_args()


def get_phase_blocks(vcf_handle, chrom, start, end):
    rec_map = defaultdict(list)
    for rec in vcf_handle.fetch(chrom, start, end):
        rec_map[rec.samples[0]['PS']].append(rec)
    phase_blocks = []
    for ps, recs in rec_map.iteritems():
        if len(recs) > 1:
            phase_blocks.append(map(str, [chrom, recs[0].POS, recs[-1].POS, ps]))
    return phase_blocks


def main():
    args = parse_args()
    vcf_handle = vcf.Reader(args.invcf)
    phase_blocks = get_phase_blocks(vcf_handle, args.chrom, args.startPos, args.endPos)
    if args.Id is not None:
        args.outbed.write('track name="{} Phase Blocks"\n'.format(args.genome))
    for rec in phase_blocks:
        args.outbed.write('\t'.join(rec) + '\n')


if __name__ == '__main__':
    main()
