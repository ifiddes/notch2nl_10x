"""
Attempt to construct unitigs out of reads in the region and map them back to the reference
"""

roi = ['chr1', 145966000, 146510000]

na12878_bam = '/hive/users/ifiddes/longranger-1.2.0/NA12878_combined/outs/phased_possorted_bam.bam'

import pysam
import numpy as np
from collections import defaultdict
from extract_linked_bams import bin_reads

bam_handle = pysam.Samfile(na12878_bam)
unitig_region_mapped_reads = bin_reads(roi[0], roi[1], roi[2], bam_handle, offset=10000)


def eval_rec(rec, read_dict):
    tags = dict(rec.tags)
    if 'BX' in tags:
        t = tags['BX']
        if t in read_dict:
            read_dict[t].append(rec)


# now bring in the unmapped reads
unmapped_bam = '/hive/users/ifiddes/longranger-1.2.0/NA12878_combined/outs/unmapped.bam'
unitig_region_mapped_and_unmapped_reads = unitig_region_mapped_reads.copy()
unmapped_handle = pysam.Samfile(unmapped_bam)
for rec in unmapped_handle.fetch(until_eof=True):
    eval_rec(rec, unitig_region_mapped_and_unmapped_reads)


from pycbio.bio.bio import reverse_complement


def bam_to_rec(read):
    """Convert pysam record to fastq
    https://www.biostars.org/p/6970/
    """
    seq = read.seq
    qual = read.qual
    if read.is_reverse:
        seq = reverse_complement(seq)
        qual = qual[::-1]
    return "@{}\n{}\n+\n{}\n".format(read.qname, seq, qual)


reads = []
for tag in unitig_region_mapped_and_unmapped_reads:
    for read in unitig_region_mapped_and_unmapped_reads[tag]:
        reads.append(bam_to_rec(read))


with open("reads_to_assemble.fq", "w") as outf:
    for read in reads:
        outf.write(read)


# fermikit assembly
fermi.kit/fermi2.pl unitig -s1m -t20 -l98 -p nudtp4.mak \
" cat /hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/reads_to_assemble.fq | fermi.kit/trimadap-mt -p4 " > nudtp4.mak
make -f nudtp4.mak

# align fermikit assembly
bwa mem ~/ifiddes_hive/chr1_fastas/full_hg38_reference/hg38.fa nudtp4.mak.mag > nudtp4.mak.mag.sam

# also try SV calls
fermi.kit/run-calling -t20 ~/ifiddes_hive/chr1_fastas/full_hg38_reference/hg38.fa nudtp4.mag.gz | sh


convert_sort_sam () {
    n=${1%.*}
    samtools view -b $1 | samtools sort -O bam -T $TMPDIR/$n - > $n.sorted.bam
    samtools index $n.sorted.bam
}

