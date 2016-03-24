import sys
sys.path.append("/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio")
from pycbio.bio.intervals import gap_merge_intervals, ChromosomeInterval

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

