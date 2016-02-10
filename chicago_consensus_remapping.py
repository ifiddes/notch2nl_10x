"""
Mapping Chicago reads to my Notch2NL consensus
"""
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


for f in *bam; do
    samtools bamshuf -Ou $f tmp | samtools bam2fq -s /dev/null - > $f.pairs.fq
done
# toss out singles, since for chicago that means no phasing information
# remap to chr1 hg38
for f in *fq; do
    bwa mem -t 4 -p /hive/users/ifiddes/chr1_fastas/hg38_chr1_repeatmasked.fasta $f | samtools view -b - | samtools sort -T tmp -O bam - > $f.hg38.sorted.bam
done

# remap also to notch2nl consensus
for f in *fq; do
    bwa mem -t 10 -p /hive/users/ifiddes/notch2nl_suns/notch2_aligned_consensus.fasta $f | samtools view -b - | samtools sort -T tmp -O bam - > $f.notch2nl_consensus.sorted.bam
done

# rather than sorting by read name, just pull in all reads and make a dict
bam_handle = pysam.Samfile("/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/chicago_data/CP499_chicago_500kb_notch2nl.bam.pairs.fq.notch2nl_consensus.sorted.bam")
read_dict = defaultdict(list)
for read in bam_handle:
    if read.is_secondary:
        continue
    read_dict[read.qname].append(read)


# interesting read pairs: 1) match SNP and SUN, 2) match SNP/SUN and other read remains unmapped
# filter out more complicated alignments
filtered_read_dict = {"BothMapped": {}, "OneMapped": {}}
for name, reads in read_dict.iteritems():
    if len(reads) != 2:
        continue
    l, r = reads
    if (l.is_unmapped and not r.is_unmapped) or (not l.is_unmapped and r.is_unmapped):
        filtered_read_dict["OneMapped"][name] = [l, r]
    elif not(l.is_unmapped and r.is_unmapped):
        filtered_read_dict["BothMapped"][name] = [l, r]


# use code from remap_reads_consensus to find interesting read pairs

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


def find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
    sun_intersections = [[sun_recs[i], x] for i, x in enumerate(sun_intervals) if x.proper_subset(read_interval)]
    snp_intersections = [[snp_recs[i], x] for i, x in enumerate(snp_intervals) if x.proper_subset(read_interval)]
    sun_results = list(compare_sun_to_reference(sun_intersections, read))
    snp_results = list(compare_snp_to_reference(snp_intersections, read))
    return snp_results, sun_results


snp_h = vcf.Reader(open("../remapped_heterozygous_snps_consensus.vcf"))
snp_recs = [x for x in snp_h if not (x.is_indel is False and x.is_deletion is True)]
v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels_ConsensusRef.vcf.gz"))
sun_recs = [x for x in v_h if not (x.is_indel is False and x.is_deletion is True)]
sun_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in sun_recs]
snp_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in snp_recs]


analyzed_onemapped = {"snp": [], "sun": []}
read_holder = {}
for name, (read1, read2) in filtered_read_dict["OneMapped"].iteritems():
    read = read1 if not read1.is_unmapped else read2
    if read.mapq < 30:
        continue
    assert read1.is_unmapped or read2.is_unmapped
    snp, sun = find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle)
    if len(snp) > 0:
        analyzed_onemapped["snp"].append([name, snp])
    if len(sun) > 0:
        analyzed_onemapped["sun"].append([name, sun])
    if (len(snp) > 0 and len(sun) == 0) or (len(snp) == 0 and len(sun) > 0):
        read_holder[name] = read


# 1220 snp hits, 915 sun hits.
# read holder holds the 1913 pairs that are not overlaps
good_names = set(read_holder.viewkeys())
hg38_bam_handle = pysam.Samfile("/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/chicago_data/CP499_chicago_500kb_notch2nl.bam.pairs.fq.hg38.sorted.bam")
distant_reads = {}
for read in hg38_bam_handle.fetch('chr1', 119989614 - 550000, 149469354 + 550000):  # save computation time
    if not read.is_unmapped and read.qname in good_names and read.seq != read_holder[read.qname].seq:
        distant_reads[read.qname] = read


# 1871 reads mapped elsewhere.
# how many match phased SNPs based on 10x phasing?
hg38_snp_vcf = vcf.Reader(open("/hive/users/ifiddes/longranger-1.2.0/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz"))
hg38_snp_recs = hg38_snp_vcf.fetch('chr1', 119989614 - 550000, 149469354 + 550000)
hg38_snp_recs = [rec for rec in hg38_snp_recs if rec.samples[0].phased is True]


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
            yield "NoHaplotype",  vcf_rec, read
        else:
            yield gt_map[read_seq], vcf_rec, read


def find_read_snp_intersections(read, hg38_snp_intervals, hg38_snp_recs, hg38_bam_handle):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    read_interval = ChromosomeInterval(hg38_bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
    snp_intersections = [[hg38_snp_recs[i], x] for i, x in enumerate(hg38_snp_intervals) if x.proper_subset(read_interval)]
    snp_results = list(compare_snp_to_reference(snp_intersections, read))
    return snp_results


# now we bin these reads using the same logic as before
hg38_snp_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in hg38_snp_recs]
read_tag_holder = defaultdict(list)
for n, read in distant_reads.iteritems():
    snp_results = find_read_snp_intersections(read, hg38_snp_intervals, hg38_snp_recs, hg38_bam_handle)
    if len(snp_results) > 0:
        for h in snp_results:
            tag, vcf_rec, read = h
            read_tag_holder[tag].append([vcf_rec, read, read_holder[read.qname]])


# we now have two reciprocal sets of reads, each with 1 SUN hit and 1 SNP hit elsewhere.
# we need to visualize this. First, we need to re-establish the pairing information and remap the consensus read to hg38
# the challenge here is to maintain mismatch information - we can try IGV paired-read mode? or separate bins



# this didn't really work at all. My reads are mostly poorly mapped.
# I have another idea - simply loop over hg38 aligned reads and look for general SUN/SNP matches.

# i need to sort the BAM by read tag
bam_handle = pysam.Samfile("CP499_chicago_500kb_notch2nl.bam.pairs.fq.hg38.read_sorted.bam")
read_dict = defaultdict(list)
for read in bam_handle:
    read_dict[read.qname].append(read)


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
            yield gt_map[read_seq], vcf_rec, read


def find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle):
    """
    Given a set of binned reads, do they intersect with unique positions? if so, determine if they also match the
    unique base(s) at that position.
    """
    # todo: get insertions to work
    read_interval = ChromosomeInterval(bam_handle.getrname(read.tid), read.positions[0], read.positions[-1], None)
    sun_intersections = [[sun_recs[i], x] for i, x in enumerate(sun_intervals) if x.proper_subset(read_interval)]
    snp_intersections = [[snp_recs[i], x] for i, x in enumerate(snp_intervals) if x.proper_subset(read_interval)]
    sun_results = list(compare_sun_to_reference(sun_intersections, read))
    snp_results = list(compare_snp_to_reference(snp_intersections, read))
    return snp_results, sun_results


snp_h = vcf.Reader(open("/hive/users/ifiddes/longranger-1.2.0/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz"))
snp_recs = [x for x in snp_h.fetch('chr1', 119989614 - 550000, 149469354 + 550000) if not (x.is_indel is False and x.is_deletion is True)]
v_h = vcf.Reader(open("/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels.vcf.gz"))
sun_recs = [x for x in v_h if not (x.is_indel is False and x.is_deletion is True)]
sun_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in sun_recs]
snp_intervals = [ChromosomeInterval(x.CHROM, x.start - 1, x.end - 1, None) for x in snp_recs]


analyzed_reads = defaultdict(lambda: defaultdict(list))
for qname, reads in read_dict.iteritems():
    if len(reads) != 2 or any([x.is_unmapped for x in reads]):
        continue
    for read in reads:
        snp_results, sun_results = find_read_snp_sun_intersections(read, snp_intervals, snp_recs, sun_intervals, sun_recs, bam_handle)
        if len(snp_results) > 0:
            analyzed_reads[qname]['snp'].append(snp_results)
        if len(sun_results) > 0:
            analyzed_reads[qname]['sun'].append(sun_results)


