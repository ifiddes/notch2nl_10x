bamfile = "/hive/users/ifiddes/longranger/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam"
vcf_path = "/hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels.vcf.gz"
region_bed = "/cluster/home/ifiddes/public_html/notch_beds/notch2nl_hg38_no_trackname.bed"
bam_handle = pysam.Samfile(bamfile)
v_h = vcf.Reader(open(vcf_path))
regions = list(parse_bed(region_bed))



# for each gem that joins A and B, build a BED track showing each as a separate color. See if you can find haplotypes this way.

gems = defaultdict(set)
for para, c in counts.iteritems():
    if para not in ["Notch2NL-A", "Notch2NL-B"]:
        continue
    for gem_id, (num_support, total_intersections, num_reads) in c.iteritems():
        if num_support > 0:
            gems[para].add(gem_id)


a_and_b = gems["Notch2NL-A"] & gems["Notch2NL-B"]

# need to reverse the alignment pos - seq pos thing as before
start_pos = {"Notch2": 119989614, "Notch2NL-A": 146149503, "Notch2NL-B": 148600445, "Notch2NL-C": 149469354, "Notch2NL-D": 120799019}
# which of these are backwards?
backwards = ["Notch2NL-A", "Notch2NL-B", "Notch2"]



# build a map of alignment positions to sequence positions
r = {name: seq for name, seq in fastaRead("/hive/users/ifiddes/notch2nl_suns/notch2_aligned.fasta")}
r_sort = sorted(r.iteritems(),key=lambda x: x[0])
names, seqs = zip(*r_sort)
tgt_is = {n: 0 for n in names if n not in backwards}
for n in backwards:
    tgt_is[n] = len(r[n].replace("-", "")) + 1


pos_map = defaultdict(dict)
for ref_i, cs in enumerate(zip(*seqs)):
    for name, tgt_i in tgt_is.iteritems():
        pos_map[name][ref_i] = tgt_i
    for name, c in zip(*[names, cs]):
        if c != "-":
            if name not in backwards:
                tgt_is[name] += 1
            else:
                tgt_is[name] -= 1


# map interesting reads
read_holder = defaultdict(list)
for para in ["Notch2NL-A", "Notch2NL-B"]:
    for gem in a_and_b:
        for sup_read in interesting_read_holder[para][gem][0]:
            read_pos, ref_pos = zip(*sup_read.aligned_pairs)
            ref_pos = [x for x in ref_pos if x is not None]
            left = pos_map[para][min(ref_pos) - start_pos[para]]
            right = pos_map[para][max(ref_pos) - start_pos[para]]
            sup_read = [para, left, right]
            read_holder[gem].append(sup_read)
        for read in interesting_read_holder[para][gem][1]:
            if read.is_unmapped:
                continue
            read_pos, ref_pos = zip(*read.aligned_pairs)
            ref_pos = [x for x in ref_pos if x is not None]
            try:
                left = pos_map[para][min(ref_pos) - start_pos[para]]
                right = pos_map[para][max(ref_pos) - start_pos[para]]
            except:
                continue  # outside plot window anyways
            read = ["None", left, right]
            read_holder[gem].append(read)


sorted_read_holder = OrderedDict()
for gem, reads in read_holder.iteritems():
    reads = sorted(reads, key=lambda x: x[1])
    sorted_read_holder[gem] = reads



colors = {"Notch2NL-A": "#90C3D4", "Notch2NL-B": "#D490C3", "None": "grey"}
alphas = {"Notch2NL-A": 1.0, "Notch2NL-B": 1.0, "None": 0.5}
heights = {"Notch2NL-A": 1.0, "Notch2NL-B": 1.0, "None": 0.2}
bottom = {"Notch2NL-A": 0, "Notch2NL-B": 0, "None": 0.4}


import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
sns.set_style("white")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('{}_disco_haplotypes.pdf'.format(sample)) as pdf:
    for j in xrange(0, len(read_holder), 32):
        tmp_reads = OrderedDict(list(sorted_read_holder.items()[j: j + 32]))
        fig, plots = plt.subplots(min(32, len(read_holder) - j), sharey=True, sharex=True)
        plt.yticks([])
        plt.ylim((0, 1))
        plt.xticks((0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))
        plt.xlim((0, 100143))
        plt.xlabel("Notch2NL Consensus Position")
        for i, (p, gem) in enumerate(zip(plots, tmp_reads)):
            p.hlines(color="black", y=0.5, xmin=0, xmax=100143, linewidth=0.5)
            for para, left, right in tmp_reads[gem]:
                box = {"xy": (left, bottom[para]), "width": right - left, "height": heights[para], "color": colors[para], "alpha": alphas[para]}
                p.add_patch(patches.Rectangle(**box))
            if i == len(read_holder) - 1:
                sns.despine(top=True, left=True, right=True, ax=p)
            else:
                sns.despine(top=True, left=True, right=True, bottom=True, ax=p)
        if j == 0:
            plt.suptitle("Visualization of linked reads in {} gems which support both paralogs".format(len(read_holder)))
        a_patch = patches.Patch(color=colors["Notch2NL-A"], label="Notch2NL-A")
        b_patch = patches.Patch(color=colors["Notch2NL-B"], label="Notch2NL-B")
        fig.legend(handles=[a_patch, b_patch], labels=["Notch2NL-A", "Notch2NL-B"])
        plt.tight_layout(pad=2.5, h_pad=0.2)
        pdf.savefig()


# now build a version showing gems that support a specific paralog above a certain level
read_holder = {"Notch2NL-A": defaultdict(list), "Notch2NL-B": defaultdict(list)}
for para in ["Notch2NL-A", "Notch2NL-B"]:
    for gem in gems[para]:
        if gem in a_and_b:
            continue
        for sup_read in interesting_read_holder[para][gem][0]:
            read_pos, ref_pos = zip(*sup_read.aligned_pairs)
            ref_pos = [x for x in ref_pos if x is not None]
            left = pos_map[para][min(ref_pos) - start_pos[para]]
            right = pos_map[para][max(ref_pos) - start_pos[para]]
            sup_read = [para, left, right]
            read_holder[para][gem].append(sup_read)
        for read in interesting_read_holder[para][gem][1]:
            if read.is_unmapped:
                continue
            read_pos, ref_pos = zip(*read.aligned_pairs)
            ref_pos = [x for x in ref_pos if x is not None]
            try:
                left = pos_map[para][min(ref_pos) - start_pos[para]]
                right = pos_map[para][max(ref_pos) - start_pos[para]]
            except:
                continue  # outside plot window anyways
            read = ["None", left, right]
            read_holder[para][gem].append(read)



sorted_read_holder = OrderedDict((("Notch2NL-A", defaultdict(list)), ("Notch2NL-B",defaultdict(list))))
for para in ["Notch2NL-A", "Notch2NL-B"]:
    for gem, reads in read_holder[para].iteritems():
        reads = sorted(reads, key=lambda x: x[1])
        sup = Counter([x[0] for x in reads])
        if sup[para] > 3:
            sorted_read_holder[para][gem] = reads



for para in ["Notch2NL-A", "Notch2NL-B"]:
    tmp_read_holder = sorted_read_holder[para]
    with PdfPages('{}_{}_haplotypes.pdf'.format(sample, para)) as pdf:
        for j in xrange(0, len(tmp_read_holder), 23):
            tmp_reads = OrderedDict(list(tmp_read_holder.items()[j: j + 23]))
            fig, plots = plt.subplots(min(23, len(read_holder) - j), sharey=True, sharex=True)
            plt.yticks([])
            plt.ylim((0, 1))
            plt.xticks((0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000))
            plt.xlim((0, 100143))
            plt.xlabel("Notch2NL Consensus Position")
            for i, (p, gem) in enumerate(zip(plots, tmp_reads)):
                p.hlines(color="black", y=0.5, xmin=0, xmax=81507, linewidth=0.5)
                for t, left, right in tmp_reads[gem]:
                    box = {"xy": (left, bottom[t]), "width": right - left, "height": heights[t], "color": colors[t], "alpha": alphas[t]}
                    p.add_patch(patches.Rectangle(**box))
                if i == len(tmp_read_holder) - 1:
                    sns.despine(top=True, left=True, right=True, ax=p)
                else:
                    sns.despine(top=True, left=True, right=True, bottom=True, ax=p)
            if j == 0:
                plt.suptitle("Visualization of linked reads in {} gems which support only {}\nwith support level > 3".format(len(tmp_read_holder), para))
            plt.tight_layout(pad=2.5, h_pad=0.2)
            pdf.savefig()
    plt.close('all')