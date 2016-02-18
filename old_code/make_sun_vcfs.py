"""
Reformat jvarkit SUN calls
java -jar ~/jvarkit/dist-1.133/biostar94573.jar -R Notch2 notch2_aligned.fasta  > aligned_notch2_raw.vcf
"""
import sys, os, vcf, string
from collections import defaultdict, Counter, OrderedDict
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from sonLib.bioio import fastaRead
from lib.seq_lib import ChromosomeInterval


def rev_comp(s, rc=string.maketrans("ATGC", "TACG")):
    return string.translate(s, rc)[::-1]


# We need the start positions for each in the actual genome, used a browser BLAT
start_pos = {"Notch2": 120087516, "Notch2NL-A": 146248223, "Notch2NL-B": 148698969, "Notch2NL-C": 149374496, "Notch2NL-D": 120707775}
# which of these are backwards?
backwards = {"Notch2NL-C", "Notch2NL-D"}

names = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D'] # same as in VCF
header = "##fileformat=VCFv4.1"
fields = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)
rec_template = "Notch2\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"


# build a map of alignment positions to sequence positions
r = {name: seq for name, seq in fastaRead("notch2_aligned.fasta")}
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


# now lets restructure the vcf to always make Notch2 the reference allele
# as well as turning it into a homozygous call, and removing the depth.
# finally, we filter for unique calls
recs = [x.split() for x in open("aligned_notch2_raw.vcf") if not x.startswith("#")]
# filtered indels file has no header  - I verified that the names are indeed in sort order as above
with open("filtered_records_alignment_reference.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for rec in sorted(recs, key=lambda x: int(x[1])):
        aln_pos = rec[1]
        ref = rec[3]
        alts = rec[4].split(",")
        combined_gts = [ref] + alts
        gts = OrderedDict()
        for i, name in enumerate(names, 9):
            h = rec[i].split("/")[0]
            if h == ".":
                # skip unaligned regions
                continue
            else:
                h = int(h)
            gts[name] = combined_gts[h]
        # remove non-unique variants here
        if 1 not in Counter(gts.itervalues()).values():
            continue
        # remove ridiculously long indels
        if any([len(x) > 20 for x in gts.itervalues()]):
            continue
        if gts["Notch2"] == ref:
            # no need to fix anything, just write
            this_rec = map(str, [combined_gts.index(x) for x in gts.itervalues()])
            this_ref = ref
            this_alt = rec[4]
        else:
            # we need to make this Notch2 reference
            new_gts = [gts["Notch2"]] + list(set(gts.itervalues()) - set([gts["Notch2"]]))
            this_rec = map(str, [new_gts.index(x) for x in gts.itervalues()])
            this_ref = new_gts[0]
            this_alt = ",".join(new_gts[1:])
        outf.write(rec_template.format(pos=aln_pos, ref=this_ref, alt=this_alt, gts="\t".join(this_rec)))


from pyfasta import Fasta
chr1_fasta = Fasta("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.fa")['chr1']

# now we need to translate this to each paralog, creating a separate VCF record from the perspective of that paralog
recs = [x.split() for x in open("filtered_records_alignment_reference.vcf") if not x.startswith("#")]
rec_template = "chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
new_recs = []
for rec in recs:
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    for i, name in enumerate(names, 9):
        h = int(rec[i].split("/")[0])
        gts[name] = combined_gts[h]
    counts = Counter(gts.itervalues())
    for name in names:
        if counts[gts[name]] == 1:
            # this is unique to this paralog
            new_gts = [gts[name]] + list(set(gts.itervalues()) - set([gts[name]]))
            if name not in backwards:
                pos = start_pos[name] - pos_map[name][aln_pos]
                this_ref = rev_comp(new_gts[0])
                this_alt = ",".join([rev_comp(x) for x in new_gts[1:]])
                if len(this_ref) > 1:
                    pos -= len(this_ref) - 1
            else:
                pos = pos_map[name][aln_pos] + start_pos[name]
                this_ref = new_gts[0]
                this_alt = ",".join(new_gts[1:])
            assert chr1_fasta[pos: pos + len(this_ref)].upper() == this_ref
            this_rec = map(str, [new_gts.index(x) for x in gts.itervalues()])
            new_recs.append(rec_template.format(pos=pos + 1, ref=this_ref, alt=this_alt, gts="\t".join(this_rec)))


with open("Notch2NL_SUN_UniqueIndels.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(new_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)


# make a VCF based off the consensus sequence
recs = [x.split() for x in open("aligned_notch2_raw.vcf") if not x.startswith("#")]
consensus = list(fastaRead("notch2_aligned_consensus.fasta"))[0][1]
rec_template = "Notch2NL_consensus\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
new_recs = []
for rec in sorted(recs, key=lambda x: int(x[1])):
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    for i, name in enumerate(names, 9):
        h = rec[i].split("/")[0]
        if h == ".":
            # skip unaligned regions
            continue
        else:
            h = int(h)
        gts[name] = combined_gts[h]
    # remove non-unique indels here
    if 1 not in Counter(gts.itervalues()).values():
        continue
    # remove ridiculously long indels
    if any([len(x) > 20 for x in gts.itervalues()]):
        continue
    consensus_ref = consensus[aln_pos - 1: aln_pos - 1 + len(ref)]
    if consensus_ref == ref:
        # no need to fix anything, just write
        this_rec = map(str, [combined_gts.index(x) for x in gts.itervalues()])
        this_ref = ref
        this_alt = rec[4]
    else:
        # we need to make this consensus reference
        all_gts = set(gts.values())
        new_gts = [consensus_ref] + list(all_gts - set([consensus_ref]))
        this_rec = map(str, [new_gts.index(x) for x in gts.itervalues()])
        this_ref = new_gts[0]
        this_alt = ",".join(new_gts[1:])
    new_recs.append(rec_template.format(pos=aln_pos, ref=this_ref, alt=this_alt, gts="\t".join(this_rec)))


with open("Notch2NL_SUN_UniqueIndels_ConsensusRef.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(new_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)


#repeat the process, but without offsetting, in BED format
recs = [x.split() for x in open("filtered_records_alignment_reference.vcf") if not x.startswith("#")]
rec_template = "Notch2NL_consensus\t{pos}\t{pos2}\t{ref}_{name}\n"
new_recs = []
for rec in recs:
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    for i, name in enumerate(names, 9):
        h = int(rec[i].split("/")[0])
        gts[name] = combined_gts[h]
    counts = Counter(gts.itervalues())
    for name in names:
        if counts[gts[name]] == 1:
            new_gts = [gts[name]] + list(set(gts.itervalues()) - set([gts[name]]))
            if name not in backwards:
                this_ref = rev_comp(new_gts[0])
            else:
                this_ref = new_gts[0]
            this_alt = ",".join(new_gts[1:])
            this_rec = map(str, [new_gts.index(x) for x in gts.itervalues()])
            new_recs.append(rec_template.format(pos=aln_pos - len(this_ref), pos2=aln_pos, ref=this_ref, name=name))


with open("Notch2NL_SUN_UniqueIndels.bed", "w") as outf:
    for r in sorted(new_recs, key=lambda x: (x.split()[0], int(x.split()[1]))):
        outf.write(r)



# in this version, a VCF record is produced for every SUN on every target.
recs = [x.split() for x in open("filtered_records_alignment_reference.vcf") if not x.startswith("#")]
rec_template = "chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
new_recs = []
for rec in recs:
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    for i, name in enumerate(names, 9):
        h = int(rec[i].split("/")[0])
        gts[name] = combined_gts[h]
    counts = Counter(gts.itervalues())
    if 1 not in counts.values():
        continue
    for ref_para in names:
        new_gts = [gts[ref_para]] + list(set(gts.itervalues()) - set([gts[ref_para]]))
        if ref_para not in backwards:
            pos = start_pos[ref_para] - pos_map[ref_para][aln_pos]
            this_ref = rev_comp(new_gts[0])
            this_alt = ",".join([rev_comp(x) for x in new_gts[1:]])
            if len(this_ref) > 1:
                pos -= len(this_ref) - 1
        else:
            pos = pos_map[ref_para][aln_pos] + start_pos[ref_para]
            this_ref = new_gts[0]
            this_alt = ",".join(new_gts[1:])
        assert chr1_fasta[pos: pos + len(this_ref)].upper() == this_ref
        this_gt = map(str, [new_gts.index(x) for x in gts.itervalues()])
        new_recs.append(rec_template.format(pos=pos + 1, ref=this_ref, alt=this_alt, gts="\t".join(this_gt)))



with open("Notch2NL_SUN_UniqueIndels_AllPerspectives.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(new_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)
