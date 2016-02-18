"""
I totally fucked this up and deleted the code. This is the salvage attempt from the log.
"""
import sys, os, vcf, string
from collections import defaultdict, Counter, OrderedDict
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from sonLib.bioio import fastaRead
from lib.seq_lib import ChromosomeInterval
sys.path.append("/hive/users/ifiddes/test_pipeline/pipeline/submodules/pycbio")
from pycbio.sys.procOps import callProcLines
from pyfasta import Fasta

def rev_comp(s, rc=string.maketrans("ATGC", "TACG")):
    return string.translate(s, rc)[::-1]


# We need the start positions for each in the actual genome, used a browser BLAT
start_pos = {"Notch2": 120163697, "Notch2NL-A": 146328264, "Notch2NL-B": 148779287, "Notch2NL-C": 149350805, "Notch2NL-D": 120705667}
backwards = {"Notch2NL-C", "Notch2NL-D"}
names = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D'] # same as in VCF
header = "##fileformat=VCFv4.1"
fields = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + names)
rec_template = "Notch2\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"


# dumb map
m = {"Notch2NL-C_Notch2NL-D": ("Notch2", "Notch2NL-A", "Notch2NL-B"),
    "Notch2NL-D": ("Notch2", "Notch2NL-A", "Notch2NL-B", "Notch2NL-C"),
    "all": ("Notch2", "Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D")}

regions = {frozenset(["Notch2NL-D"]): [[0, 15866], [74917, 81068], [162369, 165396]],
           frozenset(["Notch2NL-D", "Notch2NL-C"]): [[15867, 74916]],
           frozenset(): [[81069, 162368], [165397, 2000000]]}



f = Fasta("stitched_alignment.fa")
results = {}
for exclude in [frozenset(), frozenset(["Notch2NL-D"]), frozenset(["Notch2NL-D", "Notch2NL-C"])]:
    t = open("tmp.fasta", "w")
    for para in sorted(set(f.keys()) - exclude):
        t.write(">{}\n{}\n".format(para, f[para]))
    t.close()
    n = '_'.join(sorted(exclude)) if len(exclude) > 0 else 'all'
    cmd = ['java', '-jar', '/cluster/home/ifiddes/jvarkit/dist-1.133/biostar94573.jar', '-R', n,
           'tmp.fasta']
    r = callProcLines(cmd)
    recs = [x.split() for x in r if not x.startswith("#")]
    results[exclude] = recs


raw_recs = []
for exclude, region in regions.iteritems():
   for start, stop in region:
    raw_recs.extend([x for x in results[exclude] if start < int(x[1]) <= stop])


# region with poor alignment
exclude_regions = [[28574, 31093]]
exclude_regions = [ChromosomeInterval('a', x[0], x[1], '.') for x in exclude_regions]
recs = []
for r in raw_recs:
    i = ChromosomeInterval('a', int(r[1]), int(r[1]) + 1, '.')
    if not any([i.overlap(x) for x in exclude_regions]):
        recs.append(r)

# build a map of alignment positions to sequence positions
r = {name: seq for name, seq in fastaRead("stitched_alignment.fa")}
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


rec_template = "Notch2\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
consensus_rec_template = "Notch2NL_extended_consensus\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
notch2_ref_recs = []
consensus_ref_recs = []
consensus_fa = Fasta("consensus.fa")["Notch2NL_extended_consensus"]
for rec in sorted(recs, key=lambda x: int(x[1])):
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    these_names = m[rec[0]]
    i = 9
    for name in m['all']:
        if name in these_names:
            h = rec[i].split("/")[0]
            i += 1
            try:
                h = int(h)
            except:
                h = "."
        else:
            h = "."
        if h != ".":
            gts[name] = combined_gts[h]
        else:
            gts[name] = "."
    # remove non-unique variants here
    q = {x: y for x, y in gts.iteritems() if y != "."}
    if 1 not in Counter(q.itervalues()).values():
        continue
    # remove ridiculously long indels
    if any([len(x) > 20 for x in gts.itervalues()]):
        continue
    if gts["Notch2"] == ref:
        # no need to fix anything, just write
        this_rec = map(str, [combined_gts.index(x) for x in gts.itervalues() if x != "."])
        this_ref = ref
        this_alt = rec[4]
    else:
        # we need to make this Notch2 reference
        new_gts = [gts["Notch2"]] + list(set(gts.itervalues()) - set([gts["Notch2"]]))
        new_gts = [x for x in new_gts if x != "."]
        this_rec = map(str, [new_gts.index(x) if x != "." else "." for x in gts.itervalues()])
        this_ref = new_gts[0]
        this_alt = ",".join(new_gts[1:])
    # do consensus ref
    con_ref = consensus_fa[aln_pos - 1: aln_pos - 1 + len(ref)].upper()
    all_gts = set(gts.values())
    con_gts = [con_ref] + list(all_gts - set([con_ref]))
    con_gts = [x for x in con_gts if x != "."]
    con_rec = map(str, [con_gts.index(x) if x != "." else "." for x in gts.itervalues()])
    con_alt = ",".join(con_gts[1:])
    notch2_ref_recs.append(rec_template.format(pos=aln_pos, ref=this_ref, alt=this_alt, gts="\t".join(this_rec)))
    consensus_ref_recs.append(consensus_rec_template.format(pos=aln_pos, ref=con_ref, alt=con_alt, gts="\t".join(con_rec)))


with open("filtered_records_extended_alignment_reference.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(notch2_ref_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)



with open("Notch2NL_SUN_UniqueIndels_ExtendedConsensusRef.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(consensus_ref_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)


chr1_fasta = Fasta("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.fa")['chr1']
# now we need to translate this to each paralog, creating a separate VCF record from the perspective of that paralog
filt_recs = [x.split() for x in open("filtered_records_extended_alignment_reference.vcf") if not x.startswith("#")]
rec_template = "chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n"
new_recs = []
for rec in filt_recs:
    aln_pos = int(rec[1])
    ref = rec[3]
    alts = rec[4].split(",")
    combined_gts = [ref] + alts
    gts = OrderedDict()
    for i, name in enumerate(names, 9):
        try:
            h = int(rec[i].split("/")[0])
            gts[name] = combined_gts[h]
        except:
            gts[name] = "."
    counts = Counter(gts.itervalues())
    for name in names:
        if counts[gts[name]] == 1 and gts[name] != '.':
            # this is unique to this paralog
            new_gts = [gts[name]] + list(set(gts.itervalues()) - set([gts[name]]))
            new_gts = [x for x in new_gts if x != "."]
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
            this_rec = map(str, [new_gts.index(x) if x != "." else "." for x in gts.itervalues()])
            new_recs.append(rec_template.format(pos=pos + 1, ref=this_ref, alt=this_alt, gts="\t".join(this_rec)))


with open("Notch2NL_SUN_UniqueIndels_extended.vcf", "w") as outf:
    outf.write(header + "\n")
    outf.write(fields + "\n")
    for r in sorted(new_recs, key=lambda x: int(x.split()[1])):
        outf.write(r)
