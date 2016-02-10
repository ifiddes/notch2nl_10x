"""
Functions to build mapping of consensus positions to genome positions
"""

import sys
from collections import defaultdict
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
from sonLib.bioio import fastaRead


# We need the start positions for each in the actual genome, used a browser BLAT
start_pos = {"Notch2": 120087516, "Notch2NL-A": 146248223, "Notch2NL-B": 148698969, "Notch2NL-C": 149374496, "Notch2NL-D": 120707775}
# which of these are backwards?
backwards = {"Notch2NL-C", "Notch2NL-D"}
names = ['Notch2', 'Notch2NL-A', 'Notch2NL-B', 'Notch2NL-C', 'Notch2NL-D'] # same as in VCF


def build_pos_map():
# build a map of alignment positions to sequence positions
    r = {name: seq for name, seq in fastaRead("/hive/users/ifiddes/notch2nl_suns/notch2_aligned.fasta")}
    r_sort = sorted(r.iteritems(), key=lambda x: x[0])
    names, seqs = zip(*r_sort)
    tgt_is = {n: 0 for n in names}
    pos_map = defaultdict(dict)
    for ref_i, cs in enumerate(zip(*seqs)):
        for name, tgt_i in tgt_is.iteritems():
            pos_map[name][ref_i] = tgt_i
        for name, c in zip(*[names, cs]):
            if c != "-":
                tgt_is[name] += 1
    return pos_map


def build_inverted_pos_map(pos_map=None):
    if pos_map is None:
        pos_map = build_pos_map()
    # invert pos_map
    pos_map_inverted = defaultdict(dict)
    for para, vals in pos_map.iteritems():
        if para not in backwards:
            vals = {start_pos[para] - y: x for x, y in vals.iteritems()}
        else:
            vals = {start_pos[para] + y: x for x, y in vals.iteritems()}
        pos_map_inverted[para] = vals
    return pos_map_inverted