"""
Accessory code for analyzing pickled counts datasets from extract_linked_bams.py
"""
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import cPickle as pickle

base_dir = "/hive/users/ifiddes/longranger-1.2.0/notch2nl_10x/linked_bam_analysis"
sample = "NA12878"


def find_unique_gems(assigned_gems_by_para):
    gem_counter = Counter()
    for gems in assigned_gems_by_para.itervalues():
        for gem in gems:
            gem_counter[gem] += 1
    return {x for x, y in gem_counter.iteritems() if y == 1}


def assign_gems(counts):
    r = {"GemId": [] ,"SUN concordance rate": [], "SUN overlaps": [], "# reads in gem": []}
    total_gems = set()
    assigned_gems = set()
    assigned_gems_by_para = {x: set() for x in counts.keys()}
    for para, c in counts.iteritems():
        for gem_id, (num_support, total_intersections, num_reads) in c.iteritems():
            if total_intersections != 0:
                r["GemId"].append(gem_id)
                r["SUN concordance rate"].append(num_support)
                r["SUN overlaps"].append(total_intersections)
                r["# reads in gem"].append(num_reads)
                assigned_gems.add(gem_id)
                if num_support != 0:
                    assigned_gems_by_para[para].add(gem_id)
            total_gems.add(gem_id)
    return r, total_gems, assigned_gems, assigned_gems_by_para


def do_pairplots(counts, base_dir, sample):
    """
    Produces three pairplots - one for each group and a joint plot.
    """
    markers = ["o", "s"]
    r, total_gems, assigned_gems, assigned_gems_by_para = assign_gems(counts)
    df = pd.DataFrame.from_dict(r)
    unique_gems = find_unique_gems(assigned_gems_by_para)
    num_unique = len(unique_gems)
    num_not_unique = len(df) - num_unique
    unique_bins = ["{:,} unique".format(num_unique) if x in unique_gems else "{:,} not unique".format(num_not_unique) for x in df["GemId"]]
    df["Unique mappings"] = unique_bins
    sns_plot = sns.pairplot(df, hue="Unique mappings", markers=markers, plot_kws=dict(s=10))
    sns_plot.fig.text(0.87, 0.6, "{:,} Total Gems".format(len(total_gems)))
    sns_plot.savefig(os.path.join(base_dir, "{}_combined_plot.pdf".format(sample)), format="pdf")
    # now re-label to simply unique/not unique and make separate pairplots
    unique_simple_bins = ["Unique" if x in unique_gems else "Not Unique" for x in df["GemId"]]
    df["Unique mappings"] = unique_simple_bins
    for i, subset in enumerate(["Unique", "Not Unique"]):
        df2 = df[df["Unique mappings"] == subset]
        color = sns.color_palette()[i]
        cmap = sns.light_palette(color, as_cmap=True)
        sns_plot = sns.pairplot(df2, markers=markers[i], plot_kws=dict(color=color, s=10))
        sns_plot.map_lower(sns.kdeplot, cmap=cmap, n_levels=50)
        p = subset.replace(" ", "_").lower()
        sns_plot.savefig(os.path.join(base_dir, "{}_{}_combined_plot.pdf".format(sample, p)), format="pdf")
    plt.close('all')


def do_bin_plot(counts, base_dir, sample):
    r, total_gems, assigned_gems, assigned_gems_by_para = assign_gems(counts)
    unique_gems = find_unique_gems(assigned_gems_by_para)
    r = defaultdict(set)
    for para, c in counts.iteritems():
        for gem_id, (num_support, total_intersections, num_reads) in c.iteritems():
                if gem_id not in unique_gems and num_support > 0:
                    r[gem_id].add(para)
    r2 = Counter()
    for paras in r.itervalues():
        r2[frozenset(paras)] += 1
    data = {}
    for names, val in r2.iteritems():
        n = ",".join(sorted([x[-1] if "-" in x else x[0] for x in names]))
        data[n] = val
    colors = sns.color_palette("cubehelix", len(data))
    colors = dict([(k, color) for k, color in zip(*[data.keys(), colors])])
    cats = ["A", "B", "C", "D", "N"]
    for i, cat in enumerate(sorted(cats)):
        y = 0
        for key, val in data.items():
            if cat in key:
                plt.bar(i, val, 0.6, bottom=y, color=colors[key])
                if val > 10:
                    keys = key.split(",")
                    p = keys.index(cat)
                    del keys[p]
                    new_key = ",".join(keys)
                    plt.text(i + 0.25, y, ' '.join(new_key), color="black")
                y += val
    plt.xlabel("Paralog")
    plt.ylabel("Paralogs also having supported reads in gem")
    plt.title("GEMs that have one or more reads supporting more than one paralog")
    plt.xticks(np.arange(len(cats)) + 0.3, cats)
    plt.savefig(os.path.join(base_dir, "{}_groups_plot.pdf".format(sample)), format="pdf")
    plt.close('all')
