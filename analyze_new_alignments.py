import vcf
import subprocess
from collections import defaultdict
import pysam

loci = [x.split()[:3] for x in open("/hive/users/ifiddes/notch2nl_suns/Notch2NL.bed")]
loci = [[x[0], int(x[1]), int(x[2])] for x in loci[:-1]]
loci = sorted(loci, key=lambda x: x[1])


for g in ["NA12878", "H9"]:
    r = pysam.Samfile("/hive/users/ifiddes/longranger-1.2.0/{}/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam".format(g))
    ps_map = defaultdict(list)
    for x, y in zip(*[loci, ["Notch2", "Notch2NL-D", "Notch2NL-A", "Notch2NL-B", "Notch2NL-C"]]):
        x[1] = x[-1] - 100000
        x[2] = x[2] + 100000
        for read in r.fetch(*x):
            tags = dict(read.tags)
            if "PS" not in tags:
                continue
            ps_map[(tags["PS"], tags["HP"])].append(read)
    for tag, reads in ps_map.iteritems():
        assert len(reads) > 0
        tag = "_".join(map(str, tag))
        with pysam.Samfile(os.path.join(g, tag + ".bam"), "wb", template=r) as outf:
            for read in reads:
                outf.write(read)


# create bigwigs for each
find . | grep bam$ | xargs -I{} -n 1 -P 10 sh -c "samtools sort -O bam -T '{}'.tmp '{}' > '{}.sorted'"
find . | grep bam.sorted$ | xargs -I{} -n 1 -P 10 sh -c "genomeCoverageBed -ibam '{}' -g /hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.fa -bg > '{}.bedGraph'"
for f in /bin/ls `pwd`/*/*_1.bam.sorted.bedGraph; do
cat $f | awk -F "\t" '{printf("%s\t%s\t%s\t-%s\n",$1,$2,$3,$4)}' > $f.negative
mv -f $f.negative $f
done
# merge all of same haplotype
for g in ["NA12878", "H9"]:
    wigs = {x.split(".")[0]: os.path.join(g, x) for x in os.listdir(g) if x.endswith("bedGraph")}
    combined_wigs = defaultdict(list)
    for x, y in wigs.iteritems():
        combined_wigs[x.split("_")[1]].append(y)
    for key, vals in combined_wigs.iteritems():
        subprocess.call("cat {} > {}".format(" ".join(vals), os.path.join(g, key + ".combined.bedGraph")), shell=True)

    
find . | grep .combined.bedGraph$ | xargs -n 1 -P 10 -I{} bedSort {} {}
find . | grep .combined.bedGraph$ | xargs -I{} -n 1 -P 10 sh -c "bedGraphToBigWig '{}' /hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.chrom.sizes '{}'.bw"


# add in SNP calls
zcat /hive/users/ifiddes/longranger-1.2.0/NA12878/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz | bedtools intersect -a stdin -b /hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels.vcf.gz > hub/hg38/NA12878.SNP_SUN.vcf

zcat /hive/users/ifiddes/longranger-1.2.0/H9/PHASER_SVCALLER_CS/PHASER_SVCALLER/_SNPINDEL_PHASER/ANALYZE_SNPINDEL_CALLS/fork0/files/varcalls.vcf.gz | bedtools intersect -a stdin -b /hive/users/ifiddes/notch2nl_suns/Notch2NL_SUN_UniqueIndels.vcf.gz > hub/hg38/H9.SNP_SUN.vcf

find . | grep vcf$ | xargs bgzip
find . | grep vcf.gz$ | xargs tabix -p vcf


# trackDb entries

base = """track {0}
container multiWig
shortLabel {0} Phased Coverage
longLabel {0} Phased Coverage
type bigWig

    track {0}_h_1
    parent {0}
    type bigWig
    bigDataUrl {0}.1.bw
    shortLabel {0} haplotype A
    longLabel {0} haplotype A
    yLineMark 0
    yLineOnOff on
    autoScale off
    visibility full
    viewLimits -50:50

    track {0}_h_2
    parent {0}
    type bigWig
    bigDataUrl {0}.2.bw
    shortLabel {0} haplotype B
    longLabel {0} haplotype B
    yLineMark 0
    yLineOnOff on
    autoScale off
    visibility full
    viewLimits -50:50

track {0}_vcf
type vcfTabix
bigDataUrl {0}.SNP_SUN.vcf.gz
shortLabel {0} SNP_SUNs
longLabel {0} SNP_SUNs

track {0}_full_vcf
type vcfTabix
bigDataUrl {0}.SNP.vcf.gz
shortLabel {0} SNPs
longLabel {0} SNPs

"""

with open("hub/hg38/trackDb.txt", "w") as outf:
    for g in ["NA12878", "H9"]:
        outf.write(base.format(g))


# why does phasing not span the whole locus? Why are PS numbers different between BAM and VCF?

# going to try local assembly using fermikit of phased reads from each haplotype separately

find . | grep sorted$ | xargs -n 1 -P 20 -I{} sh -c "samtools bam2fq '{}' > '{}'.fq"

for f in  /hive/users/ifiddes/longranger-1.2.0/separated_phased_bams/NA12878/*.fq; do
n=`basename $f`
fermi.kit/fermi2.pl unitig -s100k -l100 -p $n $f > $n.mak
done

# ran all makefiles in fermkit folder, moved them back.
find /hive/users/ifiddes/longranger-1.2.0/separated_phased_bams/fermi-assembly_NA12878 | grep mag.gz$ | xargs -n 1 -P 20 -I{} sh -c "fermi.kit/run-calling -t4 /hive/users/ifiddes/chr1_fastas/hg38_chr1_repeatmasked.fasta '{}' | sh"

cd /hive/users/ifiddes/longranger-1.2.0/separated_phased_bams/fermi-assembly_NA12878
find . | grep bam$ | xargs -n 1 -P 10 samtools index
# merge same haplotypes
find . -name "*_1*bam" | xargs samtools merge -h 40149366446_1.bam.sorted.fq.srt.bam h1.bam 
find . -name "*_2*bam" | xargs samtools merge -h 40149366446_2.bam.sorted.fq.srt.bam h2.bam 
samtools index h1.bam
samtools index h2.bam

# i need to redo this, making use of alignment information from LARIAT to create sub-references dynamically

import vcf
import subprocess
from collections import defaultdict
import pysam
from pyfasta import Fasta
import tempfile
import shutil

loci = [x.split()[:3] for x in open("/hive/users/ifiddes/notch2nl_suns/Notch2NL.bed")]
loci = [[x[0], int(x[1]), int(x[2])] for x in loci[:-1]]
loci = sorted(loci, key=lambda x: x[1])


g = "NA12878"
r = pysam.Samfile("/hive/users/ifiddes/longranger-1.2.0/{}/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam".format(g))
ps_map = defaultdict(list)
for x, y in zip(*[loci, ["Notch2", "Notch2NL-D", "Notch2NL-A", "Notch2NL-B", "Notch2NL-C"]]):
    x[1] = x[-1] - 100000
    x[2] = x[2] + 100000
    for read in r.fetch(*x):
        tags = dict(read.tags)
        if "PS" not in tags:
            continue
        ps_map[(tags["PS"], tags["HP"])].append(read)


sorted_reads = []
for x, y in ps_map.iteritems():
    sorted_reads.append(["_".join(map(str, x)), sorted(y, key=lambda x: x.pos)])


padding = 5000

def parse_map(name, reads):
    f = Fasta("/hive/groups/recon/projs/gorilla_eichler/pipeline_data/assemblies/susie_3_2/human.fa")
    start, stop = reads[0].pos - padding, reads[-1].pos + padding
    seq = f['chr1'][start:stop]
    tmp_ref = tempfile.mkstemp()[1]
    with open(tmp_ref, "w") as outf:
        outf.write(">chr1_{}-{}\n{}\n".format(start, stop, seq))
    tmp_fq = tempfile.mkstemp()[1]
    with open(tmp_fq, "w") as outf:
        for i, read in enumerate(reads):
            outf.write("@{}\n{}\n+\n{}\n".format(read.qname + str(i), read.seq, read.qual))
    subprocess.call("/cluster/home/ifiddes/fermikit/fermi.kit/fermi2.pl unitig -s {}k -l100 -p {} {} > {}".format((stop - start) / 1000, tmp_fq, tmp_fq, tmp_fq + ".mak"), shell=True)
    subprocess.call("make -f {}".format(tmp_fq + ".mak"), shell=True)
    subprocess.call("bwa index {}".format(tmp_ref), shell=True)
    subprocess.call("/cluster/home/ifiddes/fermikit/fermi.kit/run-calling {} {} | sh".format(tmp_ref, tmp_fq + ".mag.gz"), shell=True)
    header = {"HD": {"VN": "1.3"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}
    with pysam.Samfile(os.path.join("/hive/users/ifiddes/longranger-1.2.0/separated_phased_bams/fermi-assembly_NA12878", name + ".bam"), "wb", header=header) as outf:
        for read in pysam.Samfile(tmp_fq + ".srt.bam"):
            read.pos = read.pos + start
            outf.write(read)


for name, reads in sorted_reads:
    parse_map(name, reads)


find . -name "*_1*bam" | xargs samtools merge -h 40120000161_1.bam h1.bam 
find . -name "*_2*bam" | xargs samtools merge -h 40120000161_2.bam h2.bam 
samtools index h1.bam
samtools index h2.bam