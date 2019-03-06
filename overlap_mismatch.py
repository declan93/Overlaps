# coding=utf-8
__author__ = "Declan Bennett"
__email__ = "d.bennett1@nuigalway.ie"
__version__ = "0.1"
__status__ = "development"

import os
import re
import sys
import time
import math
from collections import Counter
from collections import defaultdict
from matplotlib.offsetbox import AnchoredText
import pysam
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


strt = time.time()

# Files ###########
BAM = sys.argv[1]  # input BAM,CRAM will need to check + index
SNP = sys.argv[2]  # Known polymorphs maybe revisit this with maf ~1-2%? should exclude any rare variants that may
REF = sys.argv[3]
# be called wrong leaving true pop variants. contig = sys.argv[3] # chromosome start = sys.argv[4] # chunk start end
# = sys.argv[5] # chunk end di = sys.argv[6] if not os.path.exists(sys.argv[6]): os.makedirs(sys.argv[6])
# ##############

nucs = {'A', 'C', 'G', 'T'}  #
snp_loc = defaultdict(int)
sub_dict = {}  # This needs to be emptied after every pair processed
sub_dict_read2 = {}  # This also needs to be emptied
mismatch = {}
mismatch_ref = {}
mismatch_r1_f = {}
mismatch_r1_r = {}
mm_quals = defaultdict(list)
m_quals = defaultdict(list)
overlap_names = defaultdict(list)

# Thresholds #####
bs_q = 20  # base quality thrshold
cvg = 0  # 'coverage
mp_q = 40  # read mapping quality

print(f'Base Quality Threshold is {bs_q}', flush=True)
print(f'Putative Mutation Sequencing Depth Threshold is {cvg}', flush=True)
print(f'Mapping Quality is {mp_q}\n', flush=True)


# reverse complement not currently used will need it for strand asymmetry
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna])  # change dna[::-1] for sequence reversing


def read_r1_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
    if mpos in range(pos, pos + rlen) and '^' not in md1 and '^' not in md2 and str(rlen) != md1:
        if str(mlen) != md2:
            if md1 is not None:
                if md2 is not None:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    else:
        return False


def read_r2_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
    if pos in range(mpos, mpos + mlen) and '^' not in md1 and '^' not in md2 and str(rlen) != md1:
        if md2 != str(mlen):
            if md1 is not None:
                if md2 is not None:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    else:
        return False


cig_g = ["D", "d", "I", "i", "S", "s", "H", "h"]


def check_cigar(cig_pat, cigar_string):
    """
    :param cig_pat: basestring
    :type cigar_string: basestring
    """
    if any(ele in cigar_string for ele in cig_pat):
        return True
    else:
        return False


def ref_matches(align_type, idx, rcig, mcig, mdz):  # remember to MDZ=MDZ1 etc
    if align_type.isdigit() and 0 != int(align_type) and check_cigar(cig_g, rcig) is False and \
            check_cigar(cig_g, mcig) is False: # and len(mdz) - 1 != idx   # sub out matching bases and progress through read
        return True
    else:
        return False


def ref_mismatches(align_type, nucs, rcig, mcig):
    if align_type.isalpha() and align_type in nucs and len(align_type) == 1 and check_cigar(cig_g, rcig) is False and \
            check_cigar(cig_g, mcig) is False:
        return True
    else:
        return False


rcig = None
def read_pair_generator(bam, cigar_check, r_cig, cig, region_string=None):
    # Generate read pairs in a BAM file or within a region string.
    # Reads are added to read_dict until a pair is found.

    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        rcig = read.cigarstring
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate or \
                read.mapq < mp_q or cigar_check(cig_g, rcig) is True:
            continue
        try:
            md_tag = read.get_tag('MD')
        except KeyError:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]



#######

num_reads = 999999 # progression output
read_counts = 0
ovrlp_seq = 0
gtr1 = 0
car1 = 0
gtr2 = 0
car2 = 0
C_sum_r1 = 0
G_sum_r1 = 0
C_sum_r2 = 0
G_sum_r2 = 0
r1f_c = 0
r1r_c = 0
name = BAM.split('.')[0].split('/')[-1]  # for files with absolute paths. files should be named sample_name.bam
index = BAM+".crai"
#print(index)
file_type = BAM.split('.')[-1].split('/')[-1]
if file_type == "cram":
    compres = "rc"
elif file_type == "bam":
    compres = "rb"
elif file_type =="sam":
    compres = "r"

if not os.path.exists(name):
    os.makedirs(name)


print(f'Detected {file_type} file type\nRunning get_overlaps to estimate mismatches for sample {name}\n', flush=True)


with pysam.AlignmentFile(BAM, compres) as SAM:
    if SAM.is_cram:
        print(f'File: {BAM} is Cram formatted', flush=True)
    elif SAM.is_bam:
        print(f'File: {BAM} is Bam formatted', flush=True)
    elif SAM.is_sam:
        print(f'File: {BAM} is Sam formatted', flush=True)

    for read, read2 in read_pair_generator(bam=SAM, cigar_check=check_cigar,r_cig=rcig, cig=cig_g, region_string=None):
        if read is None or read2 is None:
            continue

        pos = read.pos
        rlen = read.infer_query_length()  # read length **** cProfiler showed that 1/8 of the time was spent on
        # read.get_reference_sequence() going to switch to infer query length. This uses the Cigar and skips
        # clipped bases which is good because i am not considering complex alignments yet
        mlen = read2.infer_query_length()  # read2 length
        rref = read.reference_name  #
        mpos = read2.pos
        md1 = read.get_tag('MD')  # MD:z tag
        md2 = read2.get_tag('MD')
        rcig = read.cigarstring
        mcig = read2.cigarstring
        MDZ1 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md1)
        # need to make sure that this works for indels that are not in cigar. ie move len(align_type)-1 in casee^AC
        # indel THIS IS IMPORTANT NEEDS TO BE 100% CORRECT PLUS ROBUST TO EDGECASES. Pretty sure the middle clause is
        # not required ? means before is optional so we know MD only returns integers. but it is currently working
        # so come back later.
        MDZ2 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md2)

        if read_r1_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
            ovl = (int(pos) + int(rlen) - int(mpos))
#            if ovl >= 100:
 #               print(ovl)
            ovrlp_seq += ovl
            for i, j in enumerate(read.query_alignment_sequence):
                if j == 'C' and read.query_qualities[i] >= bs_q:
                    C_sum_r1 += 1
                elif j == 'G' and read.query_qualities[i] >= bs_q:
                    G_sum_r1 += 1

            for i, j in enumerate(read2.query_alignment_sequence):
                if j == 'C' and read2.query_qualities[i] >= bs_q:
                    G_sum_r2 += 1
                elif j == 'G' and read2.query_qualities[i] >= bs_q:
                    C_sum_r2 += 1
            mm_cnt = 0
            read_pos = 0  # used to catch substitution position
            read2_pos = 0
            m1_pos = int(pos)
            m2_pos = int(mpos)
            ##todo hash overlapping reads with positions and qualities this way we can tell if a read has
            # multiple overlapping posiions also take read_pos as I want to see where in the reads are we seeing the mismatches.
            for idx, align_type in enumerate(MDZ1):
                rmm_q = 0
                if ref_matches(align_type, idx, rcig, mcig, MDZ1):
                    read_pos = read_pos + int(align_type)
                    m1_pos = int(m1_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig):  # and read.query_qualities[(int(read_pos))] >= bs_q:
                    # read base call quality coverage funct will check for read2
                    ref = align_type
                    mm_cnt = mm_cnt + 1  # mutation per read count
                    m1_pos = m1_pos + 1
                    mm_b = str(read.query_alignment_sequence)[read_pos]  # mismatching base
                    rmm_q = read.query_qualities[read_pos]
                    m1_key = str(rref) + ':' + str(m1_pos)
                    if int(rmm_q) >= bs_q:
                        sub_dict[m1_key] = ref + mm_b
                    read_pos += 1
                    if read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                        car1 += 1
                    elif not read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                        gtr1 += 1

                    if read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                        gtr1 += 1
                    elif not read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                        car1 += 1
            for idx2, align_type in enumerate(MDZ2):
                rmm_q = 0
                if ref_matches(align_type, idx2, rcig, mcig, MDZ2):
                    read2_pos = read2_pos + int(align_type)
                    m2_pos = int(m2_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig):
                    ref = align_type
                    mm_cnt = mm_cnt + 1
                    m2_pos = m2_pos + 1
                    rmm_b = str(read2.query_alignment_sequence)[read2_pos]  # mismatching base
                    rmm_q = read2.query_qualities[read2_pos]
                    if read2.is_reverse and (ref + rmm_b) == 'GT' and rmm_q >= bs_q:
                        car2 += 1
                    elif not read2.is_reverse and (ref + rmm_b) == 'GT' and rmm_q >= bs_q:
                        gtr2 += 1

                    if read2.is_reverse and (ref + rmm_b) == 'CA' and rmm_q >= bs_q:
                        gtr2 += 1
                    elif not read2.is_reverse and (ref + rmm_b) == 'CA' and rmm_q >= bs_q:
                        car2 += 1
                    m2_key = str(rref) + ':' + str(m2_pos)  # + ':' + read2.query_name
                    if int(rmm_q) >= bs_q:
                        sub_dict_read2[m2_key] = ref + rmm_b
                    if m2_key in sub_dict and m2_key in sub_dict_read2 and sub_dict[m2_key] == \
                            sub_dict_read2[m2_key]:  # and not m2_key in snp_loc:
                        if read.is_reverse:
                            mismatch[m2_key] = reverse_complement(ref + rmm_b)
                            r1r_c += 1
                            mismatch_ref[m2_key] = ref + rmm_b
                            mismatch_r1_r[m2_key] = ref + rmm_b
                        else:
                            mismatch[m2_key] = ref + rmm_b
                            mismatch_ref[m2_key] = ref + rmm_b
                            mismatch_r1_f[m2_key] = ref + rmm_b
                            r1f_c += 1
                        if m2_key in overlap_names.keys():
                            overlap_names[m2_key].append(str(ref+rmm_b))
                            overlap_names[m2_key].append(read.query_name)
                            overlap_names[m2_key].append(str(ovl))
                        else:
                            overlap_names[m2_key] = [(ref+rmm_b), read.query_name, str(ovl)]
                    read2_pos += 1
            sub_dict = {}
            sub_dict_read2 = {}

        elif read_r2_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
            ovl = (int(mpos) + int(mlen) - int(pos))
           # if ovl >= 100:
            #    print(ovl)
            ovrlp_seq += ovl
            for i, j in enumerate(read.query_alignment_sequence):
                if j == 'C' and read.query_qualities[i] >= bs_q:
                    C_sum_r1 += 1
                elif j == 'G' and read.query_qualities[i] >= bs_q:
                    G_sum_r1 += 1

            for i, j in enumerate(read2.query_alignment_sequence):
                if j == 'C' and read2.query_qualities[i] >= bs_q:
                    G_sum_r2 += 1
                elif j == 'G' and read2.query_qualities[i] >= bs_q:
                    C_sum_r2 += 1
            mm_cnt = 0
            read_pos = 0  # used to catch substitution position
            read2_pos = 0
            m1_pos = int(pos)
            m2_pos = int(mpos)

            for idx3, align_type in enumerate(MDZ1):
                rmm_q = 0
                if ref_matches(align_type, idx3, rcig, mcig, MDZ1):
                    read_pos = read_pos + int(align_type)
                    m1_pos = int(m1_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig): # and read.query_qualities[(int(read_pos))] >= bs_q:
                    ref = align_type
                    mm_cnt = mm_cnt + 1  # mutation per read count
                    m1_pos = m1_pos + 1
                    mm_b = str(read.query_alignment_sequence)[read_pos]  # mismatching base
                    if read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                        car1 += 1
                    elif not read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                        gtr1 += 1

                    if read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                        gtr1 += 1
                    elif not read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                        car1 += 1
                    rmm_q = read.query_qualities[read_pos]
                    m1_key = str(rref) + ':' + str(m1_pos)  # + ':' + read.query_name
                    if int(rmm_q) >= bs_q:
                        sub_dict[m1_key] = ref + mm_b
                    read_pos += 1

            for idx4, align_type in enumerate(MDZ2):
                rmm_q = 0
                if ref_matches(align_type, idx4, rcig, mcig, MDZ2):
                    read2_pos = read2_pos + int(align_type)
                    m2_pos = int(m2_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig):
                    ref = align_type
                    mm_cnt = mm_cnt + 1
                    m2_pos = m2_pos + 1
                    rmm_b = str(read2.query_alignment_sequence)[read2_pos]  # mismatching base
                    rmm_q = read2.query_qualities[read2_pos]
                    if read2.is_reverse and (ref + rmm_b) == 'GT' and rmm_q >= bs_q:
                        car2 += 1
                    elif not read2.is_reverse and (ref + rmm_b) == 'GT' and rmm_q >= bs_q:
                        gtr2 += 1

                    if read2.is_reverse and (ref + rmm_b) == 'CA' and rmm_q >= bs_q:
                        gtr2 += 1
                    elif not read2.is_reverse and (ref + rmm_b) == 'CA' and rmm_q >= bs_q:
                        car2 += 1
                    m2_key = str(rref) + ':' + str(m2_pos)  # + ':' + read2.query_name
                    if int(rmm_q) >= bs_q:
                        sub_dict_read2[m2_key] = ref + rmm_b
                    if m2_key in sub_dict and m2_key in sub_dict_read2 and sub_dict[m2_key] == \
                            sub_dict_read2[m2_key]:  # and not m2_key in snp_loc:
                        if read.is_reverse:
                            mismatch[m2_key] = reverse_complement(ref + rmm_b)
                            mismatch_ref[m2_key] = ref + rmm_b
                            mismatch_r1_r[m2_key] = ref + rmm_b
                            r1r_c += 1
                        else:
                            mismatch[m2_key] = ref + rmm_b
                            r1f_c += 1
                            mismatch_ref[m2_key] = ref + rmm_b
                            mismatch_r1_f[m2_key] = ref + rmm_b
                        if m2_key in overlap_names.keys():
                            overlap_names[m2_key].append(str(ref+rmm_b))
                            overlap_names[m2_key].append(read.query_name)
                            overlap_names[m2_key].append(str(ovl))
                        else:
                            overlap_names[m2_key] = [(ref+rmm_b), read.query_name, str(ovl)]
                    read2_pos += 1

            sub_dict = {}
            sub_dict_read2 = {}
        else:
            continue

        read_counts = read_counts + 1

        if int(read_counts) > int(num_reads):
            print(f'Number of paired reads processed that satisfy thresholds is {read_counts} ', end="\r", flush=True)
            num_reads += 1000000

    init_mism = len(mismatch_ref)
    with open(SNP, "r") as snps:
        for line in snps:
            if line.strip() in mismatch_ref:
                del mismatch_ref[line.strip()]
                del overlap_names[line.strip()]

    print(gtr1, gtr2, car1, car2, G_sum_r1, G_sum_r2, C_sum_r1, C_sum_r2)
    G_iv = str((((gtr1 + car2) / (G_sum_r1 + C_sum_r2)) / ((car1 + gtr2) / (C_sum_r1 + G_sum_r2))))
    G_ivl = math.log2((((gtr1 + car2) / (G_sum_r1 + C_sum_r2)) / ((car1 + gtr2) / (C_sum_r1 + G_sum_r2))))
    print(f'\nThe G -> T imbalance is {G_iv}')
    print(f'\nThe G -> T log2 imbalance is {G_ivl}\n\nThe initial number of mismatches was {init_mism} while the number of found snps was {init_mism - len(mismatch_ref)}')

    ##todo Tidy this into functions
    for key, value in mismatch_ref.items():
        fields = key.strip().split(':')
        pos_coverage = {}
        for cov in SAM.pileup(fields[0], int(fields[1])-1, int(fields[1]), min_base_quality=0, nofilter=True):
            if int(cov.pos) == int(fields[1]):
                for pileupread in cov.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.mapq >= mp_q:
                        if int(fields[1]) == int(pileupread.alignment.pos + pileupread.query_position):
                            if pileupread.alignment.query_sequence[pileupread.query_position - 1] == value[1]:
                                pos_coverage[key] = cov.n
                                if key in mm_quals.keys():
                                    mm_quals[key].append(
                                        pileupread.alignment.query_qualities[pileupread.query_position - 1])
                                else:
                                    mm_quals[key] = [cov.n, pileupread.alignment.query_qualities[
                                        pileupread.query_position - 1]]

                            elif pileupread.alignment.query_sequence[pileupread.query_position - 1] == value[0]:
                                if key in m_quals.keys():
                                    m_quals[key].append(
                                        pileupread.alignment.query_qualities[pileupread.query_position - 1])
                                else:
                                    m_quals[key] = [cov.n, pileupread.alignment.query_qualities[pileupread.query_position - 1]]
                            else:
                                continue
                        else:
                            continue

print('Writing mutations to file \n', flush=True)

with open(name + '/' + name + '_substitutions.out', 'w') as out:
    out.write('pos\tref_alt\tmismatching_qualities\tmatching_qualities\n')
    for key, value in mismatch_ref.items():
        outstring = str(key) + '\t' + str(value) + '\t' + ",".join([str(x) for x in mm_quals[key]]) + '\t' + ",".join(
            [str(y) for y in m_quals[key]]) + '\n'
        out.write(outstring)

print('Calculating Mutation counts where Pair agree \n', flush=True)

with open(name + '/' + name + '_sub_counts.out', 'w') as cnts:
    subs = Counter(mismatch_ref.values())
    cnts.write('Sub Count Overlap\n')
    for key, value in subs.items():
        out = str(key) + ' ' + str(value) + ' ' + str(ovrlp_seq) + '\n'
        cnts.write(out)

print(f'making VCF for {name}', flush=True)
print(overlap_names)
with open(name + '/' + name + '_substitutions.out', 'r') as IN:
    next(IN)
    with open(name + '/' + name + '.vcf', 'w') as out:
        out.write(
            '##fileformat=VCFv4.0\n##contig=<ID=1,assembly=b37,length=249250621>\n##contig=<ID=2,assembly=b37,'
            'length=243199373>\n##contig=<ID=3,assembly=b37,length=198022430>\n##contig=<ID=4,assembly=b37,'
            'length=191154276>\n##contig=<ID=5,assembly=b37,length=180915260>\n##contig=<ID=6,assembly=b37,'
            'length=171115067>\n##contig=<ID=7,assembly=b37,length=159138663>\n##contig=<ID=8,assembly=b37,'
            'length=146364022>\n##contig=<ID=9,assembly=b37,length=141213431>\n##contig=<ID=10,assembly=b37,'
            'length=135534747>\n##contig=<ID=11,assembly=b37,length=135006516>\n##contig=<ID=12,assembly=b37,'
            'length=133851895>\n##contig=<ID=13,assembly=b37,length=115169878>\n##contig=<ID=14,assembly=b37,'
            'length=107349540>\n##contig=<ID=15,assembly=b37,length=102531392>\n##contig=<ID=16,assembly=b37,'
            'length=90354753>\n##contig=<ID=17,assembly=b37,length=81195210>\n##contig=<ID=18,assembly=b37,'
            'length=78077248>\n##contig=<ID=19,assembly=b37,length=59128983>\n##contig=<ID=20,assembly=b37,'
            'length=63025520>\n##contig=<ID=21,assembly=b37,length=48129895>\n##contig=<ID=22,assembly=b37,'
            'length=51304566>\n##contig=<ID=GL000191.1,assembly=b37,length=106433>\n##contig=<ID=GL000192.1,'
            'assembly=b37,length=547496>\n##contig=<ID=GL000193.1,assembly=b37,'
            'length=189789>\n##contig=<ID=GL000194.1,assembly=b37,length=191469>\n##contig=<ID=GL000195.1,'
            'assembly=b37,length=182896>\n##contig=<ID=GL000196.1,assembly=b37,length=38914>\n##contig=<ID=GL000197.1,'
            'assembly=b37,length=37175>\n##contig=<ID=GL000198.1,assembly=b37,length=90085>\n##contig=<ID=GL000199.1,'
            'assembly=b37,length=169874>\n##contig=<ID=GL000200.1,assembly=b37,'
            'length=187035>\n##contig=<ID=GL000201.1,assembly=b37,length=36148>\n##contig=<ID=GL000202.1,assembly=b37,'
            'length=40103>\n##contig=<ID=GL000203.1,assembly=b37,length=37498>\n##contig=<ID=GL000204.1,assembly=b37,'
            'length=81310>\n##contig=<ID=GL000205.1,assembly=b37,length=174588>\n##contig=<ID=GL000206.1,assembly=b37,'
            'length=41001>\n##contig=<ID=GL000207.1,assembly=b37,length=4262>\n##contig=<ID=GL000208.1,assembly=b37,'
            'length=92689>\n##contig=<ID=GL000209.1,assembly=b37,length=159169>\n##contig=<ID=GL000210.1,assembly=b37,'
            'length=27682>\n##contig=<ID=GL000211.1,assembly=b37,length=166566>\n##contig=<ID=GL000212.1,assembly=b37,'
            'length=186858>\n##contig=<ID=GL000213.1,assembly=b37,length=164239>\n##contig=<ID=GL000214.1,'
            'assembly=b37,length=137718>\n##contig=<ID=GL000215.1,assembly=b37,'
            'length=172545>\n##contig=<ID=GL000216.1,assembly=b37,length=172294>\n##contig=<ID=GL000217.1,'
            'assembly=b37,length=172149>\n##contig=<ID=GL000218.1,assembly=b37,'
            'length=161147>\n##contig=<ID=GL000219.1,assembly=b37,length=179198>\n##contig=<ID=GL000220.1,'
            'assembly=b37,length=161802>\n##contig=<ID=GL000221.1,assembly=b37,'
            'length=155397>\n##contig=<ID=GL000222.1,assembly=b37,length=186861>\n##contig=<ID=GL000223.1,'
            'assembly=b37,length=180455>\n##contig=<ID=GL000224.1,assembly=b37,'
            'length=179693>\n##contig=<ID=GL000225.1,assembly=b37,length=211173>\n##contig=<ID=GL000226.1,'
            'assembly=b37,length=15008>\n##contig=<ID=GL000227.1,assembly=b37,length=128374>\n##contig=<ID=GL000228.1,'
            'assembly=b37,length=129120>\n##contig=<ID=GL000229.1,assembly=b37,length=19913>\n##contig=<ID=GL000230.1,'
            'assembly=b37,length=43691>\n##contig=<ID=GL000231.1,assembly=b37,length=27386>\n##contig=<ID=GL000232.1,'
            'assembly=b37,length=40652>\n##contig=<ID=GL000233.1,assembly=b37,length=45941>\n##contig=<ID=GL000234.1,'
            'assembly=b37,length=40531>\n##contig=<ID=GL000235.1,assembly=b37,length=34474>\n##contig=<ID=GL000236.1,'
            'assembly=b37,length=41934>\n##contig=<ID=GL000237.1,assembly=b37,length=45867>\n##contig=<ID=GL000238.1,'
            'assembly=b37,length=39939>\n##contig=<ID=GL000239.1,assembly=b37,length=33824>\n##contig=<ID=GL000240.1,'
            'assembly=b37,length=41933>\n##contig=<ID=GL000241.1,assembly=b37,length=42152>\n##contig=<ID=GL000242.1,'
            'assembly=b37,length=43523>\n##contig=<ID=GL000243.1,assembly=b37,length=43341>\n##contig=<ID=GL000244.1,'
            'assembly=b37,length=39929>\n##contig=<ID=GL000245.1,assembly=b37,length=36651>\n##contig=<ID=GL000246.1,'
            'assembly=b37,length=38154>\n##contig=<ID=GL000247.1,assembly=b37,length=36422>\n##contig=<ID=GL000248.1,'
            'assembly=b37,length=39786>\n##contig=<ID=GL000249.1,assembly=b37,length=38502>\n##contig=<ID=MT,'
            'assembly=b37,length=16569>\n##contig=<ID=NC_007605,assembly=b37,length=171823>\n##contig=<ID=X,'
            'assembly=b37,length=155270560>\n##contig=<ID=Y,assembly=b37,length=59373566>\n##contig=<ID=hs37d5,'
            'assembly=b37,length=35477943>\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
            '##INFO=<ID=AC,Number=1,Type=Integer,Description="Nonrefernece Count">\n'
            '##INFO=<ID=OV,Number=1,Type=Integer,Description="Overlapping mismatch Count">\n'
            '##INFO=<ID=OVL_LEN,Number=1,Type=String,Description="Overlapping sequencelength">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype tag">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(name))
        # move to new file hash headers and only write headers that have mutations
        for line in IN:
            fields = line.strip().split()
            ref = list(fields[1])[0]
            mut = list(fields[1])[1]
            location = fields[0].strip().split(':')
            try:
                key = mm_quals[fields[0]][0]
            except IndexError:
                try:
                    key = m_quals[fields[0]][0]
                except IndexError:
                    key = "NA"
            out_str = location[0] + '\t' + location[1] + "\t" + fields[
                0] + '\t' + ref + '\t' + mut + f'\t.\t.\tDP={key};AC={len(mm_quals[fields[0]])-1};' \
                f'OV={(len(overlap_names[fields[0]][0::3]))};' \
                f'OVL_LEN={",".join(overlap_names[fields[0]][2::3])}\tGT\t1/0\n'
            out.write(out_str)

tme = (time.time() - strt) / 3600

## plot

cvg_mism = {}
for key, value in mismatch.items():
    if mm_quals[key]:
        if mm_quals[key][0] >= cvg:
            cvg_mism[key] = value
    elif m_quals[key]:
        if m_quals[key][0] >= cvg:
            cvg_mism[key] = value
    else:
        continue

cvg_mism_ref = {}
for key, value in mismatch_ref.items():
    if mm_quals[key]:
        if mm_quals[key][0] >= cvg:
            cvg_mism_ref[key] = value
    elif m_quals[key]:
        if m_quals[key][0] >= cvg:
            cvg_mism_ref[key] = value
    else:
        continue
        
cvg_mism_r1f = {}
for key, value in mismatch_r1_f.items():
    if mm_quals[key]:
        if mm_quals[key][0] >= cvg:
            cvg_mism_r1f[key] = value
    elif m_quals[key]:
        if m_quals[key][0] >= cvg:
            cvg_mism_r1f[key] = value
    else:
        continue

cvg_mism_r1r = {}
for key, value in mismatch_r1_r.items():
    if mm_quals[key]:
        if mm_quals[key][0] >= cvg:
            cvg_mism_r1r[key] = value
    elif m_quals[key]:
        if m_quals[key][0] >= cvg:
            cvg_mism_r1r[key] = value
    else:
        continue

subs1 = Counter(cvg_mism.values())
subs2 = Counter(cvg_mism_ref.values())
subs3 = Counter(cvg_mism_r1f.values())
subs4 = Counter(cvg_mism_r1r.values())

print(f'\n % of Read 1 mapping to forward strand is {(r1f_c/(r1f_c+r1r_c))}\n')

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(13,13))

### Read 1 aligned

per_CTR1 = (subs1["CA"]/sum(subs1.values()), subs1["GT"]/sum(subs1.values()),
          subs1["CG"]/sum(subs1.values()), subs1["GC"]/sum(subs1.values()),
          subs1["CT"]/sum(subs1.values()), subs1["GA"]/sum(subs1.values()),
          subs1["TA"]/sum(subs1.values()), subs1["AT"]/sum(subs1.values()),
          subs1["TC"]/sum(subs1.values()), subs1["AG"]/sum(subs1.values()),
          subs1["TG"]/sum(subs1.values()), subs1["AC"]/sum(subs1.values()))

# create plot
bar_width = 1
opacity = 0.8

y_pos = [0,1,4,5,8,9,12,13,16,17,20,21]
ax1.bar(y_pos, per_CTR1, bar_width,
                 alpha=opacity,
                 color=['b','g','b','g','b','g','b','g','b','g','b','g'])

ax1.set_xlabel('Substitutions')
ax1.set_ylabel('% Relative contribution')
ax1.set_title(f'Mutational frequencies where read 1 have no bearing\n to the reference for {name}')
ax1.set_xticks(y_pos)
ax1.set_xticklabels(('CA', 'GT', 'CG', 'GC', 'CT', 'GA', 'TA', 'AT', 'TC', 'AG', 'TG', 'AC'))

## Reference orientated

per_CTR1F = (subs2["CA"]/sum(subs2.values()), subs2["GT"]/sum(subs2.values()),
          subs2["CG"]/sum(subs2.values()), subs2["GC"]/sum(subs2.values()),
          subs2["CT"]/sum(subs2.values()), subs2["GA"]/sum(subs2.values()),
          subs2["TA"]/sum(subs2.values()), subs2["AT"]/sum(subs2.values()),
          subs2["TC"]/sum(subs2.values()), subs2["AG"]/sum(subs2.values()),
          subs2["TG"]/sum(subs2.values()), subs2["AC"]/sum(subs2.values()))

# create plot

ax2.bar(y_pos, per_CTR1F, bar_width,
                 alpha=opacity,
                 color=['b','g','b','g','b','g','b','g','b','g','b','g'])

ax2.set_xlabel('Substitutions')
ax2.set_ylabel('% Relative contribution')
ax2.set_title(f'Mutational frequencies with respect to\n the reference for {name}')
ax2.set_xticks(y_pos)
ax2.set_xticklabels(('CA', 'GT', 'CG', 'GC', 'CT', 'GA', 'TA', 'AT', 'TC', 'AG', 'TG', 'AC'))

### Read1 forward strand

R1F = (subs3["CA"]/sum(subs3.values()), subs3["GT"]/sum(subs3.values()),
          subs3["CG"]/sum(subs3.values()), subs3["GC"]/sum(subs3.values()),
          subs3["CT"]/sum(subs3.values()), subs3["GA"]/sum(subs3.values()),
          subs3["TA"]/sum(subs3.values()), subs3["AT"]/sum(subs3.values()),
          subs3["TC"]/sum(subs3.values()), subs3["AG"]/sum(subs3.values()),
          subs3["TG"]/sum(subs3.values()), subs3["AC"]/sum(subs3.values()))

# create plot

ax3.bar(y_pos, R1F, bar_width, alpha=opacity,color=['b','g','b','g','b','g','b','g','b','g','b','g'])

ax3.set_xlabel('Substitutions')
ax3.set_ylabel('% Relative contribution')
ax3.set_title(f'Mutational frequencies where read 1\n maps to the + for {name}')
ax3.set_xticks(y_pos)
ax3.set_xticklabels(('CA', 'GT', 'CG', 'GC', 'CT', 'GA', 'TA', 'AT', 'TC', 'AG', 'TG', 'AC'))

### Read1 maps to the reverse
R1R = (subs4["CA"]/sum(subs4.values()), subs4["GT"]/sum(subs4.values()),
          subs4["CG"]/sum(subs4.values()), subs4["GC"]/sum(subs4.values()),
          subs4["CT"]/sum(subs4.values()), subs4["GA"]/sum(subs4.values()),
          subs4["TA"]/sum(subs4.values()), subs4["AT"]/sum(subs4.values()),
          subs4["TC"]/sum(subs4.values()), subs4["AG"]/sum(subs4.values()),
          subs4["TG"]/sum(subs4.values()), subs4["AC"]/sum(subs4.values()))


box_text = f'G->T log2 damage is {G_ivl:.3f}\nThe read1 + strand % = {(r1f_c/(r1f_c+r1r_c)):.3f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3.text(0.5, 0.95, box_text, transform=ax1.transAxes, fontsize=14)
# create plot

ax4.bar(y_pos, R1R, bar_width,
                 alpha=opacity,
                 color=['b','g','b','g','b','g','b','g','b','g','b','g'])

ax4.set_xlabel('Substitutions')
ax4.set_ylabel('% Relative contribution')
ax4.set_title(f'Mutational frequencies where read 1\n maps to the - strand for {name}')
ax4.set_xticks(y_pos)
ax4.set_xticklabels(('CA', 'GT', 'CG', 'GC', 'CT', 'GA', 'TA', 'AT', 'TC', 'AG', 'TG', 'AC'))
#ax4.savefig(name + "/" + name + '_read1_only_reverse_strand.png')

f.tight_layout()
f.savefig(name + '/' + name + '_.png')
print(f'overlap_mismatch has completed. The number of overlapping bases is {ovrlp_seq}\n'
      f'The number of putative mismatches is {sum(subs.values())}\nThe time taking to analyse {name} was {tme} hrs', flush=True)
