import os
import re
import sys
import time
import math
from collections import Counter
from collections import defaultdict
import pysam

strt = time.time()

# Files ###########
BAM = sys.argv[1]  # input BAM,CRAM will need to check + index
REF = sys.argv[2]
#SNP = sys.argv[2]  # Known polymorphs maybe revisit this with maf ~1-2%? should exclude any rare variants that may
# be called wrong leaving true pop variants. contig = sys.argv[3] # chromosome start = sys.argv[4] # chunk start end
# = sys.argv[5] # chunk end di = sys.argv[6] if not os.path.exists(sys.argv[6]): os.makedirs(sys.argv[6])
# ##############

nucs = {'A', 'C', 'G', 'T'}  #
snp_loc = defaultdict(int)
sub_dict = {}  # This needs to be emptied after every pair processed
sub_dict_read2 = {}  # This also needs to be emptied
mismatch = {}
mm_quals = defaultdict(list)
m_quals = defaultdict(list)
overlap_names = defaultdict(list)

# Thresholds #####
bs_q = 20  # base quality thrshold
cvg = 1  # 'coverage
mp_q = 20  # read mapping quality

print(f'Base Quality Threshold is {bs_q}', flush=True)
print(f'Putative Mutation Sequencing Depth Threshold is {cvg}', flush=True)
print(f'Mapping Quality is {mp_q}\n', flush=True)


# reverse complement not currently used will need it for strand asymmetry
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


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


def ref_matches(align_type, idx, rcig, mcig, mdz): #, check_cigar):  # remember to MDZ=MDZ1 etc
    if align_type.isdigit() and 0 != int(align_type) and len(mdz) - 1 != idx: # and check_cigar(cig_g, rcig) is False and \
        #    check_cigar(cig_g, mcig) is False:  # sub out matching bases and progress through read
        return True
    else:
        return False



def ref_mismatches(align_type, nucs, rcig, mcig):
    if align_type.isalpha() and align_type in nucs and len(align_type) == 1: # and check_cigar(cig_g, rcig) is False and \
      #      check_cigar(cig_g, mcig) is False:
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

name = BAM.split('.')[0].split('/')[-1]  # for files with absolute paths. files should be named sample_name.bam
file_type = BAM.split('.')[-1].split('/')[-1]
if file_type == "cram":
    compres = "rc"
elif file_type == "bam":
    compres = "rb"
elif file_type =="sam":
    compres = "r"

#if not os.path.exists(name):
 #   os.makedirs(name)


print(f'Running get_overlaps to estimate mismatches for sample {name}\nDetected {file_type} file\n', flush=True)

r1f = 0
r2f = 0
r1r = 0
r2r = 0

with pysam.AlignmentFile(BAM, compres, reference_filename=REF) as SAM:
    print(f'File {BAM} opened', flush=True)
    for read, read2 in read_pair_generator(bam=SAM, cigar_check=check_cigar,r_cig=rcig, cig=cig_g, region_string=None):
        if read is None or read2 is None:
            continue
        if not read.is_reverse and read2.is_reverse:
            r1f += 1
            r2r += 1
        elif read.is_reverse and not read2.is_reverse:
            r1r += 1
            r2f += 1

        pos = read.pos
        rlen = read.infer_query_length()
        mlen = read2.infer_query_length()  # read2 length
        rref = read.reference_name  #
        mpos = read2.pos
        md1 = read.get_tag('MD')  # MD:z tag
        md2 = read2.get_tag('MD')
        rcig = read.cigarstring
        mcig = read2.cigarstring
        MDZ1 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md1)
        MDZ2 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md2)
        if not read.is_read1:
            print(read)
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
                #print(read_pos, len(str(read.query_alignment_sequence)))
                mm_b = str(read.query_alignment_sequence)[read_pos]  # mismatching base
                rmm_q = read.query_qualities[read_pos]
                m1_key = str(rref) + ':' + str(m1_pos)
                if int(rmm_q) >= bs_q:
                    sub_dict[m1_key] = ref + mm_b
                if read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                    car1 += 1
                elif not read.is_reverse and (ref + mm_b) == 'GT' and rmm_q >= bs_q:
                    gtr1 += 1

                if read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                    gtr1 += 1
                elif not read.is_reverse and (ref + mm_b) == 'CA' and rmm_q >= bs_q:
                    car1 += 1
                read_pos += 1
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
                read2_pos += 1
    print(gtr1, gtr2, car1, car2, G_sum_r1, G_sum_r2, C_sum_r1, C_sum_r2)
    G_iv = str((((gtr1 + car2) / (G_sum_r1 + C_sum_r2)) / ((car1 + gtr2) / (C_sum_r1 + G_sum_r2))))
    G_ivl = math.log2((((gtr1 + car2) / (G_sum_r1 + C_sum_r2)) / ((car1 + gtr2) / (C_sum_r1 + G_sum_r2))))
    print(f'\nThe G -> T imbalance is {G_iv}\nlog is {G_ivl}')
    tme = (time.time() - strt) / 3600
    print(f'Time to run: {tme}')
    print(f'read 1 maps to the forward strand {r1f} times and to the reverse strand {r1r}\nread 2 maps to the forward strand {r2f} times and to the reverse strand {r2r}')