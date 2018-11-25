import os
import pickle
import pysam
import re
import sys
import time
from collections import Counter
from collections import defaultdict

strt = time.time()

# Code to calculate mismatches to the reference genome in overlapping read-pairs. This Version differs from the
# initial implementation that used random access to get the mate pair. I am now hashing the reads until its mate is
# found. This is substantially faster and allows for the analysis of biobank scale datasets


# Also need to write an implementation of Lawrence Etwillers DNA damage extimator. this basically just looks for
# reverse complement agreement in the read2 read. # Need to think about this more i,e does the position of the
# substitution matter in the read or is it just we see more G->T in r1 and less C-A in R2 irrespective if they occur
# at the same loci. I'n not convinced DNA damage will be that critical we are requiring the pair to agree on the
# exact location of the mismatch # Also is there Bias in this as there is strand specific (Transcribed) asymmetry in
# Exome data. Will test with WGS if there is data available with enough sequencing depth . have downloaded the "high
# coverage" HG0096 WGS and currently testing to see any pronounced differences

#  genome = pysam.Fastafile('/home/declan/phd/Homo_sapiens.GRCh37.dna.toplevel.fa.gz)
# this can be used for getting +1-1 bases from reference fasta if choosing to not take from alignment sequence

# # Will also need to derive a qc plan. WGS data for G1K doesnt have MD-Ztags Try to access samtools calmd function
# which can compute md tag. Heng-li has talked about a new field (cs) which should be generally easy to parse. "MD
# done right but too late"

# #Current code now works to calculate exact GRCh37 co-ordinates Mon 29th Oct @23:49 #Have a filtering process set up
# which removes variants with a VAF > 0.4 giving the benefit of doubt that .4 and up is a snp. (at some point store
# these could be interesting)

###################### Files ###########
BAM = sys.argv[1]  # input BAM,CRAM will need to check + index
# SNP = sys.argv[2] # Known polymorphs maybe revisit this with maf ~1-2%? should exclude any rare variants that may
# be called wrong leaving true pop variants. contig = sys.argv[3] # chromosome start = sys.argv[4] # chunk start end
# = sys.argv[5] # chunk end di = sys.argv[6] if not os.path.exists(sys.argv[6]): os.makedirs(sys.argv[6])
# ##############

nucs = {'A', 'C', 'G', 'T'}  #
homo = {'AA', 'CC', 'GG', 'TT'}  # Can be removed was used to flag mistakes in getting mutated nucleotide.
snp_loc = {}
sub_dict = {}  # This needs to be emptied after every pair processed
sub_dict_read2 = {}  # This also needs to be emptied
mismatch = {}
tmismatch = {}

############### Thresholds ##### 
bs_q = 40  # base quality thrshold
cvg = 20  # 'coverage
mp_q = 60  # read mapping quality
################################


###########################
print("Base Quality Threshold is {}".format(bs_q))
print("Putative Mutation Sequencing Depth Threshold is {}".format(cvg))
print("Mapping Quality is {}\n".format(mp_q))


###########################

# with open(SNP,"r") as snps:
# 	print("Hashing Known SNP locations\n")
# 	for pos in snps:
# 		col = pos.strip().split()
# 		key = col[0]+":"+col[1]
# 		snp_loc[key] = 1
# 		print("Hashing Done \n") # this works for the most part. I'm calling every read that covers each
# 		"Mutation" when i calculate coverage. getting the fraction of covering reads that have the mutation
# 		should be trivial to do but cost on timing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# The Idea here is to write the hashed snps as a cPickle object. I have to check if this is feasable interms of
# Snps.pkl size and if any substantial gain is seen in loading dict over reading file and writing to dict ~
# Any small gain will be translated over the full dataset so will be useful 1 second quicker saves 40 mins cpu
# time on G1K data

# print "Loading Known Snps"
# snp_loc = pickle.load( open( "Snps.pkl", "rb" ) )
# print "Snps Loaded/n see"
# print snp_loc
##########################

## reverse complement not currently used will need it for strand asymmetry
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])  # [::-1] reads a string right to left


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


def ref_matches(align_type, idx, rcig, mcig, MDZ):  # remember to MDZ=MDZ1 etc
    if align_type.isdigit() and 0 != int(align_type) and len(MDZ) - 1 != idx and 'D' not in str(rcig) and 'D' not in str(mcig) and 'I' not in str(rcig) and 'I' not in str(mcig) and 'H' not in str(mcig) and 'H' not in str(rcig) and 'S' not in str(mcig) and 'S' not in str(rcig):  # sub out matching bases and progress through read
        return True
    else:
        return False


def filter1(read, SAM, mp_q):
    if read.mapq >= mp_q and read2.mapq >= mp_q and read.tid == read2.tid:
        return True
    else:
        return False


def chk_cov_snp_ql(read, cov, m1_pos, read_pos, bs_q, rref, snp_loc):
    if int(cov.pos) == int(m1_pos) and int(cov.n) >= cvg and read.query_qualities[(int(read_pos))] >= bs_q and (str(rref) + ":" + str(m1_pos)) not in snp_loc:  # and returns coverage for every base in reads also filtering out variants here (snp_loc)
        return True
    else:
        return False


def ref_mismatches(align_type, nucs, rcig, mcig):
    if align_type.isalpha() and align_type in nucs and len(align_type) == 1 and 'D' not in str(rcig) and 'D' not in str(mcig) and 'I' not in str(rcig) and 'I' not in str(mcig) and 'H' not in str(mcig) and 'H' not in str(rcig) and 'S' not in str(mcig) and 'S' not in str(rcig):
        return True
    else:
        return False


def read_pair_generator(bam, region_string=None):

# Generate read pairs in a BAM file or within a region string.
# Reads are added to read_dict until a pair is found.

    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate or read.mapq < mp_q:
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


#######


num_reads = 1499999  # progression output
read_counts = 0

name = BAM.split(".")[0]
if not os.path.exists(name):
    os.makedirs(name)
print("Running get_overlaps to estimate mismatches for sample {}\n".format(name))
ovrlp_seq = 0
with pysam.AlignmentFile(BAM, "rb") as SAM:
    print("File {} opened".format(BAM))
    for read, read2 in read_pair_generator(SAM):
        if read is None or read2 is None:
            continue
        pos = read.pos
        rlen = read.infer_query_length()  # read length **** cProfiler showed that 1/8 of the time was spent on read.get_reference_sequence() going to switch to infer query length. This uses the Cigar and skips clipped bases which is good because i am not considering complex alignments yet
        mlen = read2.infer_query_length()  # read2 length
        rref = read.reference_name  #
        mpos = read2.pos
        md1 = read.get_tag("MD")  # MD:z tag
        md2 = read2.get_tag("MD")
        rcig = read.cigarstring
        mcig = read2.cigarstring
        MDZ1 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md1)  # need to make sure that this works for indels that are not in cigar. ie move len(align_type)-1 in casee^AC indel THIS IS IMPORTANT NEEDS TO BE 100% CORRECT PLUS ROBUST TO EDGECASES. Pretty sure the middle clause is not required ? means before is optional so we know MD only returns integers. but it is currently working so come back later.
        MDZ2 = re.findall(r'[A-Za-z]|-?\d+\.\d+|\d+', md2)
        if read_r1_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
            ovrlp_seq = ovrlp_seq + (int(pos) + int(rlen) - int(mpos))
            mm_cnt = 0
            read_pos = 0  # used to catch substitution position
            read2_pos = 0
            m1_pos = int(pos)
            m2_pos = int(mpos)
            for idx, align_type in enumerate(MDZ1):
                if ref_matches(align_type, idx, rcig, mcig, MDZ1):
                    read_pos = read_pos + int(align_type)
                    m1_pos = int(m1_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig) and read.query_qualities[
                    (int(read_pos))] >= bs_q:  # read base call quality coverage funct will check for read2
                    ref = align_type
                    mm_cnt = mm_cnt + 1  # mutation per read count
                    m1_pos = m1_pos + 1
                    mm_b = str(read.query_alignment_sequence)[read_pos]  # mismatching base
                    m1_key = str(rref) + ":" + str(m1_pos)
                    sub_dict[m1_key] = ref + mm_b
                    read_pos = read_pos + 1
            for idx2, align_type in enumerate(MDZ2):
                if ref_matches(align_type, idx2, rcig, mcig, MDZ2):
                    read2_pos = read2_pos + int(align_type)
                    m2_pos = int(m2_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig):
                    ref = align_type
                    mm_cnt = mm_cnt + 1
                    m2_pos = m2_pos + 1
                    rmm_b = str(read2.query_alignment_sequence)[read2_pos]  # mismatching base
                    for cov in SAM.pileup(str(rref), int(m2_pos), int(m2_pos) + 1):  # this seems to be quite slow as it pulls all the reads containing this position
                        if chk_cov_snp_ql(read2, cov, m2_pos, read2_pos, bs_q, rref, snp_loc):
                            nucl = [pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in cov.pileups if m2_pos == cov.n]
                            m2_key = str(rref) + ":" + str(m2_pos)
                            sub_dict_read2[m2_key] = ref + rmm_b
                        else:
                            continue
                        if m2_key in sub_dict and m2_key in sub_dict_read2 and sub_dict[m2_key] == sub_dict_read2[m2_key]:
                            mismatch[m2_key] = ref + rmm_b
                            count = 2
                            depth = int(cov.n)
                            for pileupread in cov.pileups:
                                if not pileupread.is_del and not pileupread.is_refskip and m2_pos == cov.n and not pileupread.indel and rmm_b == \
                                        pileupread.alignment.query_sequence[pileupread.query_position] and \
                                        pileupread.alignment.query_sequence[pileupread.query_position] != ref:
                                    #       print ref,alt,pileupread.alignment.query_sequence[pileupread.query_position]
                                    count += 1
                                if 0 < float(float(count) / float(depth)) < 0.4:
                                    if float(float(count) / float(depth)) > 1:
                                        tmismatch[m2_key] = 1
                                    else:
                                        tmismatch[m2_key] = ref + rmm_b
                            del mismatch[m2_key]

        elif read_r2_1st_overlap(pos, mpos, rlen, mlen, md1, md2):
            ovrlp_seq = ovrlp_seq + (int(mpos) + int(mlen) - int(pos))
            mm_cnt = 0
            read_pos = 0  # used to catch substitution position
            read2_pos = 0
            m1_pos = int(pos)
            m2_pos = int(mpos)

            for idx3, align_type in enumerate(MDZ1):
                if ref_matches(align_type, idx3, rcig, mcig, MDZ1):
                    read_pos = read_pos + int(align_type)
                    m1_pos = int(m1_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig) and read.query_qualities[(int(read_pos))] >= bs_q:
                    ref = align_type
                    mm_cnt = mm_cnt + 1  # mutation per read count
                    m1_pos = m1_pos + 1
                    mm_b = str(read.query_alignment_sequence)[read_pos]  # mismatching base
                    m1_key = str(rref) + ":" + str(m1_pos)
                    sub_dict[m1_key] = ref + mm_b
                    read_pos = read_pos + 1

            for idx4, align_type in enumerate(MDZ2):
                if ref_matches(align_type, idx4, rcig, mcig, MDZ2):
                    read2_pos = read2_pos + int(align_type)
                    m2_pos = int(m2_pos) + int(align_type)
                elif ref_mismatches(align_type, nucs, rcig, mcig):
                    ref = align_type
                    mm_cnt = mm_cnt + 1
                    m2_pos = m2_pos + 1
                    rmm_b = str(read2.query_alignment_sequence)[read2_pos]  # mismatching base
                    for cov in SAM.pileup(str(rref), int(m2_pos) - 1, int(m2_pos)):  # this seems to be quite slow as it pulls all the reads containing this position
                        if chk_cov_snp_ql(read2, cov, m2_pos, read2_pos, bs_q, rref, snp_loc):
                            nucl = [pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in cov.pileups if m2_pos == cov.n]
                            m2_key = str(rref) + ":" + str(m2_pos)
                            sub_dict_read2[m2_key] = ref + rmm_b
                        else:
                            continue

                        if m2_key in sub_dict and m2_key in sub_dict_read2 and sub_dict[m2_key] == sub_dict_read2[m2_key]:
                            mismatch[m2_key] = ref + rmm_b
                            count = 2
                            depth = int(cov.n)
                            for pileupread in cov.pileups:
                                if not pileupread.is_del and m2_pos == cov.n and not pileupread.is_refskip and not pileupread.indel and rmm_b == \
                                        pileupread.alignment.query_sequence[pileupread.query_position] and \
                                        pileupread.alignment.query_sequence[pileupread.query_position] != ref:
                                    count += 1
                                if (float(float(count) / float(depth))) > 0 and float(
                                        float(count) / float(depth)) < 0.4:
                                    if (float(float(count) / float(depth))) > 1:
                                        tmismatch[m2_key] = 1
                                    else:
                                        tmismatch[m2_key] = ref + rmm_b
                            del mismatch[m2_key]
                        read2_pos += 1
            else:
                continue

        sub_dict_read2 = {}
        sub_dict = {}
        read_counts = read_counts + 1

        if int(read_counts) > int(num_reads):
            print("Number of paired reads processed that satisfy thresholds is {} ".format(read_counts), end="\r")
            sys.stdout.flush()
            num_reads += 1500000

print("Writing mutations to file \n")
with open(name + "/" + name + "_substitutions.out", "w") as out:
    for key, value in tmismatch.items():
        outstring = str(key) + "\t" + str(value) + "\n"
        out.write(outstring)

print("Calculating Mutation counts where Pair agree \n")
with open(name + "/" + name + "_sub_counts.out", "w") as cnts:
    subs = Counter(tmismatch.values())
    cnts.write("Sub Count Overlap\n")
    for key, value in subs.items():
        out = str(key) + " " + str(value) + " " + str(ovrlp_seq) + "\n"
        cnts.write(out)

print("making VCF for {}".format(name))
with open(name + "/" + name + "_substitutions.out", "r") as IN:
    with open(name + "/" + name + ".vcf", "w") as out:
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
            'assembly=b37,length=35477943>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(name))
        for line in IN:
            fields = line.strip().split()
            ref = list(fields[1])[0]
            mut = list(fields[1])[1]
            location = fields[0].strip().split(":")
            out_str = location[0] + "\t" + location[1] + "\t" + fields[
                0] + "\t" + ref + "\t" + mut + "\t.\t.\t.\t.\t1/0\n"
            out.write(out_str)
tme = (time.time() - strt) / 3600
print('overlap_mismatch has completed. The number of overlapping bases is {}\n \
    The number of putative mismatches is {}\n \
	The time taking to analyse {} was {} hrs'.format(ovrlp_seq, sum(subs.values()), name, tme))