import os
import re
import glob
import sys
import time
import math
from collections import Counter
from collections import defaultdict
import pysam

strt = time.time()

# Files ###########

# Running get_overlaps to estimate mismatches for sample
# The G -> T imbalance is 1.094878523269532
#  % of Read 1 mapping to forward strand is

## File structure is PUR/
#                             |
#                             logs/
#                                    |
#                                    slurm.xxxxx.log
#                             |
#                             |
#                             HGX1/
#                             HGX2/


#def avg_ovl():


for filepath in glob.iglob('/data/dbennett/Overlap/1KG/PUR/logs/slurm*'):
    with open(filepath, 'r') as IN:
        for line in IN:
            if line.startswith('Running get_overlaps to estimate mismatches for sample'):
                sample_id = line.strip().split()[-1]
                num_frags = 0
                ovl = 0
                with open('/data/dbennett/Overlap/1KG/PUR/' + sample_id+ '/' + sample_id + '.vcf','r') as VCF:
                    for entry in VCF:
                        if entry.startswith('#'):
                            continue
                        else:
                            fields=entry.strip().split("\t")
                            DP = fields[7].split(";")[0].split("=")[1]
                            OVL_len = fields[7].split(";")[3].split("=")[1].split(",")
                            num_frags += len(OVL_len)
                            ovl += sum(OVL_len)
            elif line.startswith('The G -> T imbalance is '):
                Giv = line.strip().split()[-1]
            elif line.startswith(' % of Read 1 mapping to forward strand is'):
                mp_imb = line.strip().split()[-1]
        print(sample_id,Giv,mp_imb,sep="\t")
