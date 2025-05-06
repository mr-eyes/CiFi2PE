#!/usr/bin/python

import argparse
import os
import re
import matplotlib.pyplot as plt
import numpy as np
from Bio import Restriction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='FiC2PE', description='Converts HiFi FiC/CiFi reads to HiC paired end reads')

parser.add_argument("FastqFile", help="fastq file to split")
parser.add_argument("RestrictionEnzyme", help="restriction enzyme")
parser.add_argument("--out", default="out")
parser.add_argument("--version", action='version', version='%(prog)s 1.0')

args = parser.parse_args()

def fic2pe(ficread, restrictionenzyme):
    re = restrictionenzyme
    read = ficread
    lengthOfConcatenatedFragments = []
    if re == "HindIII":
        digested_record = list(Restriction.HindIII.catalyse(read.seq))
        for sequence in digested_record:
            lengthOfConcatenatedFragments.append(len(sequence))
        re_indices = Restriction.HindIII.search(read.seq)
        numberOfConcatenatedFragments = len(re_indices)
    elif re == "NlaIII":
        digested_record = list(Restriction.NlaIII.catalyse(read.seq))
        for sequence in digested_record:
            lengthOfConcatenatedFragments.append(len(sequence))
        re_indices = Restriction.NlaIII.search(read.seq)     
        numberOfConcatenatedFragments = len(re_indices)
    q_scores = read.letter_annotations
    digested_q_scores = []
    start = 0
    last_item = re_indices[-1]  
    for x,y in q_scores.items():
        for i in re_indices:
            if i != last_item:
                end = i-1
                sliced = y[start:end]
                digested_q_scores.append(sliced)
                start = i-1
            else:
                end = i-1
                sliced = y[start:end]
                digested_q_scores.append(sliced)
                sliced = y[end:]
                digested_q_scores.append(sliced)   
    record_id_seq = []
    record_seq_pairs = []
    record_q_score_pairs = []
    for i1, seq1 in enumerate(digested_record):
        for i2, seq2 in enumerate(digested_record[i1 + 1:]):
            New_record_id_seq = str(record.id + "_" + str(i1) + "_" + str(i2))
            R1_seq = Seq(seq1[-150:])
            R2_seq = Seq(seq2[5:155])
            if (len(R1_seq) == 150) and (len(R2_seq) == 150):
                record_id_seq.append(New_record_id_seq)
                record_seq_pairs.append([R1_seq, R2_seq])
            else:
                continue
    for j1, q1 in enumerate(digested_q_scores):
        for j2, q2 in enumerate(digested_q_scores[j1 + 1:]):
            R1_q = q1[-150:]
            R2_q = q2[5:155]
            if (len(R1_q) == 150) and (len(R2_q) == 150):
                record_q_score_pairs.append([R1_q, R2_q])
            else:
                continue
    list_of_R1_dicts = []
    list_of_R2_dicts = []
    for k, q_score in enumerate(record_q_score_pairs):
        R1_q_dict = {}
        R2_q_dict = {}
        R1_q_dict['phred_quality'] = q_score[0]
        R2_q_dict['phred_quality'] = q_score[1]
        list_of_R1_dicts.append(R1_q_dict)
        list_of_R2_dicts.append(R2_q_dict)
    R1s = []
    R2s = []
    for k, new_record in enumerate(record_id_seq):
        current_R1_seq_record = SeqRecord(Seq(record_seq_pairs[k][0]), id = new_record + " 1", letter_annotations = list_of_R1_dicts[k])
        current_R2_seq_record = SeqRecord(Seq(record_seq_pairs[k][1]), id = new_record + " 2", letter_annotations = list_of_R2_dicts[k])
        R1s.append(current_R1_seq_record)
        R2s.append(current_R2_seq_record)
    return(R1s, R2s, lengthOfConcatenatedFragments, numberOfConcatenatedFragments)

all_R1s = []
all_R2s = []
all_lengthOfConcatenatedFragments = []
all_numberOfConcatenatedFragments = []

sequencenumber=0

for record in SeqIO.parse(args.FastqFile, "fastq"):
    if args.RestrictionEnzyme == "HindIII":
        digested_record = list(Restriction.HindIII.catalyse(record.seq))
        re_indices = Restriction.HindIII.search(record.seq)
    elif args.RestrictionEnzyme == "NlaIII":
        digested_record = list(Restriction.NlaIII.catalyse(record.seq))
        re_indices = Restriction.NlaIII.search(record.seq)
    sequencenumber += 1
    # print("We're on " + str(sequencenumber))
    if len(digested_record) > 3:
        current_R1s, current_R2s, current_lengths, current_number = fic2pe(record, args.RestrictionEnzyme)
        all_R1s += current_R1s
        all_R2s += current_R2s
        all_lengthOfConcatenatedFragments += current_lengths
        all_numberOfConcatenatedFragments.append(current_number)
    else:
        continue

# Write paired files
R1_fastq_file = args.out + "_HiC_R1.fastq"
R2_fastq_file = args.out + "_HiC_R2.fastq"

with open(R1_fastq_file, "w") as output_handle_R1:
     SeqIO.write(all_R1s, output_handle_R1, "fastq")
with open(R2_fastq_file, "w") as output_handle_R2:
     SeqIO.write(all_R2s, output_handle_R2, "fastq") 

# Fragment counts histogram
bins = np.arange(0,3000,100)
plt.hist(np.clip(all_lengthOfConcatenatedFragments, bins[0], bins[-1]), bins=20)  # Adjust the number of bins as needed
plt.xlabel('Lengths of concatenated fragments (bp)')
plt.ylabel('Frequency')
plt.title(args.RestrictionEnzyme)

fraglenhistfile = args.out + "_fraglenhist.png"
plt.savefig(fraglenhistfile)

fraglenfile = args.out + "_fraglens.txt"
with open(fraglenfile, 'w') as f:
    for item in all_lengthOfConcatenatedFragments:
        f.write("%s\n" % item)

# Fragment counts histogram
bins = np.arange(0,100,5)
plt.hist(np.clip(all_lengthOfConcatenatedFragments, bins[0], bins[-1]), bins=20)  # Adjust the number of bins as needed
plt.xlabel('Number of concatenated fragments per HiFi read')
plt.ylabel('Frequency')
plt.title(args.RestrictionEnzyme)

fragcounthistfile = args.out + "_fragcounthist.png"
plt.savefig(fragcounthistfile)

fragcountfile = args.out + "_fragcounts.txt"
with open(fragcountfile, 'w') as f:
    for item in all_numberOfConcatenatedFragments:
        f.write("%s\n" % item)
