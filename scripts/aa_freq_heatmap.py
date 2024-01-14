#!/usr/bin/python

import glob
import sys
import csv

outp = sys.argv[1]

# script to generate heatmap data -- each row will be of the form input aa designed aa percent designed (20 rows per each input aa)

# create dictionary of input sequences, identified by bb id
core_pos = [6, 8, 10, 12, 14, 22, 24, 26, 28, 30, 36, 38, 40, 42, 50, 52, 54, 56, 58, 60, 66, 68, 70, 72, 74, 82, 84, 86, 88, 90, 92, 98, 100, 102, 104, 112, 114, 116, 118, 120]
core_pos_py = [i - 1 for i in core_pos]

surf_pos = [7, 9, 11, 13, 15, 21, 23, 25, 27, 29, 31, 35, 37, 39, 41, 43, 49, 51, 53, 55, 57, 59, 61, 65, 67, 69, 71, 73, 75, 81, 83, 85, 87, 89, 91, 93, 97, 99, 101, 103, 105, 111, 113, 115, 117, 119, 121]
surf_pos_py = [i - 1 for i in surf_pos]

inp_dict = {}
for inp_bb in glob.glob('input_fastas/*.fasta'):
    bb = (str(inp_bb).split('bb_'))[1]
    bb = bb.split('.')[0]

    with open(inp_bb, 'r') as inp:
        lines = inp.readlines()
        seq = lines[1][:-1]

    inp_dict[bb] = seq

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
amino_counts = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
freq_dict = {}
# dict with 20 aas as keys, each containing a dict of aa counts
freq_dict = {aa: amino_counts.copy() for aa in aas}
freq_dict_core = {aa: amino_counts.copy() for aa in aas}
freq_dict_surf = {aa: amino_counts.copy() for aa in aas}

# iterate over all designs
for des in glob.glob('*.fasta'):
    bb = (str(des).split('align'))[1]
    bb = bb.split('.')[0]

    with open(des, 'r') as des_file:
        lines = des_file.readlines()
        if len(lines) < 2:
            continue
        else:
            des_seq = lines[1][:-1]

    des_seq = list(des_seq)
    ref_seq = list(inp_dict[bb])

    for i in range(len(ref_seq)):
        ref_aa = ref_seq[i]
        des_aa = des_seq[i]

        freq_dict[str(ref_aa)][str(des_aa)] += 1

        if i in core_pos_py:
            freq_dict_core[str(ref_aa)][str(des_aa)] += 1

        elif i in surf_pos_py:
            freq_dict_surf[str(ref_aa)][str(des_aa)] += 1

for aa in freq_dict.keys():
    freqs = freq_dict[aa]
    count_sum = 0.0
    for aa2 in freqs.keys():
        count_sum += float(freqs[aa2])

    for aa2 in freqs.keys():
        try:
            freqs[aa2] = round(float((float(freqs[aa2]) / float(count_sum))), 3)
        except:
            freqs[aa2] = float(0.0)

    freq_dict[aa] = freqs

for aa in freq_dict_core.keys():
    freqs_core = freq_dict_core[aa]
    count_sum = 0.0
    for aa2 in freqs_core.keys():
        count_sum += float(freqs_core[aa2])

    for aa2 in freqs_core.keys():
        try:
            freqs_core[aa2] = round(float((float(freqs_core[aa2]) / float(count_sum))), 3)
        except:
            freqs_core[aa2] = float(0.0)

    freq_dict_core[aa] = freqs_core

for aa in freq_dict_surf.keys():
    freqs_surf = freq_dict_surf[aa]
    count_sum = 0.0
    for aa2 in freqs_surf.keys():
        count_sum += float(freqs_surf[aa2])

    for aa2 in freqs_surf.keys():
        try:
            freqs_surf[aa2] = round(float((float(freqs_surf[aa2]) / float(count_sum))), 3)
        except:
            freqs_surf[aa2] = float(0.0)

    freq_dict_surf[aa] = freqs_surf

out_lines = [['ref_aa', 'des_aa', 'des_freq']]
out_lines_core = [['ref_aa', 'des_aa', 'des_freq']]
out_lines_surf = [['ref_aa', 'des_aa', 'des_freq']]

for aa in freq_dict.keys():
    freqs = freq_dict[aa]
    freqs_core = freq_dict_core[aa]
    freqs_surf = freq_dict_surf[aa]

    for aa2 in aas:
        aa_freq = freqs[aa2]
        out_line = [aa, aa2, aa_freq]
        out_lines.append(out_line)

        aa_freq_core = freqs_core[aa2]
        out_line_core = [aa, aa2, aa_freq_core]
        out_lines_core.append(out_line_core)

        aa_freq_surf = freqs_surf[aa2]
        out_line_surf = [aa, aa2, aa_freq_surf]
        out_lines_surf.append(out_line_surf)

with open(outp, 'w') as freq_file:
    writer = csv.writer(freq_file, delimiter='\t', lineterminator='\n')
    for out_line in out_lines:
        writer.writerow(out_line)

core_outp = outp[:-4] + '_core.csv'
surf_outp = outp[:-4] + '_surf.csv'

with open(core_outp, 'w') as core_file:
    writer = csv.writer(core_file, delimiter='\t', lineterminator='\n')
    for out_line in out_lines_core:
        writer.writerow(out_line)


with open(surf_outp, 'w') as surf_file:
    writer = csv.writer(surf_file, delimiter='\t', lineterminator='\n')
    for out_line in out_lines_surf:
        writer.writerow(out_line)


