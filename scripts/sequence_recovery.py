#!/usr/bin/python

import glob
import sys
import csv

outp = sys.argv[1]

# script to calculate the sequence recovery at different regions (core, surface, all)

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

recovs = [['seq_recov']]
core_recovs = [['seq_recov']]
surf_recovs = [['seq_recov']]

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
    
    all_rec = 0 
    core_rec = 0
    surf_rec = 0
    
    for i in range(len(ref_seq)):
        ref_aa = ref_seq[i]
        des_aa = des_seq[i]

        if des_aa == ref_aa:
            all_rec += 1
        else:
            continue

        if i in core_pos_py:
            if des_aa == ref_aa:
                core_rec += 1
            else: 
                continue

        elif i in surf_pos_py:
            if des_aa == ref_aa:
                surf_rec += 1
            else:
                continue

    all_rec = [((all_rec) / (len(ref_seq))) * 100]
    recovs.append(all_rec)

    core_rec = [((core_rec) / (len(core_pos_py))) * 100]
    core_recovs.append(core_rec)

    surf_rec = [((surf_rec) / (len(surf_pos_py))) * 100]
    surf_recovs.append(surf_rec)

with open(outp, 'w') as freq_file:
    writer = csv.writer(freq_file, delimiter='\t', lineterminator='\n')
    for out_line in recovs:
        writer.writerow(out_line)

core_outp = outp[:-4] + '_core.csv'
surf_outp = outp[:-4] + '_surf.csv'

with open(core_outp, 'w') as core_file:
    writer = csv.writer(core_file, delimiter='\t', lineterminator='\n')
    for out_line in core_recovs:
        writer.writerow(out_line)


with open(surf_outp, 'w') as surf_file:
    writer = csv.writer(surf_file, delimiter='\t', lineterminator='\n')
    for out_line in surf_recovs:
        writer.writerow(out_line)


