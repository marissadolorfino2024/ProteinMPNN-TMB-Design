#!/bin/usr/python3

import sys
import glob
import subprocess
import os, itertools
import csv
import math

# same as tm_to_inp.py but doesn't average scores 

# script to calculate tm align score between af folded proteinmpnn sequences and input backbone

# should calculate mean tm align to original backbone for each (should be 10 designs per input backbone)

# 9_mpnn_align46_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000_0001.pdb

# input_bbs/fulldes_bb_28.pdb

# make dict of bb_ids and associated pdb file
inputs = {}
input_seqs = {}
for bb in glob.glob("input_bbs/*.pdb"):
    bb_id = str(bb)[21:-4]
    inputs[bb_id] = bb
    fa_file = str(bb)[:-3]+'fasta'
        
    with open(fa_file, 'r') as fa:
        lines = fa.readlines()
        seq = lines[1][:-1]

    input_seqs[bb] = seq

# make dict of bb_ids and associated list of designed pdb files (list of 10)
designs = {}
for des in glob.glob("*unrelaxed*"):
    bb_id = str(des)[12:-60]
    if "n" in bb_id:
        bb_id = bb_id[1:]
    else:
	bb_id = bb_id

    if bb_id in designs.keys():
	designs[bb_id].append(des)
    else:
	designs[bb_id] = [des]

to_write = [['des', 'tmscore', 'RMSD', 'avg_plddt', 'seq_diversity']]

for bb_id in inputs.keys():
    inp = inputs[bb_id]
    tms = []
    RMSDs = []
    inp_seq = input_seqs[bb_id]

    for des in designs[bb_id]:
	des_fasta = des.split('_unrelaxed')[0]
        des_fasta = str(des_fasta) + '.fasta'
                
        with open(des_fasta, 'r') as des_fa:
            lines = des_fa.readlines()
            if len(lines) <= 1:
                continue
                    
            elif len(lines) > 1:
                des_seq = lines[1][:-1]

                diversity = 0

                inp_seq = list(inp_seq)
                des_seq = list(des_seq)
                for i in range(len(inp_seq)):
                    inp_res = inp_seq[i]
                    des_res = des_seq[i]
                    if inp_res == des_res:
                        continue
                    elif inp_res != des_res:
                        diversity += 1

                    seq_diversity = (diversity / len(inp_seq))

                    cmd = "/data/brussel/vo/000/bvo00014/bin/TMalign " + inp + " " + des
                    output = subprocess.check_output(["/data/brussel/vo/000/bvo00014/bin/TMalign", inp, des])
                    
                    output_lines = output.decode().split("\n")
                    longer_chain = 0
                    RMSD = 0
                    tmscore = 0
                    
                    for line in output_lines:
	        	if "Length of Chain_1" in line:
	        	    chain_1_L = int(line.split()[3])
	        	elif "Length of Chain_2" in line:
	        	    chain_2_L = int(line.split()[3])
	        	elif "TM-score" in line:
	        	    if "Chain_1" in line:
	        		chain_1_tm = float(line.split()[1])
	        	    elif "Chain_2" in line:
	        		chain_2_tm = float(line.split()[1])
	        	elif "RMSD" in line:
	        	    RMSD = float(line.split(',')[1].replace("RMSD=","").lstrip())


                    
                    if chain_1_L >= chain_2_L:
                    	longer_chain = "Chain_1"
                    else:
                    	longer_chain = "Chain_2"
                    
                    if longer_chain == "Chain_1":
                    	tmscore = chain_1_tm
                    else:
                    	tmscore = chain_2_tm


                    
                    tms.append(tmscore)
                    RMSDs.append(RMSD)
                    
                    json = str(des[:-9]) + '.json'
                    json = json.replace('unrelaxed', 'scores')
                    with open(json, 'r') as json:
                        score_dict = json.readlines()
                        score_dict = score_dict[0]
                        plddts = score_dict.split('],')[0]
                        plddts = plddts.split('[')[1]
                        plddts = plddts.replace(' ', '')
                        plddts = plddts.split(',')
                        p_sum = 0
                        for plddt in plddts:
                        	plddt = float(plddt)
                        	p_sum += plddt
                        
                        avg_plddt = (p_sum / (len(plddts)))

	            to_write.append([des, tmscore, RMSD, avg_plddt, seq_diversity])
	

with open('fullydesign_seqsim_vs_tm.csv', 'w') as outp:
    writer = csv.writer(outp, delimiter = '\t', lineterminator='\n')
    for out in to_write:
	writer.writerow(out)

