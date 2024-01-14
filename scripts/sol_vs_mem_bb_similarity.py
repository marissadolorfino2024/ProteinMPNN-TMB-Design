#!/bin/usr/python3

import sys
import glob
import subprocess
import os, itertools
import csv

# script to calculate the RMSD and TM align scores between soluble and transmembrane input sets

sol_backbones = []
for sol in glob.glob("../water_soluble_raw_backbone/sol_*"):
    sol_backbones.append(sol)
    
mem_backbones = []
for mem in glob.glob("../transmembrane_raw_backbone/MR_inputs/*.pdb"):
    mem_backbones.append(mem)
    
to_write = [['rmsd', 'tmscore']]
tms = []

for sol_bb in sol_backbones:
    for mem_bb in mem_backbones:
        cmd = "/data/brussel/vo/000/bvo00014/bin/TMalign " + sol_bb + " " + mem_bb

        output = subprocess.check_output(["/data/brussel/vo/000/bvo00014/bin/TMalign", sol_bb, mem_bb])

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
                if "chain_1" in line:
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

        out_line = [RMSD, tmscore]
        to_write.append(out_line)
        tms.append(tmscore)

with open('MR_sol_v_mem_bbs.csv', 'w') as out_file:
    writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')

    for o_line in to_write:
        writer.writerow(o_line)

