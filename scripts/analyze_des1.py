import subprocess
import glob
import os, itertools
import csv
import sys

# define hydropathy scale
hydro_scale = {"A":1.800, "C":2.500, "D":-3.500, "E":-3.500, "F":2.800, "G":-0.400, "H":-3.200, "I":4.500, "K":-3.900, "L":3.800, "M":1.900, "N":-3.500, "P":-1.600, "Q":-3.500, "R":-4.500, "S":-0.800, "T":-0.700, "V":4.200, "W":-0.900, "Y":-1.300}

# define function to compute pore hydropathy of sequence 
def GRAVY_core(sequence):
	tot_score = 0
    	n_resi = 0.0
    	positions = [6,8,10,12,14,22,24,26,28,30,36,38,40,42,50,52,54,56,58,60,66,68,70,72,74,82,84,86,88,90,92,98,100,102,104,112,114,116,118,120]
    	for pos in positions:
        	tot_score += hydro_scale[sequence[pos-1]]
        	n_resi +=1.0
    	hydro_score = tot_score/n_resi
    	return(hydro_score)

# define function to compute beta sheet prop of designs
def beta_prop(ss3):
        ref = 'LLLLLLEEEEEEELLLLLLLEEEEEEEEEEELLLEEEEELLLLLLLLLLLEEEEEEEEEEELLLEEEEEEEEEEEEEEEEEEEEEEEEEEEEELLLEEEEELLLLLLLLLLLEEEEEEEEELL'
        len_ref = len(ref)

        ref_beta = 0
        for i in ref:
                if i == 'L':
                        pass
                else:
                        ref_beta += 1

        ss3_beta = 0
        ss3_pred = []
        with open(ss3, 'r') as raptor:
                for line in raptor:
                        if line[0] == '#':
                                pass
                        else:
                                pred = line[7]
                                ss3_pred.append(pred)
        for i in range(len_ref):
                ss3 = ss3_pred[i]
                ref_pred = ref[i]
                if ref_pred == 'E':
                        if ss3 == 'E':
                                ss3_beta += 1
                        else:
                                pass
                else:
                        pass

        propensity = (((float(ss3_beta))/(float(ref_beta)))) * 100
        return(propensity)

analysis_dict = {}

# run pore hydropathy calculations on all fastas
for file in glob.glob("*.fasta"):
	design_seq = []
  	with open(file, 'r') as fasta:
    		for line in fasta:
      			if line[0] == '>':
        			pass
      			else:
        			for i in line:
          				if i == '\n' or i == ' ':
            					pass
          				else:
            					design_seq.append(i)
	pore = GRAVY_core(design_seq)

# run beta sheet propensity on corresponding ss3 files
	ss3_file = file[:-6] + ".ss3"
	prop = beta_prop(ss3_file)
	
	analysis_dict[file] = [pore, prop, design_seq]

# append pore hydropathy and beta sheet prop to csv
with open('bb_analysis.csv', 'w') as csv:
	for key, value in analysis_dict.items():
		csv.write(str(key)+';'+str(value[0])+';'+str(value[1])+';'+str(value[2])+'\n')
