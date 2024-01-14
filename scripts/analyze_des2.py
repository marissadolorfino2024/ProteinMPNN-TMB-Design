import subprocess
import glob
import os, itertools
import csv
import sys


env = sys.argv[1]

# define hydropathy scale
hydro_scale = {"A":1.800, "C":2.500, "D":-3.500, "E":-3.500, "F":2.800, "G":-0.400, "H":-3.200, "I":4.500, "K":-3.900, "L":3.800, "M":1.900, "N":-3.500, "P":-1.600, "Q":-3.500, "R":-4.500, "S":-0.800, "T":-0.700, "V":4.200, "W":-0.900, "Y":-1.300}


# overall sequence hydropathy
def GRAVY(sequence):
	tot_score = 0 
	for aa in sequence:
		tot_score += hydro_scale[aa]
	hydro_score = tot_score/len(sequence)

	return hydro_score

# define function to compute pore hydropathy of sequence 
def GRAVY_core(sequence, env):
	tot_score = 0
	n_resi = 0.0
	if env == 'sol':
		positions = [9,11,13,15,23,25,27,29,31,37,39,41,43,51,53,55,57,59,65,67,69,77,79,81,83,89,91,93,101,103,105,107]
	else:
		positions = [6,8,10,12,14,22,24,26,28,30,36,38,40,42,50,52,54,56,58,60,66,68,70,72,74,82,84,86,88,90,92,98,100,102,104,112,114,116,118,120]
	for pos in positions:
        	tot_score += hydro_scale[sequence[pos-1]]
        	n_resi +=1.0
	hydro_score = tot_score/n_resi
	
	return hydro_score

def GRAVY_surf(sequence, env):
	tot_score = 0
	n_resi = 0.0
	if env == 'sol':
		positions = [10,12,14,16,22,24,26,28,30,32,36,38,40,42,44,50,52,54,56,58,60,64,66,68,70,76,78,80,82,84,88,90,92,94,100,102,104,106,108]
	else:
		positions = [7,9,11,13,15,21,23,25,27,29,31,35,37,39,41,43,49,51,53,55,57,59,61,65,67,69,71,73,75,81,83,85,87,89,91,93,97,99,101,103,105,111,113,115,117,119,121]
	for pos in positions:
		tot_score += hydro_scale[sequence[pos-1]]
		n_resi += 1.0
	hydro_score = tot_score/n_resi

	return hydro_score 

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

# run pore hydropathy calcs after checking if both files exist
for file in glob.glob("*.fasta"):
	file = file
	if env == 'mem':
		ss3_file = file[:-6] + ".ss3"
		if os.path.exists(ss3_file):
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
			pore = GRAVY_core(design_seq, env)
			surf = GRAVY_surf(design_seq, env)
			full = GRAVY(design_seq)	
			# beta sheet prop
			prop = beta_prop(ss3_file)
			
			analysis_dict[file] = [pore, surf, full, prop]
	else:
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
		pore = GRAVY_core(design_seq, env)
		surf = GRAVY_surf(design_seq, env)
		full = GRAVY(design_seq)
		prop = 'na'
						
		analysis_dict[file] = [pore, surf, full, prop]

# append pore hydropathy and beta sheet prop to csv
with open('medref_soluble_analysis_11012023.csv', 'w') as csv:
	csv.write('design;'+'pore_hydro;'+'surf_hydro;'+'all_hydro;'+'beta_prop'+';'+'category'+'\n')
	for key, value in analysis_dict.items():
		csv.write(str(key)+';'+str(value[0])+';'+str(value[1])+';'+str(value[2])+';'+str(value[3])+';'+'medref_soluble'+'\n')
