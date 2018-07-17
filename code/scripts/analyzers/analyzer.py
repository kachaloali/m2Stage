#!/usr/bin/python2
# <!- -*- coding: utf-8 -*- ->

#	Perform SNP analysis of the ec numbers predicted by the automatic enzyme annotation 
#	programs. An ec number is considered :

#	True Positive : when for a given protein this ec number appears in this protein both 
#	at the reference level and also at the level of the predictions made by the prediction 
#	program.

#	False Positive : when for a given protein this ec number appears in this protein only 
#	at the level of the predictions made by the prediction program.

#	False Negative : when for a given protein this ec number appears in this protein only 
#	at the reference level.

#	True Negative : when, for a given protein, this ec number does not appear in this 
#	protein neither at the reference level nor at the level of the predictions made by 
#	the prediction program. Predictables ec numbers are those from the enzyme database

#	Once this is done, for each ec number the sensitivity, precision and specificity are 
#	calculated.

import subprocess, os, pandas as pd
from scipy.stats import zscore
from Bio.ExPASy import Enzyme
from utils import *
import numpy as np, sys

def get_enzyme_ecs(level): 
	'''
	Reads a Expasy Enzyme .dat file and writes a tab separated file where the first column is 
	EC number, the second column is the reaction description, the third column is the associated 
	uniprot ids separated by '|', and the fourth column indicates whether the reactions described 
	by this EC have been transferred to other EC numbers.
	
	It return the list of the whole ecs numbers of ENZYME database if the latter is complete
	'''
	if not os.path.exists(os.path.join("database", "enzyme", "enzyme.dat")):
		curl_enzyme = os.path.join("ftp://ftp.expasy.org", "databases","enzyme", "enzyme.dat")
		subprocess.check_output("wget -cq -P database/enzyme " + curl_enzyme, shell = True)
	if not os.path.exists(os.path.join("database", "enzyme", "enzyme.dat")): 
		print ("%s\n", "Missing enzyme database!")
		exit(0)
	input_name = os.path.join("database", "enzyme", "enzyme.dat")
	output_name = os.path.join("database", "enzyme", "enzyme.tsv")
	records = Enzyme.parse(open(input_name))

	out = dict() # dict of dicts, first key: EC number, second key: field
	transferred = dict() #dict of lists
	for record in records:
		if 'Transferred entry:' in record['DE']:
			record['DE'] = record['DE'].rstrip('.')
			record['DE'] = record['DE'].replace('Transferred entry:',' ')
			record['DE'] = record['DE'].replace(',',' ')
			record['DE'] = record['DE'].replace('and',' ')
			point_to = record['DE'].split()
			transferred[record['ID']] = point_to
		else:
			out[record['ID']] = dict()
			out[record['ID']]['uniprot'] = '|'.join([x[0] for x in record['DR']])
			out[record['ID']]['description'] = record['DE']
			out[record['ID']]['transferred'] = False
	for id in transferred:
		out[id] = dict()
		out[id]['uniprot'] = '|'.join([out[x]['uniprot'] for x in transferred[id]])
		out[id]['description'] = 'Transferred entry: ' + ' '.join(transferred[id])
		out[id]['transferred'] = True
	df = pd.DataFrame.from_dict(out, orient = 'index')
	df.index.name = 'EC'
	
	# write all data in a enzyme.csv file
	df.to_csv(output_name, sep = '\t')
	# ignore EC numbers with no uniprot ids associated
	df.dropna(subset = ['uniprot'], inplace = True)
	# ignore EC numbers that are obsolete due to transfer 
	df = df[df.transferred == False]
	all_ECs = list(set(df.index.values))
	if 4 - int(level) == 0:
		all_ECs = [ec for ec in all_ECs if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		all_ECs = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in all_ECs \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	return list(set(all_ECs))

def get_dic_ec_score(predictibles, prog_pred_dic, ref_dic, level):
	"""
	This function makes it possible to perform SNP analysis of the ec numbers predicted by the program
	For each ec number, we enumerates all the cases it appeared in the proteins annotation. Then we try to 
	distinguish how many times it appears True Positive, False Positive, False Negative and finally 
	True Negative. Once this is done, we then calculate sensitivity, precision and the f1-score
	associated to each ec number.
	
	It return a dictionnary of the kind: key: ec number, value: [TP, FP, FN, sens, prec, f1score]
	"""
	dic_ecNumbers = {}
	if (4 - int(level) > 0):
		predictibles_ecs = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in predictibles \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	else:
		predictibles_ecs = [ec for ec in predictibles 
		if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	
	for protID in ref_dic:
		if protID in prog_pred_dic and prog_pred_dic[protID] != ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_prog = prog_pred_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
			elif (4 - int(level) == 0):
				list_ECs_prog_cl = [ec for ec in list_ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				list_ECs_uniprot_cl = [ec for ec in list_ECs_uniprot 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			
			TP, FP, FN = 0, 0, 0
			characterized_ecs = set()
			for ec in list_ECs_prog_cl:
				if ec in list_ECs_uniprot_cl:
					if ec in dic_ecNumbers: 
						dic_ecNumbers[ec][0]+= 1
						TP+= 1
					else: 
						dic_ecNumbers[ec] = [1,0,0,0]
						TP+= 1
				else:
					if ec in dic_ecNumbers: 
						dic_ecNumbers[ec][1]+= 1
						FP+= 1
					else: 
						dic_ecNumbers[ec] = [0,1,0,0]
						FP+= 1
				characterized_ecs.add(ec)
			
			for ec in list_ECs_uniprot_cl:
				if ec not in list_ECs_prog_cl and ec in predictibles_ecs:
					if ec in dic_ecNumbers: 
						dic_ecNumbers[ec][2]+= 1
						FN+= 1
					else: 
						dic_ecNumbers[ec] = [0,0,1,0]
						FN+= 1
					characterized_ecs.add(ec)
			# computing TN
			for ec in characterized_ecs:
				if ec in predictibles_ecs: 
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP + FN)
				else:
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP) + FN
			
		elif protID in prog_pred_dic and prog_pred_dic[protID] != ['NA'] and ref_dic[protID] == ['NA']:
			list_ECs_prog = prog_pred_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
			else: 
				list_ECs_prog_cl = [ec for ec in list_ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				
			TP, FP, FN = 0, 0, 0
			characterized_ecs = set()
			for ec in list_ECs_prog_cl:
				if ec in dic_ecNumbers: 
					dic_ecNumbers[ec][1]+= 1
					FP+= 1
				else: 
					dic_ecNumbers[ec] = [0,1,0,0]
					FP+= 1
				characterized_ecs.add(ec)	
			# computing TN
			for ec in characterized_ecs:
				if ec in predictibles_ecs: 
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP + FN)
				else:
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP) + FN
					
		elif protID in prog_pred_dic and prog_pred_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl))
			else: 
				list_ECs_uniprot_cl = [ec for ec in list_ECs_uniprot 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			
			TP, FP, FN = 0, 0, 0
			characterized_ecs = set()
			for ec in list_ECs_uniprot_cl:
				if ec in predictibles_ecs:
					if ec in dic_ecNumbers: 
						dic_ecNumbers[ec][2]+= 1
						FN+= 1
					else: 
						dic_ecNumbers[ec] = [0,0,1,0]
						FN+= 1	
					characterized_ecs.add(ec)
			# computing TN
			for ec in characterized_ecs:
				if ec in predictibles_ecs: 
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP + FN)
				else:
					dic_ecNumbers[ec][3]+= len(predictibles_ecs) - (TP + FP) + FN
	dic = dict()	
	for ec, l in dic_ecNumbers.items():
		# compute sensitivity: TP / (TP + FN) 
		num, den = l[0], l[0] + l[2]
		if den != 0: sens = round(float(num)/den, 8)
		else: sens = float(0)

		# compute specificity: TN / (TN + FP)
		#num = l[3]
		#den = l[3] + l[1]
		#if den != 0: spec = round(float(num)/den, 8)
		#else: spec = float(0)

		# compute precision: TP / (TP + FP)
		num, den = l[0], l[0] + l[1]
		if den != 0: prec = round(float(num)/den, 8)
		else: prec = float(0)

		# compute f1-score: 2*( sens * prec) / (sens + prec)
		num, den = 2*( sens * prec), sens + prec
		if den != 0: f1score = round(float(num)/den, 8)
		else: f1score = float(0)
		TP, FP, FN = l[0], l[1], l[2]
		dic[ec] = [TP, FP, FN, sens, prec, f1score]
	return dic				

def write_ecs_to_pick(dic_global_ec, ref_dic, level, outfile):
	"""
	This function writes the file of the ecs numbers to pick during the process of the annotation of the new 
	data of sequences
	"""
	dic_ecNumbers = {}
	ref_ec_list = [ref_dic[x] for x in ref_dic if "NA" not in ref_dic[x]]
	ref_ec_list = list(set(reduce(lambda x, y: x + y, ref_ec_list)))
	ref_ec_list = [ec for ec in ref_ec_list	if len([x for x in ec.split(".") if x != "-"]) >= int(level)]
	if (4-int(level)> 0): ref_ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ref_ec_list]
	for prog, dic in dic_global_ec.items():				
		for ec, l in dic.items():
			if l[5] > 0: # keep only ec with f1-score greater than 0
				if ec not in dic_ecNumbers: 
					dic_ecNumbers[ec] = dict()
					dic_ecNumbers[ec]["N"] = 0
				dic_ecNumbers[ec][prog] = l
				if ec in ref_ec_list: dic_ecNumbers[ec]["N"]+= 1
	
	out = open(outfile, "w")
	out.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format("EC".ljust(12), "\tPROG".ljust(12), 
	"\tTP".ljust(06), "\tFP".ljust(06), "\tFN".ljust(06),"\tSENS".ljust(12), "\tPREC".ljust(12), 
	"\tF1SCORE".ljust(12), "\tN".ljust(06), "\tZSCORE".ljust(16)))
	for ec, dic in dic_ecNumbers.items():
		list_f1score = [float(dic[prog][5]) for prog in dic if prog != "N"]
		list_zscore = zscore(list_f1score)
		# Pour une distribution de notes de moyenne M et d'écart-type E, la note z-score correspondant au score x 
		# se calcule facilement avec la formule : Z_x = (x-M) / E
		dic_zscore, j = dict(), 0
		for prog in dic:
			if prog != "N":
				dic_zscore[prog] = list_zscore[j]
				j+= 1
		N = dic["N"]
		for prog, l in dic.items():
			if prog != "N":
				if np.isnan(dic_zscore[prog]): zscore_value = float(0)
				else: zscore_value = round(float(dic_zscore[prog]), 8) 
				out.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(str(ec).ljust(12), 
				"\t" + prog.ljust(12),"\t" +str(l[0]).ljust(06),"\t" + str(l[1]).ljust(06), 
				"\t" + str(l[2]).ljust(06), "\t"+str(l[3]).ljust(12),"\t"+str(l[4]).ljust(12),
				"\t"+str(l[5]).ljust(12), "\t"+str(N).ljust(06), "\t"+str(zscore_value).ljust(16)))
	out.close()

def write_file_for_heatmap(dic_global_ec, ref_dic, level, choice):
	"""
	This function writes the file of the mesures values for the programs to build heatmap.
	"""
	dic_ecNumbers = {}
	ref_ec_list = [ref_dic[x] for x in ref_dic if "NA" not in ref_dic[x]]
	ref_ec_list = list(set(reduce(lambda x, y: x + y, ref_ec_list)))
	ref_ec_list = [ec for ec in ref_ec_list	if len([x for x in ec.split(".") if x != "-"]) >= int(level)]
	if (4-int(level)> 0): ref_ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ref_ec_list]
	for prog, dic in dic_global_ec.items():				
		for ec, l in dic.items():
			if l[5] > 0: # keep only ec with f1-score greater than 0
				if ec not in dic_ecNumbers: 
					dic_ecNumbers[ec] = dict()
					dic_ecNumbers[ec]["N"] = 0
				dic_ecNumbers[ec][prog] = l
				if ec in ref_ec_list: dic_ecNumbers[ec]["N"]+= 1
	#			
	if choice == 1:		out = open("out.heatmap.sens."+ level +".txt", "w")
	elif choice == 2:	out = open("out.heatmap.prec."+ level +".txt", "w")
	elif choice == 3:	out = open("out.heatmap.f1score."+ level +".txt", "w")
	progs, header = list(), "ECs"
	for key in dic_global_ec: 
		header+= "\t"
		header+= str(key)
	out.write(header + "\n")
	for ec, dic in dic_ecNumbers.items():
		line = str(ec)
		for key in dic_global_ec:
			if key in dic_ecNumbers[ec]:
				line+= "\t"
				line+= str(dic_ecNumbers[ec][key][int(choice) + 2]).replace(".", ",")
			else: line+= "\t0"
		out.write(line + "\n")
	out.close()
		
if __name__ == '__main__':
	import warnings, os
	warnings.filterwarnings("ignore")
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	
	parser.add_argument("conv_files", type = str, nargs = '*',
						help = "The list of the converted programs files. Each file must be preceded by the\
						 name of the corresponded program (ex: E2P2 e2p2_file.conv PRIAM priam_file.conv)")    
	parser.add_argument("-std", dest = "std",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted std file")
	parser.add_argument("-plud", dest = "plud",
						help = "plud to the limputils directory")
	parser.add_argument("-level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("-out",
		                dest = "output", default = "ecs_to_pick.tab", help = "The output file name")
	parser.add_argument("-mes", type = int, choices = [1, 2, 3],
						help = "To produce the file of the mesure values for the programs to build heatmap."
							 " 1: Sensitivity; 2: Precision; 3: F1-Score")
	args = parser.parse_args()
	if args.std is None or args.level is None or args.plud is None or len(args.conv_files) % 2 != 0: 
		parser.print_help()
		exit(0)
	conv_files = {args.conv_files[i]: args.conv_files[i+1] for i in range(0,len(args.conv_files), 2)}
	for f in conv_files.values(): is_valid_file(parser, f)
	#
	predictibles_ecs = get_enzyme_ecs(args.level)
	ref_dic = get_std_data(args.std)
	dic_global_ec = dict()
	for key, value in conv_files.items():
		prog = get_prog_data(value, args.plud)
		dic_global_ec[key] = get_dic_ec_score(predictibles_ecs, prog, ref_dic, args.level)
	
	# writing the file of the ecs numbers to pick when annotating new sub family dataset.
	write_ecs_to_pick(dic_global_ec, ref_dic, args.level, args.output)
	# writing the file for constructing heatmap 
	if args.mes is not None: write_file_for_heatmap(dic_global_ec, ref_dic, args.level, args.mes)
