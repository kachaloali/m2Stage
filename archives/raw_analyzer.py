#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
"""
	Perform SNP analysis of the ec numbers predicted by the automatic enzyme annotation 
	programs. An ec number is considered :

	True Positive : when for a given protein this ec number appears in this protein both 
	at the reference level and also at the level of the predictions made by the prediction 
	program.

	False Positive : when for a given protein this ec number appears in this protein only 
	at the level of the predictions made by the prediction program.

	False Negative : when for a given protein this ec number appears in this protein only 
	at the reference level.

	True Negative : when, for a given protein, this ec number does not appear in this 
	protein neither at the reference level nor at the level of the predictions made by 
	the prediction program. Predictables ec numbers are those from the enzyme database

	Once this is done, for each ec number the sensitivity, precision and specificity are 
	calculated.
"""
import subprocess, os, pandas as pd
from scipy.stats import zscore
from Bio.ExPASy import Enzyme
from utils import *
import numpy as np

def get_enzyme_ecs(level): 
	'''
	Reads a Expasy Enzyme .dat file and writes a tab separated file where the first column is 
	EC number, the second column is the reaction description, the third column is the associated 
	uniprot ids separated by '|', and the fourth column indicates whether the reactions described 
	by this EC have been transferred to other EC numbers.
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
	True Negative. Once this is done, we then calculate the sensitivities, the specificities, and precision 
	associated to each ec number.
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
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
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
		num = l[0]
		den = l[0] + l[2]
		if den != 0: sens = round(float(num)/den, 8)
		else: sens = float(0)

		# compute specificity: TN / (TN + FP)
		#num = l[3]
		#den = l[3] + l[1]
		#if den != 0: spec = round(float(num)/den, 8)
		#else: spec = float(0)

		# compute precision: TP / (TP + FP)
		num = l[0]
		den = l[0] + l[1]
		if den != 0: prec = round(float(num)/den, 8)
		else: prec = float(0)

		# compute f1-score: 2*( sens * prec) / (sens + prec)
		num = 2*( sens * prec)
		den = sens + prec
		if den != 0: f1score = round(float(num)/den, 8)
		else: f1score = float(0)
		TP = l[0]
		FP = l[1]
		FN = l[2]
		dic[ec] = [TP, FP, FN, sens, prec, f1score]
	return dic				

def write_annot_ec_file(dic_global_ec, merged_dic, ref_dic, predictibles, level, output):
	dic_ecNumbers = {}
	ref_ec_list = [ref_dic[x] for x in ref_dic if "NA" not in ref_dic[x]]
	#print ref_ec_list
	ref_ec_list = list(set(reduce(lambda x, y: x + y, ref_ec_list)))
	ref_ec_list = [ec for ec in ref_ec_list	if len([x for x in ec.split(".") if x != "-"]) >= int(level)]
	ref_ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ref_ec_list]
	for prog, dic in dic_global_ec.items():				
		for ec, l in dic.items():
			#if l[5] > 0: # keep only ec with f1-score greater than 0
			if ec not in dic_ecNumbers: 
				dic_ecNumbers[ec] = dict()
				dic_ecNumbers[ec]["N"] = 0
			dic_ecNumbers[ec][prog] = l
			if ec in ref_ec_list: dic_ecNumbers[ec]["N"]+= 1
	new_ref_dic = dict()
	for protID, ec_list in ref_dic.items():
		if 4 - int(level) == 0:
			ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(level)]
		else:
			ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ec_list \
			if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
		new_ref_dic[protID] = ec_list		
				
	out = open(output, "w")
	out.write("{0} {1} {2} {3} {4} {5}\n".format("ProtID".ljust(12), "\tPROG".ljust(12), "\tECs".ljust(12), "\tPRED".ljust(06), "\tF1SCORE".ljust(12), "\tSP".ljust(12)))
	for protID in new_ref_dic:
		if protID in merged_dic:
			for prog in merged_dic[protID]:
				if len(merged_dic[protID][prog]) > 0:
					for ec in merged_dic[protID][prog]:
						if ec in new_ref_dic[protID]:
							f1score = dic_ecNumbers[ec][prog][5]
							out.write("{0} {1} {2} {3} {4} {5}\n".format(protID.ljust(12),"\t"+str(prog).ljust(12),"\t"+str(ec).ljust(12),"\t1".ljust(06),"\t"+str(f1score).ljust(12), "\t1".ljust(12)))
						else:
							f1score = dic_ecNumbers[ec][prog][5]
							out.write("{0} {1} {2} {3} {4} {5}\n".format(protID.ljust(12), "\t"+str(prog).ljust(12), "\t"+str(ec).ljust(12), "\t1".ljust(06), "\t"+str(f1score).ljust(12), "\t0".ljust(12)))
						for ec in new_ref_dic[protID]:
							if ec not in merged_dic[protID][prog]:
									f1score = dic_ecNumbers[ec][prog][5]
									out.write("{0} {1} {2} {3} {4} {5}\n".format(protID.ljust(12), "\t"+str(prog).ljust(12), "\t"+str(ec).ljust(12), "\t0".ljust(06), "\t"+str(f1score).ljust(12), "\t1".ljust(12)))
		else:
			for ec in new_ref_dic[protID]:
				if ec != "NA":
					if ec in dic_global_ec: 
						for prog in dic_global_ec[ec]:
							f1score = dic_ecNumbers[ec][prog][5]
							out.write("{0} {1} {2} {3} {4} {5}\n".format(protID.ljust(12), "\tNOT_FOUND".ljust(12), "\t"+str(ec).ljust(12), "\t0".ljust(06), "\t"+str(f1score).ljust(12), "\t1".ljust(12)))
	
	out.close()
		
if __name__ == '__main__':
	import warnings, os
	warnings.filterwarnings("ignore")
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--priam", dest = "priam",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted priam program file")
	parser.add_argument("--e2p2", dest = "e2p2",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted e2p2 program file")
	parser.add_argument("--b2g", dest = "b2g",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted Blast2Go program file")
	parser.add_argument("--kaas", dest = "kaas",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted kaas program file")
	parser.add_argument("--koala", dest = "koala",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted koala program file")
	parser.add_argument("--iprscan", dest = "iprscan",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted interproscan program file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--path", dest = "path",
						help = "Path to the limputils directory")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("--out",
		                dest = "output",
		                type = str,
		                default = "annotation_ec_file.tab",
		                help = "The output file name")
	args = parser.parse_args()
	if args.sp is None or args.level is None or args.path is None: 
		parser.print_help()
		exit(0)
	#
	#
	predictibles_ecs = get_enzyme_ecs(args.level)
	ref_dic = get_std_data(args.sp)
	dic_global_ec = dict()
	merged_dic = dict()
	# PRIAM
	if args.priam is not None:
		priam = get_prog_data(args.priam, args.path)
		dic_global_ec["PRIAM"] = get_dic_ec_score(predictibles_ecs, priam, ref_dic, args.level)
		
		for protID, ec_list in priam.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic: 
				merged_dic[protID]["PRIAM"] = set(ec_list)
			else: 
				merged_dic[protID] = dict()
				merged_dic[protID]["PRIAM"] = set(ec_list)
	# E2P2
	if args.e2p2 is not None:
		e2p2 = get_prog_data(args.e2p2, args.path)
		dic_global_ec["E2P2"] = get_dic_ec_score(predictibles_ecs, e2p2, ref_dic, args.level)
		
		for protID, ec_list in e2p2.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic:
				merged_dic[protID]["E2P2"] = set(ec_list)
			else:
				merged_dic[protID] = dict()
				merged_dic[protID]["E2P2"] = set(ec_list)
	# Blast2Go
	if args.b2g is not None:
		Blast2Go = get_prog_data(args.b2g, args.path)
		dic_global_ec["Blast2Go"] = get_dic_ec_score(predictibles_ecs, Blast2Go, ref_dic, args.level)
		
		for protID, ec_list in Blast2Go.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic:
				merged_dic[protID]["Blast2Go"] = set(ec_list)
			else:
				merged_dic[protID] = dict()
				merged_dic[protID]["Blast2Go"] = set(ec_list)
	# KAAS			
	if args.kaas is not None:
		kaas = get_prog_data(args.kaas, args.path)
		dic_global_ec["KAAS"] = get_dic_ec_score(predictibles_ecs, kaas, ref_dic, args.level)
		
		for protID, ec_list in kaas.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic:
				merged_dic[protID]["KAAS"] = set(ec_list)
			else:
				merged_dic[protID] = dict()
				merged_dic[protID]["KAAS"] = set(ec_list)
	# KOALA
	if args.koala is not None:
		koala = get_prog_data(args.koala, args.path)
		dic_global_ec["KOALA"] = get_dic_ec_score(predictibles_ecs, koala, ref_dic, args.level)
		
		for protID, ec_list in koala.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic:
				merged_dic[protID]["KOALA"] = set(ec_list)
			else:
				merged_dic[protID] = dict()
				merged_dic[protID]["KOALA"] = set(ec_list)
	# INTERPRO	
	if args.iprscan is not None:
		interproscan = get_prog_data(args.iprscan, args.path)
		dic_global_ec["INTERPRO"] = get_dic_ec_score(predictibles_ecs, interproscan, ref_dic, args.level)	
		
		for protID, ec_list in interproscan.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			if protID in merged_dic:
				merged_dic[protID]["INTERPRO"] = set(ec_list)
			else:
				merged_dic[protID] = dict()
				merged_dic[protID]["INTERPRO"] = set(ec_list)
	if len(dic_global_ec) < 1:
		parser.print_help()
		exit(0)
	ref_dic = get_std_data(args.sp)
	write_annot_ec_file(dic_global_ec, merged_dic, ref_dic, predictibles_ecs, args.level, args.output)
	
	
	
	
	
	
	
	
