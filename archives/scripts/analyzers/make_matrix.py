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
from Bio.ExPASy import Enzyme
from utils import *
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


def get_snp(predictibles, prog_pred_dic, ref_dic, level):
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
				#list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog]
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
				#list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
			elif (4 - int(level) == 0):
				#list_ECs_prog_cl = list_ECs_prog
				list_ECs_prog_cl = [ec for ec in list_ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				#list_ECs_uniprot_cl = list_ECs_uniprot
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
				#list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog]
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
			else: 
				#list_ECs_prog_cl = list_ECs_prog
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
				#list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
			else: 
				#list_ECs_uniprot_cl = list_ECs_uniprot
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
	return dic_ecNumbers
	
	
if __name__ == '__main__':
	import os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--priam", dest = "priam", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--e2p2", dest = "e2p2", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--kaas", dest = "kaas", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--koala", dest = "koala", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--iprscan", dest = "iprscan", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--path", dest = "path",
						help = "Path to the limputils directory")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	args = parser.parse_args()
	if args.sp is None or args.level is None or args.path is None: 
		parser.print_help()
		exit(0)
	
	#
	#
	ref_dic = get_sp_data(args.sp)
	predictibles_ecs = get_enzyme_ecs(args.level)
	
	l0 = []
	d0 = {}
	all_ECs = set()
	if args.priam is not None:
		l0.append("PRIAM")
		priam = get_prog_data(args.priam)
		d0["PRIAM"] = get_snp(predictibles_ecs, priam, ref_dic, args.level)
		all_ECs.update(d0["PRIAM"].keys())
	if args.e2p2 is not None:
		l0.append("E2P2")
		e2p2 = get_prog_data(args.e2p2)
		d0["E2P2"] = get_snp(predictibles_ecs, e2p2, ref_dic, args.level)
		all_ECs.update(d0["E2P2"].keys())
	if args.kaas is not None:
		l0.append("KAAS")
		kaas = get_prog_data(args.kaas)
		d0["KAAS"] = get_snp(predictibles_ecs, kaas, ref_dic, args.level)
		all_ECs.update(d0["KAAS"].keys())
	if args.koala is not None:
		l0.append("KOALA")
		koala = get_prog_data(args.koala)
		d0["KOALA"] = get_snp(predictibles_ecs, koala, ref_dic, args.level)
		all_ECs.update(d0["KOALA"].keys())
	if args.iprscan is not None:
		l0.append("IPRSACN")
		iprscan = get_prog_data(args.iprscan)
		d0["IPRSACN"] = get_snp(predictibles_ecs, iprscan, ref_dic, args.level)
		all_ECs.update(d0["IPRSACN"].keys())
	if len(d0) == 0:
		parser.print_help()
		exit(0)
	else:
		fw = open("out.csv", "w")
		title = "EC"
		for elt in l0: 
			title+= "\t"
			title+= elt
		title+= "\n"
		fw.write(title) 
		for ec in all_ECs:
			line = ec
			for elt in l0:
				if ec in d0[elt].keys():
					line+= "\t"
					line+= str(d0[elt][ec][0])					
				else:
					line+= "\t"
					line+= "0"
			line+= "\n"
			fw.write(line)
		fw.close()
		
	
