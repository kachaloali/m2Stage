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

def write_list_ecNumbers(prog_file, sp_file, lipmutils, level):
	# try to transfer ec numbers with lipmutils program
	lipmutils_program = os.path.join(lipmutils, "bin", "lipm_valid_ec_numbers.pl")
	command = [lipmutils_program, "--infile", prog_file, "--outfile", str(prog_file)+".transf"]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# try to write the ec numbers list in a file after being transferred
	in_file = str(prog_file) + ".transf"
	ec_list = [x.split("\t")[1].strip(" \n\t\r") for x in open(in_file, "r") if "NA" not in x]
	infile_extension = prog_file.split(".")[-1]
	outfile = prog_file[:len(prog_file) - len(infile_extension) - 1] + ".ec."+ level +".list"
	
	if 4 - int(level) == 0:
		ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ec_list \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	
	file = open(outfile, "w")
	for ec in set(ec_list): file.writelines(list("%s\n" % item for item in ec.split(",")))
	file.close()
  	
  	# try to write the swiss-prot ec numbers in a file
  	ec_list = [x.split("\t")[1].strip(" \n\t\r") for x in open(sp_file, "r") if "NA" not in x]
  	infile_extension = sp_file.split(".")[-1]
	outfile = sp_file[:len(sp_file) - len(infile_extension) - 1] + ".ec."+ level +".list"
	
	if 4 - int(level) == 0:
		ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ec_list \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	
	if not os.path.exists(outfile):
	  	file = open(outfile, 'w')
		for ec in set(ec_list): file.write("%s\n" % ec)
	  	file.close()

def write_snp(predictibles, prog_pred_dic, ref_dic, level, outfile):
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
						
	out = open(outfile, "w")
	out.write("{0} {1} {2} {3}\n".format("EC".ljust(12), "\tsens".ljust(12),\
	"\tprec".ljust(12), "\tf1score".ljust(12)))					
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
		out.write("{0} {1} {2} {3}\n".format(str(ec).ljust(12), "\t" + str(sens).ljust(12),\
		"\t" + str(prec).ljust(12), "\t" + str(f1score).ljust(12)))
	out.close()
		
if __name__ == '__main__':
	import os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--prog", dest = "prog", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted prog file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--path", dest = "path",
						help = "Path to the limputils directory")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("--out",
		                dest = "output",
		                type = str,
		                default = "[Input filename without its extension].snp.level.txt",
		                help = "The output file name")
	args = parser.parse_args()
	if args.prog is None or args.sp is None or args.level is None or args.path is None: 
		parser.print_help()
		exit(0)
	#
	#
	ref_dic = get_sp_data(args.sp)
	prog_pred_dic = get_prog_data(args.prog, args.path)
	#write_list_ecNumbers(args.prog, args.sp, args.path, args.level)
	predictibles_ecs = get_enzyme_ecs(args.level)
	infile_ext = args.prog.split(".")[-1]
	#outfile = args.prog[:len(args.prog) - len(infile_ext) - 1] + ".snp."+ str(args.level) +".txt"
	default_output = "[Input filename without its extension].snp.level.txt"
	if args.output == default_output: 
		args.output = get_output_name(args.prog)+ ".snp."+ str(args.level) + ".txt"
	write_snp(predictibles_ecs, prog_pred_dic, ref_dic, args.level, args.output)
