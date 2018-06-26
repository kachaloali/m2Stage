#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from Bio.ExPASy import Enzyme
import subprocess, pandas as pd, os

def get_enzyme_ecs(level): 
	'''
	Reads a Expasy Enzyme .dat file and writes a tab separated file where the 
	first column is EC number, the second column is the reaction description, 
	the third column is the associated uniprot ids separated by '|', 
	and the fourth column indicates whether the reactions described 
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
	#df = pd.read_table(output_name)
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

def get_kaas_ecs(level):
	if not os.path.exists(os.path.join("database", "kaas", "ko2ec.xl")):
		curl_ko2ec = os.path.join("http://www.genome.jp", "kegg", "files", "ko2ec.xl")
		subprocess.check_output("wget -cq -P database/kaas " + curl_ko2ec, shell = True)
	if not os.path.exists(os.path.join("database", "kaas", "ko2ec.xl")): 
		print ("%s\n", "Missing kaas database!")
		exit(0)
	kaas = open(os.path.join("database", "kaas", "ko2ec.xl"), "r")
	all_ECs = []
	for line in kaas:
		if len(line.split("[EC:")) > 1:
			list_ecs = line.split("[EC:")[1].split("]")[0].split()
			all_ECs+= list_ecs

	if 4 - int(level) == 0:
		all_ECs = [ec for ec in all_ECs if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		all_ECs = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in all_ECs \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	return list(set(all_ECs))

def get_blast2go_ecs(level):
	if not os.path.exists(os.path.join("database", "blast2go", "ec2go")):
		curl_ko2ec = os.path.join("http://geneontology.org", "external2go", "ec2go")
		subprocess.check_output("wget -cq -P database/blast2go " + curl_ko2ec, shell = True)
	if not os.path.exists(os.path.join("database", "blast2go", "ec2go")): 
		print ("%s\n", "Missing blast2go database!")
		exit(0)
	blast2go = open(os.path.join("database", "blast2go", "ec2go"), "r")
	all_ECs = []
	for line in blast2go:
		if "EC:" in line:
			ec = line.split(">")[0].split(":")[1].strip(" \n\t\r")
			all_ECs.append(ec)
	if 4 - int(level) == 0:
		all_ECs = [ec for ec in all_ECs if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		all_ECs = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in all_ECs \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	return list(set(all_ECs))
		
def get_sp_data(sp_file):
	uniprot_ref_file = open(sp_file, "r")
	ref_dic = {}
	for line in uniprot_ref_file:
		protID = line.split()[0].strip(" \n\t\r")
		ec = line.split()[1].strip(" \n\t\r")
		if protID in ref_dic: ref_dic[protID].append(ec)
		else: ref_dic[protID] = [ec]
	uniprot_ref_file.close()
	return ref_dic
	
def get_prog_data(prog_file):
	prog_pred_file = open(prog_file, "r")
	prog_pred_dic = {}
	for line in prog_pred_file:
		protID = line.split("\t")[0].strip(" \n\t\r")
		if "|" in protID: protID = line.split("|")[1].strip(" \n\t\r")
		ec = line.split("\t")[1].strip(" \n\t\r")
		if ec != 'NA':
			if protID in prog_pred_dic:	prog_pred_dic[protID].append(ec)
			else: prog_pred_dic[protID] = [ec]
		else: prog_pred_dic[protID] = ['NA']
	prog_pred_file.close()
	return prog_pred_dic

def write_list_ecNumbers(prog_file, sp_file):
	# try to transfer ec numbers with lipmutils program
	lipmutils_program = os.path.join("lipmutils", "bin","lipm_valid_ec_numbers.pl")
	command = [lipmutils_program, "--infile", prog_file, "--outfile", str(prog_file)+".transf"]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	# try to write the ec numbers list in a file after being transferred
	in_file = str(prog_file) + ".transf"
	ec_list = [x.split("\t")[1].strip(" \n\t\r") for x in open(in_file, "r") if "NA" not in x]
	infile_extension = prog_file.split(".")[-1]
	outfile = prog_file[:len(prog_file) - len(infile_extension) - 1] + ".ec.list"
	file = open(outfile, "w")
	for ec in set(ec_list): file.writelines(list("%s\n" % item for item in ec.split(",")))
	file.close()
  	# try to write the swiss-prot ec numbers in a file
  	ec_list = [x.split("\t")[1].strip(" \n\t\r") for x in open(sp_file, "r") if "NA" not in x]
  	infile_extension = sp_file.split(".")[-1]
	outfile = sp_file[:len(sp_file) - len(infile_extension) - 1] + ".ec.list"
  	file = open(outfile, 'w')
	for ec in set(ec_list): file.write("%s\n" % ec)
  	file.close()

def write_snp(predictibles, prog_pred_dic, ref_dic, level, outfile):
	dic_ecNumbers = {}
	if (4 - int(level) > 0):
		predictibles_ecs = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in predictibles]
	else:
		predictibles_ecs = predictibles
	for protID in ref_dic:
		if protID in prog_pred_dic and prog_pred_dic[protID] != ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_prog = prog_pred_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
			elif (4 - int(level) == 0):
				list_ECs_prog_cl = list_ECs_prog
				list_ECs_uniprot_cl = list_ECs_uniprot
		
			for ec in list_ECs_prog_cl:
				if ec in list_ECs_uniprot_cl:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][0]+= 1
					else: dic_ecNumbers[ec] = [1,0,0]
				else:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][1]+= 1
					else: dic_ecNumbers[ec] = [0,1,0]	
			for ec in list_ECs_uniprot_cl:
				if ec not in list_ECs_prog_cl and ec in predictibles_ecs:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][2]+= 1
					else: dic_ecNumbers[ec] = [0,0,1]
					
		elif protID in prog_pred_dic and prog_pred_dic[protID] != ['NA'] and ref_dic[protID] == ['NA']:
			list_ECs_prog = prog_pred_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_prog_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog]
				list_ECs_prog_cl = list(set(list_ECs_prog_cl))     # to avoid redondancies
			else: list_ECs_prog_cl = list_ECs_prog
			
			for ec in list_ECs_prog_cl:
				if ec in dic_ecNumbers: dic_ecNumbers[ec][1]+= 1
				else: dic_ecNumbers[ec] = [0,1,0]
					
		elif protID in prog_pred_dic and prog_pred_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot_cl = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_cl = list(set(list_ECs_uniprot_cl)) # to avoid redondancies
			else: list_ECs_uniprot_cl = list_ECs_uniprot

			for ec in list_ECs_uniprot_cl:
				if ec in predictibles_ecs:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][2]+= 1
					else: dic_ecNumbers[ec] = [0,0,1]	
					
		elif protID in prog_pred_dic and prog_pred_dic[protID] == ['NA'] and ref_dic[protID] == ['NA']: pass
	
	out = open(outfile, "w")
	out.write("{0} {1} {2} {3}\n".format("EC".ljust(12), "\tsens".ljust(12), 
															 "\tspec".ljust(12), "\tprec".ljust(12)))					
	for ec, l in dic_ecNumbers.items():	
		# compute sensitivity: TP / (TP + FN) 
		num = l[0]
		den = l[0] + l[2]
		if den != 0: sens = round(float(num)/den, 8)
		else: sens = float(0)
		
		# compute specificity: TN / (TN + FP)
		if ec in predictibles_ecs: 
			num = len(predictibles_ecs) - (l[0] + l[1] + l[2])
		else:
			num = len(predictibles_ecs) - (l[0] + l[1]) + l[2]
		den = num + l[1]
		if den != 0: spec = round(float(num)/den, 8)
		else: spec = float(0)
		
		# compute precision: TP / (TP + FP)
		num = l[0]
		den = l[0] + l[1]
		if den != 0: prec = round(float(num)/den, 8)
		else: prec = float(0)
		out.write("{0} {1} {2} {3}\n".format(str(ec).ljust(12), "\t"+str(sens).ljust(12), 
											"\t"+str(spec).ljust(12), "\t"+str(prec).ljust(12)))
		
	out.close()
		
if __name__ == '__main__':
	import os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	def is_valid_file(parser, arg):
		"""
		Check if arg is a valid file that already exists on the file system.

		Parameters
		----------
		parser : argparse object
		arg    : str

		Returns
		-------
		arg
		"""
		arg = os.path.abspath(arg)
		if not os.path.exists(arg):
			parser.error("The file %s does not exist!" % arg)
		else: return arg

	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--prog", dest = "prog", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted prog file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec cles")
	parser.add_argument("--kaas", dest = "kaas", action = 'store_true',
					   	help = "switch on kaas converted file")
	parser.add_argument("--b2g", dest = "b2g", action = 'store_true',
					   	help = "switch on blast2go converted file")
	args = parser.parse_args()
	if args.prog is None or args.sp is None or args.level is None: 
		parser.print_help()
		exit(0)
	#
	ref_dic = get_sp_data(args.sp)
	prog_pred_dic = get_prog_data(args.prog)
	write_list_ecNumbers(args.prog, args.sp)
	if args.kaas:
		predictibles_ecs = get_kaas_ecs(args.level)
	elif args.b2g:
		predictibles_ecs = get_blast2go_ecs(args.level)
	else:
		predictibles_ecs = get_enzyme_ecs(args.level)
	
	infile_extension = args.prog.split(".")[-1]
	outfile = args.prog[:len(args.prog) - len(infile_extension) - 1] + ".snp."+ str(args.level) +".txt"
	write_snp(predictibles_ecs, prog_pred_dic, ref_dic, args.level, outfile)
