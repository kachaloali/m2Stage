#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->

def get_sp_data(fasta, sp_file):
	all_protIDs = [line.split("|")[1] for line in open(fasta, "r") if line.startswith(">")]
	all_protIDs = list(set(all_protIDs))
	uniprot_ref_file = open(sp_file, "r")
	ref_dic = {}
	for line in uniprot_ref_file:
		protID = line.split()[0].strip(" \n\t\r")
		ec = line.split()[1].strip(" \n\t\r")
		if protID in all_protIDs:
			if protID in ref_dic: ref_dic[protID].append(ec)
			else: ref_dic[protID] = [ec]
	uniprot_ref_file.close()
	return ref_dic

def get_kaas_data(kaas_file):
	kaas_pred_file = open(kaas_file, "r")
	kaas_pred_dic = {}
	for line in kaas_pred_file:
		protID = line.split()[0].strip(" \n\t\r")
		if "|" in protID: protID = protID.split("|")[1]
		ec = line.split()[1].strip(" \n\t\r")
		if protID in kaas_pred_dic: kaas_pred_dic[protID].append(ec)
		else: kaas_pred_dic[protID] = [ec]
	kaas_pred_file.close()
	return kaas_pred_dic

def get_ref_ecNumbers(ref_dic):
	return [ x for x in set(reduce(lambda l1, l2: l1 + l2, ref_dic.values())) if x != "NA"]
	
def get_kaas_ecNumbers(kaas_pred_dic):
	return [ x for x in set(reduce(lambda l1, l2: l1 + l2, kaas_pred_dic.values())) if x != "NA"]

def write_list_ecNumbers(kaas_list_ecNumbers, ref_list_ecNumbers):
	file = open("kaas_ecNumbers.txt", 'w')
	for ec in kaas_list_ecNumbers:
	  	file.write("%s\n" % ec)
	file.close()
  	
  	file = open("sp_ecNumbers.txt", 'w')
	for ec in ref_list_ecNumbers:
  		file.write("%s\n" % ec)
  	file.close()

def get_Snp_for_ecNumbers(kaas_pred_dic, ref_dic, level):
	dic_ecNumbers = {}
	for protID in ref_dic:
		if protID in kaas_pred_dic and kaas_pred_dic[protID] != ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_kaas = kaas_pred_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_kaas_class = ['.'.join(ec[0].split('.')[:-4 +int(level)]) for ec in list_ECs_kaas]
				list_ECs_kaas_class = list(set(list_ECs_kaas_class))     # to avoid redondancies
				list_ECs_uniprot_class = ['.'.join(ec.split('.')[:-4+int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_class = list(set(list_ECs_uniprot_class)) # to avoid redondancies
			elif (4 - int(level) == 0):
				list_ECs_kaas_class = list_ECs_kaas
				list_ECs_uniprot_class = list_ECs_uniprot
			
			for ec in list_ECs_kaas_class:
				if ec in list_ECs_uniprot_class:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][0]+= 1
					else: dic_ecNumbers[ec] = [1,0,0]
				else:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][1]+= 1
					else: dic_ecNumbers[ec] = [0,1,0]			
			
			for ec in list_ECs_uniprot_class:
				if ec not in list_ECs_kaas_class:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][2]+= 1
					else: dic_ecNumbers[ec] = [0,0,1]		
		
		elif protID in kaas_pred_dic and kaas_pred_dic[protID] != ['NA'] and ref_dic[protID] == ['NA']:
			list_ECs_kaas = kaas_pred_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_kaas_class = ['.'.join(ec[0].split('.')[:-4+int(level)]) for ec in list_ECs_kaas]
				list_ECs_kaas_class = list(set(list_ECs_kaas_class))     # to avoid redondancies
			else: list_ECs_kaas_class = list_ECs_kaas
			for ec in list_ECs_kaas_class:
				if ec in dic_ecNumbers: dic_ecNumbers[ec][1]+= 1
				else: dic_ecNumbers[ec] = [0,1,0]		
		elif protID in kaas_pred_dic and kaas_pred_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot_class = ['.'.join(ec.split('.')[:-4+int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_class = list(set(list_ECs_uniprot_class)) # to avoid redondancies
			else: list_ECs_uniprot_class = list_ECs_uniprot
			for ec in list_ECs_uniprot_class:
				if ec in dic_ecNumbers: dic_ecNumbers[ec][2]+= 1
				else: dic_ecNumbers[ec] = [0,0,1]		
		elif protID in kaas_pred_dic and kaas_pred_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']: continue
		elif protID not in kaas_pred_dic and ref_dic[protID] != ["NA"]:
			list_ECs_uniprot = ref_dic[protID]
			for ec in list_ECs_uniprot:
				if ec not in list_ECs_kaas:
					if ec in dic_ecNumbers: dic_ecNumbers[ec][2]+= 1
					else: dic_ecNumbers[ec] = [0,0,1]
		
		elif protID not in kaas_pred_dic and ref_dic[protID] == ["NA"]: continue
	
	
	file = open("kaas_snp.txt", "w")
	file.write("{0} {1} {2}\n".format("EC".ljust(12), "\tsensitivity".ljust(6), "\tprecision".ljust(6)))	
	for ec, l in dic_ecNumbers.items():
		# compute sensitivity: TP / (TP + FN) 
		snp = []
		num = l[0]
		den = l[0] + l[2]
		if den != 0: op = round(float(num)/den, 2)
		else: op = float(0)
		snp.append(op)
		
		# compute precision: TP / (TP + FP)
		num = l[0]
		den = l[0] + l[1]
		if den != 0: op = round(float(num)/den, 2)
		else: op = float(0)
		snp.append(op)
		file.write("{0} {1} {2}\n".format(ec.ljust(12), ("\t" +str(snp[0])).ljust(12), ("\t" +str(snp[1])).ljust(6)))
	file.close()	
	

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
	parser.add_argument("--fasta", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")           	    
	parser.add_argument("--kaas", dest = "kaas", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted kaas file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	args = parser.parse_args()
	if args.fasta is None or args.kaas is None or args.sp is None or args.level is None: 
		parser.print_help()
		exit(0)
	#
	ref_dic = get_sp_data(args.fasta, args.sp)
	kaas_pred_dic = get_kaas_data(args.kaas)
	get_Snp_for_ecNumbers(kaas_pred_dic, ref_dic, args.level)

#	ref_ecNumbers = get_ref_ecNumbers(ref_dic)
#	kaas_ecNumbers = get_kaas_ecNumbers(kaas_pred_dic)
#	write_list_ecNumbers(kaas_ecNumbers, ref_ecNumbers)
	
