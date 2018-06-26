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

def get_priam_data(priam_file):
	priam_pred_file = open(priam_file, "r")
	priam_prediction_dic = {}
	for line in priam_pred_file:
		protID = line.split("|")[1].strip(" \n\t\r")
		ec = line.split()[-2].strip(" \n\t\r")
		if ec != 'NA':
			proba = line.split()[-1].strip(" \n\t\r")
			if protID in priam_prediction_dic: priam_prediction_dic[protID].append((ec, proba))
			else: priam_prediction_dic[protID] = [(ec, proba)]
		else: 
			priam_prediction_dic[protID] = ['NA']
	for protID in ref_dic:
		if protID not in priam_prediction_dic: priam_prediction_dic[protID] = ['NA']
	priam_pred_file.close()
	return priam_prediction_dic

def get_stats(fasta, priam_prediction_dic, ref_dic, level):
	# Stats on ECs prediction priam vs. uniprot
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	# Stats on ECs prediction priam vs. uniprot on class 1,2,3,4
	priam_uniprot = open(fasta_path +"/priam_Snp.txt", "w")
	priam_uniprot.write("{0} {1} {2} {3} {4}\n".format("protID".ljust(12), "ecNumber".ljust(12), "TP".ljust(12), 
										"FP".ljust(12),"FN".ljust(12)))
	#
	#
	nb_TP_total = 0
	nb_FP_total = 0
	nb_FN_total = 0
	for protID in ref_dic:
		nb_TP = 0
		nb_FP = 0
		nb_FN = 0
		if priam_prediction_dic[protID] != ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_priam = priam_prediction_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_priam_class = ['.'.join(ec[0].split('.')[:-4 +int(level)]) for ec in list_ECs_priam]
				list_ECs_priam_class = list(set(list_ECs_priam_class))     # to avoid redondancies
				list_ECs_uniprot_class = ['.'.join(ec.split('.')[:-4+int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_class = list(set(list_ECs_uniprot_class)) # to avoid redondancies
			elif (4 - int(level) == 0):
				list_ECs_priam_class = [ec for (ec, score) in list_ECs_priam]
				list_ECs_uniprot_class = list_ECs_uniprot
			
			for ec in list_ECs_priam_class:
				if ec in list_ECs_uniprot_class:
					priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
										str(1).ljust(12), str(0).ljust(12), str(0).ljust(12))) 
					nb_TP+= 1
				else:
					priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
										str(0).ljust(12), str(1).ljust(12), str(0).ljust(12)))  
					nb_FP+= 1
			
			for ec in list_ECs_uniprot_class:
				if ec not in list_ECs_priam_class:
					priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
										str(0).ljust(12), str(0).ljust(12), str(1).ljust(12))) 
					nb_FN+= 1
		
		elif priam_prediction_dic[protID] != ['NA'] and ref_dic[protID] == ['NA']:
			list_ECs_priam = priam_prediction_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_priam_class = ['.'.join(ec[0].split('.')[:-4+int(level)]) for ec in list_ECs_priam]
				list_ECs_priam_class = list(set(list_ECs_priam_class))     # to avoid redondancies
			else: list_ECs_priam_class = [ec for (ec, score) in list_ECs_priam]
			for ec in list_ECs_priam_class:
				priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
									str(0).ljust(12), str(1).ljust(12), str(0).ljust(12)))  
				nb_FP+= 1
		
		elif priam_prediction_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot_class = ['.'.join(ec.split('.')[:-4+int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_class = list(set(list_ECs_uniprot_class)) # to avoid redondancies
			else: list_ECs_uniprot_class = list_ECs_uniprot
			for ec in list_ECs_uniprot_class:
				priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
										str(0).ljust(12), str(0).ljust(12), str(1).ljust(12))) 	
				nb_FN+= 1
		
		elif priam_prediction_dic[protID] == ['NA'] and ref_dic[protID] != ['NA']: continue						
		#
		#
		#priam_uniprot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), str(nb_TP).ljust(12), 
		#									str(nb_FP).ljust(12), str(nb_FN).ljust(12)))
		#
		#
		nb_TP_total+= nb_TP
		nb_FP_total+= nb_FP
		nb_FN_total+= nb_FN
	priam_uniprot.close()
	

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
	parser.add_argument("--priam", dest = "priam", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted priam file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted sp file")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	args = parser.parse_args()
	if args.fasta is None or args.priam is None or args.sp is None or args.level is None: 
		parser.print_help()
		exit(0)
	#
	ref_dic 			 = get_sp_data(args.fasta, args.sp)
	priam_prediction_dic = get_priam_data(args.priam)
	get_stats(args.fasta, priam_prediction_dic, ref_dic, args.level)
