#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->

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

def comp_prog_and_sp(fasta, prog_pred_dic, ref_dic):
	"""
	This function produces a comparison file and 
	another one that contains the values ​​of specificity, sensitivity and precision
		
	TP (vrai positif): Activités enzymatiques prédites par PRIAM qui ont été rapportées dans SWISS-PROT. 
	FP (faux positif): Activités enzymatiques prédites par PRIAM qui n'ont pas été rapportées dans SP.
	FN (faux négatif): Activités enzymatiques rapportées dans SWISS-PROT qui ont été manquées par PRIAM.
	TN (vrai négatif): Activités enzymatiques non rapportées dans SWISS-PROT et aussi manquées par PRIAM.
	"""
	TP, FP, FN, TN = 0, 0, 0, 0
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	ssp = open(os.path.join(fasta_path, "priam_vs_uniprot.txt"), "w")
	ssp.write("{0} {1} {2} {3}\n".format("protID".ljust(12), "\tPriam".ljust(12), 
										"\tSwiss-Prot".ljust(12),"\tScores".ljust(12)))
	for protID in ref_dic:
		if prog_pred_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_priam = prog_pred_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			for ec in list_ECs_priam:
				if ec[0] in list_ECs_uniprot:
					ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(1).ljust(12), 
						"\t"+str(1).ljust(12), "\t"+str(ec[1]).ljust(12)))
					TP+= 1
				else:
					ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(1).ljust(12), 
						"\t"+str(0).ljust(12), "\t"+str(ec[1]).ljust(12)))
					FP+= 1
			
			list_ECs_priam = [ecNumber for (ecNumber, score) in list_ECs_priam]	
			for ec in list_ECs_uniprot:
				if ec not in list_ECs_priam:
					ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(0).ljust(12), 
						"\t"+str(1).ljust(12), "\t"+str(0).ljust(12)))
					FN+= 1
					
		elif prog_pred_dic[protID] != ["NA"] and ref_dic[protID] == ["NA"]:
			list_ECs_priam = prog_pred_dic[protID]
			for ec in list_ECs_priam:
				ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(1).ljust(12), 
					"\t"+str(0).ljust(12), "\t"+str(ec[1]).ljust(12)))
				FP+= 1
		elif prog_pred_dic[protID] == ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_uniprot = ref_dic[protID]
			for ec in list_ECs_uniprot:
				ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(0).ljust(12), 
											"\t"+str(1).ljust(12), "\t"+str(0).ljust(12)))
				FN+= 1
		
		elif prog_pred_dic[protID] == ["NA"] and ref_dic[protID] == ["NA"]:
			ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(0).ljust(12), 
										"\t"+str(0).ljust(12), "\t"+str(0).ljust(12)))
			TN+= 1
	"""
	spécificité = TN / (TN + FP)
	sensibilité = TP / (TP + FN)
	précision   = TP / (TP + FP)
	"""
	ss = open(os.path.join(fasta_path, "priam_specificity.txt"), "w")
	specificity = float(TN) / float(TN + FP)
	sensitivity = float(TP) / float(TP + FN)
	precision   = float(TP) / float(TP + FP)
	ss.write("{0} {1} {2}\n".format("specificity".ljust(20),"sensitivity".ljust(20),"precision".ljust(20)))
	ss.write("{0} {1} {2}\n".format(str(specificity).ljust(20), str(sensitivity).ljust(20), 
								str(precision).ljust(20)))
	ss.close()
	ssp.close()
	
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
	args = parser.parse_args()
	if args.fasta is None or args.priam is None or args.sp is None: 
		parser.print_help()
		exit(0)
	#
	ref_dic 			 = get_sp_data(args.fasta, args.sp)
	prog_pred_dic = get_priam_data(args.priam)
	comp_priam_and_sp(args.fasta, prog_pred_dic, ref_dic)
