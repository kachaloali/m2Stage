#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
import subprocess, commands, sys, os

def runPriam(fasta, tab, threshold, min_pro_len, max_overlap):
	"""
	Run Priam and converts its output file and the standard uniprot annotation file into an unique format.
	"""
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	try:
		print commands.getoutput("python runPriam.py -f "+ str(fasta) +" -pt "+ str(threshold) +" -mp "
								+ str(min_pro_len) +" -mo "+ str(max_overlap))
	except Exception as error: print error
	
	# Convert uniprot reference file 
	command = "python convUniprot.py --fasta "+ fasta +" --tab "+ tab 
	print commands.getoutput(command)
	
	# Convert priam output file
	command = "python ../../converters_1/convPRIAM.py -f ../resPriam/RESULTS/paj_testPriam_seqsECs.txt"
	print commands.getoutput(command)
	
	# Move priam prediction in taxon files
	print commands.getoutput("cp ../resPriam/RESULTS/paj_testPriam_seqsECs.txt "+ fasta_path)
	#
	print commands.getoutput("mv ../resPriam/RESULTS/paj_testPriam_seqsECs.conv "+ fasta_path)
	fasta_extension = fasta.split(".")[-1]
	converted_priam = fasta[:len(fasta) - len(fasta_extension) - 1] + "_pri.conv"
	print commands.getoutput("mv "+ fasta_path + "/paj_testPriam_seqsECs.conv "+  converted_priam)
	
def load_and_get_data(fasta, standard_sp_file, priam_pred_file):
	"""
	By taking both priam file and the standard uniprot one, this function returns 2 dictionnaries
	"""
	all_protIDs = subprocess.check_output("grep '^>sp' "+ fasta +" | cut -d'|' -f2", shell = True)
	all_protIDs = set(filter(None, all_protIDs.splitlines()))
	
	uniprot_ref_file = open(standard_sp_file, "r")
	ref_dic = {}
	for line in uniprot_ref_file:
		protID = line.split()[0].strip(" \n\t\r")
		ec = line.split()[1].strip(" \n\t\r")
		if protID in all_protIDs:
			if protID in ref_dic: ref_dic[protID].append(ec)
			else: ref_dic[protID] = [ec]
	uniprot_ref_file.close()
	
	priam_pred_file = open(priam_pred_file, "r")
	priam_prediction_dic = {}
	for line in priam_pred_file:
		protID = line.split("|")[1].strip(" \n\t\r")
		ec = line.split()[-2].strip(" \n\t\r")
		if ec != "NA":
			proba = line.split()[-1].strip(" \n\t\r")
			if protID in priam_prediction_dic: priam_prediction_dic[protID].append((ec, proba))
			else: priam_prediction_dic[protID] = [(ec, proba)]
		else: priam_prediction_dic[protID] = ["NA"]
	for protID in ref_dic:
		if protID not in priam_prediction_dic: priam_prediction_dic[protID] = ["NA"]
	priam_pred_file.close()
	return priam_prediction_dic, ref_dic

def comp_and_write(fasta, priam_prediction_dic, ref_dic):
	"""
	This function produces a comparison file and 
	another one that contains the values ​​of specificity, sensitivity and precision
	"""
	"""
	TP (vrai positif): Activités enzymatiques prédites par PRIAM qui ont été rapportées dans SWISS-PROT. 
	FP (faux positif): Activités enzymatiques prédites par PRIAM qui n'ont pas été rapportées dans SP.
	FN (faux négatif): Activités enzymatiques rapportées dans SWISS-PROT qui ont été manquées par PRIAM.
	TN (vrai négatif): Activités enzymatiques non rapportées dans SWISS-PROT et aussi manquées par PRIAM.
	"""
	TP = 0
	FP = 0
	FN = 0
	TN = 0
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	ssp = open(fasta_path +"/resCompAnnot.ssp", "w")
	ssp.write("{0} {1} {2} {3}\n".format("protID".ljust(12), "\tPriamPred".ljust(12), 
										"\tspRef".ljust(12),"\tScores".ljust(12)))
	for protID in ref_dic:
		if priam_prediction_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_priam = priam_prediction_dic[protID]
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
					
		elif priam_prediction_dic[protID] != ["NA"] and ref_dic[protID] == ["NA"]:
			list_ECs_priam = priam_prediction_dic[protID]
			for ec in list_ECs_priam:
				ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(1).ljust(12), 
					"\t"+str(0).ljust(12), "\t"+str(ec[1]).ljust(12)))
				FP+= 1
		elif priam_prediction_dic[protID] == ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_uniprot = ref_dic[protID]
			for ec in list_ECs_uniprot:
				ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(0).ljust(12), 
											"\t"+str(1).ljust(12), "\t"+str(0).ljust(12)))
				FN+= 1
		
		elif priam_prediction_dic[protID] == ["NA"] and ref_dic[protID] == ["NA"]:
			ssp.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+str(0).ljust(12), 
										"\t"+str(0).ljust(12), "\t"+str(0).ljust(12)))
			TN+= 1
	"""
	spécificité = TN / (TN + FP)
	sensibilité = TP / (TP + FN)
	précision   = TP / (TP + FP)
	"""
	ss = open(fasta_path +"/specificity.ss", "w")
	specificity = float(TN) / float(TN + FP)
	sensitivity = float(TP) / float(TP + FN)
	precision   = float(TP) / float(TP + FP)
	ss.write("{0} {1} {2}\n".format("specificity".ljust(20),"sensitivity".ljust(20),"precision".ljust(20)))
	ss.write("{0} {1} {2}\n".format(str(specificity).ljust(20), str(sensitivity).ljust(20), 
								str(precision).ljust(20)))
	ss.close()
	ssp.close()


def get_stats(fasta, priam_prediction_dic, ref_dic, level):
	# Stats on ECs prediction priam vs. uniprot
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	# Stats on ECs prediction priam vs. uniprot on class 1,2,3,4
	priam_uniprot = open(fasta_path +"/priam_vs_uniprot.txt", "w")
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
		if priam_prediction_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
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
		
		elif priam_prediction_dic[protID] != ["NA"] and ref_dic[protID] == ["NA"]:
			list_ECs_priam = priam_prediction_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_priam_class = ['.'.join(ec[0].split('.')[:-4+int(level)]) for ec in list_ECs_priam]
				list_ECs_priam_class = list(set(list_ECs_priam_class))     # to avoid redondancies
			else: list_ECs_priam_class = [ec for (ec, score) in list_ECs_priam]
			for ec in list_ECs_priam_class:
				priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
									str(0).ljust(12), str(1).ljust(12), str(0).ljust(12)))  
				nb_FP+= 1
		
		elif priam_prediction_dic[protID] == ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot_class = ['.'.join(ec.split('.')[:-4+int(level)]) for ec in list_ECs_uniprot]
				list_ECs_uniprot_class = list(set(list_ECs_uniprot_class)) # to avoid redondancies
			else: list_ECs_uniprot_class = list_ECs_uniprot
			for ec in list_ECs_uniprot_class:
				priam_uniprot.write("{0} {1} {2} {3} {4}\n".format(str(protID).ljust(12), ec.ljust(12),
										str(0).ljust(12), str(0).ljust(12), str(1).ljust(12))) 	
				nb_FN+= 1
		
		elif priam_prediction_dic[protID] == ["NA"] and ref_dic[protID] != ["NA"]: continue						
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
	
	print "{0} {1} {2} {3}\n".format(("On level <=> "+ str(level)).ljust(12), 
									("TP: " + str(nb_TP_total)).ljust(12), 
									("FP: " + str(nb_FP_total)).ljust(12),
									("FN: " + str(nb_FN_total)).ljust(12))
	
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
	parser.add_argument("--fas", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")           	    
	parser.add_argument("--tab", dest = "tab", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input Swiss-Prot tab file")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	args = parser.parse_args()
	if args.fasta is None or args.tab is None or args.level is None: 
		parser.print_help()
		exit(0)
	elif args.level == "0": 
		print("Oops! Level was no valid number. Try again with level in [1..4]\n")
		exit(0)
	
	# Best params values for priam prediction: threshold = 0.6, min_pro_len = 20, max_overlap = 40  	
	threshold = 0.6
	min_pro_len = 20
	max_overlap = 40
	fasta, tab, level = args.fasta, args.tab, args.level
	
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	fasta_extension = args.fasta.split(".")[-1]
	converted_sp = args.fasta[:len(args.fasta) - len(fasta_extension) - 1] + "_std.conv"
	converted_priam = args.fasta[:len(args.fasta) - len(fasta_extension) - 1] + "_pri.conv"

	#runPriam(fasta, tab, threshold, min_pro_len, max_overlap)
	priam_prediction_dic, ref_dic = load_and_get_data(fasta, converted_sp, converted_priam)
	comp_and_write(fasta, priam_prediction_dic, ref_dic)
	get_stats(fasta, priam_prediction_dic, ref_dic, level)






