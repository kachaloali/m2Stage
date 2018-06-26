#!/usr/bin/python2
# <!- -*- coding: utf-8 -*- ->
def main(args):
	"""
	Perform annotation from multiple results from multiple annotation tools. The principle is as follows:
	1. From the score threshold of f1-score received as parameter, we build a list of potentials ec numbers.
	i.e the list of ec having a f1-score greater than or equal to the threshold.

	2. The potentials ec numbers predicted by all programs are retained as valid ecs.
	3. An inventory of proteins is made, associating them with all the ecs assigned to it by the programs.
		For each protein and for each ec that has been associated with it:
			. If the ec exists in the list of valid ecs, then this protein is definitively associated with 
			this ec number.
			. Otherwise, the ec number will be removed.
	"""
	conv_files = {args.conv_files[i]: args.conv_files[i+1] for i in range(0,len(args.conv_files), 2)}
	if args.listec is not None:
#		prog_dic, intersect_dic, valid_ecs = dict(), dict(), dict()
#		for f in conv_files.values(): is_valid_file(parser, f)
#		for key, value in conv_files.items():
#			data = get_prog_data(value, args.plud)
#			ec_list = [data[x] for x in data if "NA" not in data[x]]
#			ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#			prog_dic[key] = ec_list
#						
		# Determine the list of ec numbers with f1-score greater than the thresold (potentials ec numbers)
		# The ec numbers predicted by the programs are first validated with the library of validation: lipmutils.
		annot_ec_file = open(args.listec, "r")
		dic_annot_ec = dict()
		for line in annot_ec_file.readlines()[1:]:
			line = line.split("\t")
			ec = line[0].strip(" \n\t\r")
			prog = line[1].strip(" \n\t\r")
			f1score = float(line[7].strip(" \n\t\r"))
			if ec not in dic_annot_ec: 
				dic_annot_ec[ec] = dict()
			if f1score > 0: 
				dic_annot_ec[ec][prog] = f1score
	
		potentials_ecs, rest_ecs = dict(), []			
		for ec, dic in dic_annot_ec.items():
			for prog in dic:
				f1_score = float(dic[prog])
				if f1_score >= float(args.score):
					if ec in potentials_ecs:
						potentials_ecs[ec][prog] = f1_score
						potentials_ecs[ec]["N"]+= 1
					else:
						potentials_ecs[ec] = dict()
						potentials_ecs[ec]["N"] = 1
						potentials_ecs[ec][prog] = f1_score
				else:
					if ec not in rest_ecs: rest_ecs.append(ec)
		# An inventory of proteins ID is made, associating them with all the ecs assigned to it by the programs 
		# of predictions.
		# The rule is as follow:
		# For each protein and for each ec that has been associated with it:
		#	. If the ec exists in the list of potentials ecs, then this protein is definitively associated with this ec.
		#	. Otherwise, the ec number will be removed.
		merged_dic = dict()
		for key, value in conv_files.items():
			dic = get_prog_data(value, args.plud)
			for protID, ec_list in dic.items():
				if 4 - int(args.level) == 0:
					ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
				else:
					ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
					if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
				if protID in merged_dic: 
					merged_dic[protID][key] = set(ec_list)
				else: 
					merged_dic[protID] = dict()
					merged_dic[protID][key] = set(ec_list)

		# The last step is that of the annotation, it consists in writing the result in an output file to be visualized.
		# To assign an ec number to a given protein, the fraction between the number of programs that predicted this ec 
		# number and the number of programs wich are able of predicting it must be greater than or equal to the fraction 
		# of threshold passed as a parameter (with option: -frac)
		# The programs that have predicted an ec number and are not in the list of programs capable of detecting the this 
		# ec number, and also they are not programs eliminated because of the insufficiency of their F1score, then their 
		# prediction is retained. on the other hand for this precise case, the annotation is followed by a flag 
		# indicating that the assignment is putative (the flag is: uncertified). The annotation having for flag certified,
		# means that all the conditions of the different thresholds are reached.
		annot = open(args.out, "w")
		annot.write("{0} {1} {2} {3}\n".format("ProtIDs".ljust(12),"\tECs".ljust(12),"\tScores".ljust(12),"\tFlags"))
		if args.oroc is not None: roc_outpout = open(str(args.oroc), "w")
#		for protID, dic in merged_dic.items():
#			if len(dic) != 0:
#				list_ecs = list(set(reduce(lambda x, y: list(x) + list(y), dic.values())))
#				certified_ecs = list()
#				uncertified_ecs = dict()
#				for ec in list_ecs:
#					N_before = len([prg for prg in merged_dic[protID] if ec in merged_dic[protID][prg]])
#					if ec in potentials_ecs: 
#						Ns = potentials_ecs[ec]["N"]
#						N_after = len([prg for prg in merged_dic[protID] if prg in potentials_ecs[ec]])
#						frac = float(N_after) / Ns
#						if frac >= float(args.frac):
#							ec_pred_score = float(N_before + N_after) / (len(conv_files) + Ns)
#							annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
#							"\t" + str(ec_pred_score).ljust(12), "\tT"))
#					else:
#						if ec not in rest_ecs:
#							N0 = len([prg for prg in merged_dic[protID] if ec in merged_dic[protID][prg]])
#							frac = float(N0) / len(conv_files)
#							if frac >= float(args.frac):
#								ec_pred_score = float(N0) / len(conv_files)
#								annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
#								"\t" + str(ec_pred_score).ljust(12), "\tT"))
#	 				
#	 				# when the user emits the wish to build the file to draw roc curbe, the raw data saved above are  
#					# written in an output file
#					if args.oroc is not None:
#						ec_pred_score = float(N_before + N_after) / (len(conv_files) + Ns)
#						roc_outpout.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12),"\t"+str(ec).ljust(12),
#						"\t"+str(ec_pred_score).ljust(12), "\t"))
#		annot.close()
#		time.sleep(0.5)
		for protID, dic in merged_dic.items():
			if len(dic) != 0:
				list_news_progs = dic.keys()
				certified_ecs = list()
				uncertified_ecs = dict()
				for prog, list_ecs in dic.items():
					for ec in list_ecs:
						N_before = len([prg for prg in merged_dic[protID] if ec in merged_dic[protID][prg]])
						characterized_progs = [prg for prg in list_news_progs if ec in dic[prg]]
						if ec in potentials_ecs:
							corresponded_progs = set(characterized_progs).intersection(potentials_ecs[ec].keys())
						else:
							corresponded_progs = set()
						uncorresponded_progs = set(characterized_progs).difference(corresponded_progs)
						rest_associated_progs = [prg for prg in list_news_progs if ec in potentials_ecs and prg not in potentials_ecs[ec].keys()]
						uncorresponded_progs = uncorresponded_progs.difference(rest_associated_progs)			 			
						if ec in potentials_ecs and prog in potentials_ecs[ec].keys():
							Ns = potentials_ecs[ec]["N"]
							N_after = len(corresponded_progs)
							if Ns > 0: frac = float(N_after) / Ns
							else: frac = 0
							if frac >= float(args.frac):
							 	if prog in corresponded_progs:	
							 		certified_ecs.append(ec)
							 		f1_score = potentials_ecs[ec][prog]
							else:
								N_after += len(uncorresponded_progs)
								if Ns > 0: frac = float(N_after) / len(characterized_progs)
								else: frac = 0
								if frac >= float(args.frac):
									f1_score = potentials_ecs[ec][prog]
									if ec not in uncertified_ecs:
										uncertified_ecs[ec] = dict()
										N0 = len(characterized_progs)
										uncertified_ecs[ec]["N0"] = N0
						
				for ec in list(set(certified_ecs)):
					ec_pred_score = float(N_before + N_after) / (len(conv_files) + Ns)
					annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
					"\t" + str(ec_pred_score).ljust(12), "\tT"))
				for ec in list(set(uncertified_ecs)):
					ec_pred_score = float(uncertified_ecs[ec]["N0"]) / len(conv_files)
					annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
					"\t" + str(ec_pred_score).ljust(12), "\tT"))
	else:
		# An inventory of proteins ID is made, associating them with all the ecs assigned to it by the programs 
		# of predictions.
		# The rule is as follow:
		# For each protein and for each ec that has been associated with it:
		#	. If the ec exists in the list of potentials ecs, then this protein is definitively associated with this ec.
		#	. Otherwise, the ec number will be removed.
		merged_dic = dict()
		for key, value in conv_files.items():
			dic = get_prog_data(value, args.plud)
			for protID, ec_list in dic.items():
				if 4 - int(args.level) == 0:
					ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
				else:
					ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
					if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
				if protID in merged_dic: 
					merged_dic[protID][key] = set(ec_list)
				else: 
					merged_dic[protID] = dict()
					merged_dic[protID][key] = set(ec_list)
		annot = open(args.out, "w")
		annot.write("{0} {1} {2} {3}\n".format("ProtIDs".ljust(12),"\tECs".ljust(12),"\tScores".ljust(12),"\tFlags"))
		if args.oroc is not None: roc_outpout = open(str(args.oroc), "w")
		for protID, dic in merged_dic.items():
			list_ecs = list(set(reduce(lambda x, y: list(x) + list(y), dic.values())))
			for ec in list_ecs:
				Np = len([prg for prg in merged_dic[protID] if ec in merged_dic[protID][prg]])
				ec_pred_score = float(Np) / len(conv_files)
				annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
				"\t" + str(ec_pred_score).ljust(12), "\tF"))
				
				# when the user emits the wish to build the file to draw roc curbe, the raw data saved above are written 
				# in an output file
				if args.oroc is not None:
					roc_outpout.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12),"\t"+str(ec).ljust(12),
					"\t"+str(ec_pred_score).ljust(12), "\t"))
		#print "Number of news good annotations: "+ str(nb_ec_new_good_affectation)
		#print "\tFINISHED"
		annot.close()
		time.sleep(0.5)

if __name__ == '__main__':
	import os, csv, itertools
	from utils import *
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("conv_files", type = str, nargs = '*',
						help = "The list of the converted programs files. Each file must be preceded by the\
						 name of the corresponded program (ex: E2P2 e2p2_file.conv PRIAM priam_file.conv)")	    
	parser.add_argument("-listec", dest = "listec",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The file containing the list of «ec numbers to pick» for the specified level")
	parser.add_argument("-plud", dest = "plud", help = "Path to the limputils directory")
	parser.add_argument("-score", dest = "score", default = "0.75", help = "The f1 score value")
	parser.add_argument("-frac", dest = "frac", default = "0.75", 
						help = "The minimum number of programs repporting an ec number to be considered as FP") 
	parser.add_argument("-level", dest = "level", default = "4", help = "The level of ec classe")
	parser.add_argument("-out", dest = "out", default = "annot.out.txt", help = "The output file name")
	parser.add_argument("-oroc", dest = "oroc", default = "rocfile.out.txt", help = "The output for roc displaying")	
	args = parser.parse_args()
	if args.score is None or args.plud is None or len(args.conv_files) % 2 != 0: 
		parser.print_help()
		exit(0)
	main(args = args)	
