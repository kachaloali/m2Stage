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
	prog_dic, intersect_dic, valid_ecs = dict(), dict(), dict()
	conv_files = {args.conv_files[i]: args.conv_files[i+1] for i in range(0,len(args.conv_files), 2)}
	for key, value in conv_files.items():
		data = get_prog_data(value, args.plud)
		ec_list = [data[x] for x in data if "NA" not in data[x]]
		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
		prog_dic[key] = ec_list
		
#	if args.priam is not None:
#		data = get_prog_data(args.priam, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "PRIAM"
#		prog_dic[prog_name] = ec_list
#	
#	if args.e2p2 is not None:
#		data = get_prog_data(args.e2p2, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "E2P2"
#		prog_dic[prog_name] = ec_list
#	
#	if args.b2g is not None:
#		data = get_prog_data(args.b2g, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "Blast2go"
#		prog_dic[prog_name] = ec_list
#	
#	if args.kaas is not None:
#		data = get_prog_data(args.kaas, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "KAAS"
#		prog_dic[prog_name] = ec_list
#	
#	if args.koala is not None:
#		data = get_prog_data(args.koala, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "KOALA"
#		prog_dic[prog_name] = ec_list
#	
#	if args.iprscan is not None:
#		data = get_prog_data(args.iprscan, args.plud)
#		ec_list = [data[x] for x in data if "NA" not in data[x]]
#		ec_list = list(set(reduce(lambda x, y: x + y, ec_list)))
#		prog_name = "INTERPRO"
#		prog_dic[prog_name] = ec_list
#	if len(prog_dic) < 1:
#		parser.print_help()
#		exit(0)
					
	# All possible partitions are calculated between lists of ecs numbers from prediction programs. The partitions here are 
	# exactly the same as when applying a Venn diagramm to calculate the intersections between many sets.
	#powerset1 = itertools.chain.from_iterable(itertools.combinations(prog_dic.keys(), y) for y in range(len(prog_dic) + 1))
	#powerset2 = itertools.chain.from_iterable(itertools.combinations(prog_dic.values(), y) for y in range(len(prog_dic)+1))
	#titles = [part for part in powerset1 if any(map(None, part))]
	#sources = [part for part in powerset2 if any(map(None, part))]
	#for elt in (zip(titles, sources)):
	#	key = "∩".join(str(x) for x in list(elt[0]))
	#	intersect_dic[key] = reduce(lambda x, y: set(x).intersection(y), list(elt[1])) 
	#	 
	# List of valid ec numbers are calculated according to the parameters mpi that's means: The minimum number of programs 
	# constituting an intersection to be considered.
	#for key, value in intersect_dic.items():
	#	if len(key.split("∩")) >= int(args.mpi): potentials_ecs.update(value)
	#	if len(key.split("∩")) == len(prog_dic): super_intersect_list = value
	#for key, value in intersect_dic.items():
	#	if len(key.split("∩")) == 1:
	#		extreme_ec_list = set(value) - set(super_intersect_list)
	#		potentials_ecs.update(extreme_ec_list)
	
	
	# Determine the list of ec numbers with f1-score greater than the thresold (potentials ec numbers)
	# The ec numbers predicted by the programs are first validated with the library of validation: lipmutils.
	annot_ec_file = open(args.listec, "r")
	dic_annot_ec = dict()
	for line in annot_ec_file.readlines()[1:]:
		line = line.split("\t")
		ec = line[0].strip(" \n\t\r")
		prog = line[1].strip(" \n\t\r")
		f1score = line[7].strip(" \n\t\r")
		if ec not in dic_annot_ec: 
			dic_annot_ec[ec] = dict()
		if float(f1score) > 0: 
			dic_annot_ec[ec][prog] = float(f1score)
	
	potentials_ecs, rest_ecs = dict(), dict()			
	for ec, dic in dic_annot_ec.items():
		for prog in dic:
			f1score = dic[prog]
			if float(f1score) >= float(args.score):
				if ec in potentials_ecs:
					potentials_ecs[ec][prog] = f1score
					potentials_ecs[ec]["N"]+= 1
				else:
					potentials_ecs[ec] = dict()
					potentials_ecs[ec]["N"] = 1
					potentials_ecs[ec][prog] = f1score
			else:
				if ec in rest_ecs: 
					rest_ecs[ec][prog] = f1score
					rest_ecs[ec]["N"]+= 1
				else:
					rest_ecs[ec] = dict()
					rest_ecs[ec]["N"] = 1
					rest_ecs[ec][prog] = f1score
					
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
#	# PRIAM
#	if args.priam is not None:
#		dic = get_prog_data(args.priam, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic: 
#				merged_dic[protID]["PRIAM"] = set(ec_list)
#			else: 
#				merged_dic[protID] = dict()
#				merged_dic[protID]["PRIAM"] = set(ec_list)
#	# E2P2
#	if args.e2p2 is not None:
#		dic = get_prog_data(args.e2p2, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic:
#				merged_dic[protID]["E2P2"] = set(ec_list)
#			else:
#				merged_dic[protID] = dict()
#				merged_dic[protID]["E2P2"] = set(ec_list)
#	# Blast2go
#	if args.b2g is not None:
#		dic = get_prog_data(args.b2g, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic:
#				merged_dic[protID]["Blast2go"] = set(ec_list)
#			else:
#				merged_dic[protID] = dict()
#				merged_dic[protID]["Blast2go"] = set(ec_list)
#	# KAAS
#	if args.kaas is not None:
#		dic = get_prog_data(args.kaas, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic:
#				merged_dic[protID]["KAAS"] = set(ec_list)
#			else:
#				merged_dic[protID] = dict()
#				merged_dic[protID]["KAAS"] = set(ec_list)
#	# KOALA
#	if args.koala is not None:
#		dic = get_prog_data(args.koala, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic:
#				merged_dic[protID]["KOALA"] = set(ec_list)
#			else:
#				merged_dic[protID] = dict()
#				merged_dic[protID]["KOALA"] = set(ec_list)
#	# INTERPRO
#	if args.iprscan is not None:
#		dic = get_prog_data(args.iprscan, args.plud)
#		for protID, ec_list in dic.items():
#			if 4 - int(args.level) == 0:
#				ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
#			else:
#				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
#				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
#			if protID in merged_dic:
#				merged_dic[protID]["INTERPRO"] = set(ec_list)
#			else:
#				merged_dic[protID] = dict()
#				merged_dic[protID]["INTERPRO"] = set(ec_list)

	# The last step is that of the annotation, it consists in writing the result in an output file to be visualized.
	# To assign an ec number to a given protein, the fraction between the number of programs that predicted this ec 
	# number and the number of programs wich are able of predicting it must be greater than or equal to the fraction 
	# of threshold passed as a parameter (with option: -frac)
	# The programs that have predicted an ec number and are not in the list of programs capable of detecting the this 
	# ec number, and also they are not programs eliminated because of the insufficiency of their F1score, then their 
	# prediction is retained. on the other hand for this precise case, the annotation is followed by a flag 
	# indicating that the assignment is putative (the flag is: uncertified). The annotation having for flag certified,
	# means that all the conditions of the different thresholds are reached.
	nb_ec_new_good_affectation = 0
	annot = open(args.out, "w")
	annot.write("{0} {1} {2} {3}\n".format("ProtIDs".ljust(12),"\tECs".ljust(12),"\tScores".ljust(12),"\tFlags"))
	if args.oroc is not None: roc_outpout = open(str(args.oroc), "w")		
	for protID, dic in merged_dic.items():
		if len(dic) != 0:
			list_news_progs = dic.keys()
			certified_ecs = list()
			uncertified_ecs = dict()
			dic_ec_score = dict()
			raw_ecs = dict()
			for prog, list_ecs in dic.items():
				for ec in list_ecs:
					characterized_progs = [prg for prg in list_news_progs if ec in dic[prg]]
					if ec in potentials_ecs:
						corresponded_progs = set(characterized_progs).intersection(potentials_ecs[ec].keys())
					else:
						corresponded_progs = set()
					uncorresponded_progs = set(characterized_progs).difference(corresponded_progs)
					if ec in rest_ecs:
						rest_associated_progs = set(rest_ecs[ec].keys())
					else:
						rest_associated_progs = set()
					uncorresponded_progs = uncorresponded_progs.difference(rest_associated_progs)			 			
					if ec in potentials_ecs and prog in potentials_ecs[ec].keys():
						N = potentials_ecs[ec]["N"]
						new_N = len(corresponded_progs)
						if N > 0: frac = float(new_N) / N
						else: frac = 0
						############################################################################################
						f1_score = potentials_ecs[ec][prog]			 #  To build the roc curbe we needs all the    #
						if ec in raw_ecs: 							 #  information without any filtering on the   #
							raw_ecs[ec]["f1_score"].append(f1_score) #  different thresholds. This will make it	   #
				 		else:										 #  possible to follow the effects of the	   #
				 			raw_ecs[ec] = dict()					 #  computed scores and to decide on the most  #
				 			raw_ecs[ec]["f1_score"] = [f1_score]	 #  appropriate parameters.	that's why we save #
				 			raw_ecs[ec]["frac"] = frac				 #	it in case to make roc curbe of prediction #
						############################################################################################
						if frac >= float(args.frac):
						 	if prog in corresponded_progs:	
						 		certified_ecs.append(ec)
						 		f1_score = potentials_ecs[ec][prog]
						 		if ec in dic_ec_score: dic_ec_score[ec]["f1_score"].append(f1_score)
						 		else:
						 			dic_ec_score[ec] = dict()
						 			dic_ec_score[ec]["f1_score"] = [f1_score]
						 			dic_ec_score[ec]["frac"] = frac
						else:
							new_N += len(uncorresponded_progs)
							if N > 0: frac = float(new_N) / len(characterized_progs)
							else: frac = 0
							if frac >= float(args.frac):
								f1_score = potentials_ecs[ec][prog]
								if ec in uncertified_ecs: 
									uncertified_ecs[ec]["f1_score"].append(f1_score)
						 		else:
						 			uncertified_ecs[ec] = dict()
						 			uncertified_ecs[ec]["f1_score"] = [f1_score]
						 			uncertified_ecs[ec]["frac"] = frac
					elif ec in rest_ecs and prog in rest_ecs.keys():
						N = rest_ecs[ec]["N"]
						new_N = len(corresponded_progs)
						if N > 0: frac = float(new_N) / N
						else: frac = 0
						############################################################################################
						f1_score = rest_ecs[ec][prog]			 	 #  To build the roc curbe we needs all the    #
						if ec in raw_ecs: 							 #  information without any filtering on the   #
							raw_ecs[ec]["f1_score"].append(f1_score) #  different thresholds. This will make it	   #
				 		else:										 #  possible to follow the effects of the	   #
				 			raw_ecs[ec] = dict()					 #  computed scores and to decide on the most  #
				 			raw_ecs[ec]["f1_score"] = [f1_score]	 #  appropriate parameters.	that's why we save #
				 			raw_ecs[ec]["frac"] = frac				 #	it in case to make roc curbe of prediction #
						############################################################################################
			for ec in list(set(certified_ecs)):
				mean_f1_score = float(sum(dic_ec_score[ec]["f1_score"])) / len(dic_ec_score[ec]["f1_score"])
				frac_ec = round(dic_ec_score[ec]["frac"], 8)
				ec_pred_score = round(2 *(mean_f1_score * frac_ec) / (mean_f1_score + frac_ec), 8)
				annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t"+ str(ec).ljust(12),
				"\t" + str(ec_pred_score).ljust(12), "\tcertified"))
				nb_ec_new_good_affectation+= 1
			
			for ec in list(set(uncertified_ecs)):
				mean_f1_score=float(sum(uncertified_ecs[ec]["f1_score"]))/len(uncertified_ecs[ec]["f1_score"])
				frac_ec = round(uncertified_ecs[ec]["frac"], 8)
				den = mean_f1_score + frac_ec
				if den > 0:	ec_pred_score = round(2 *(mean_f1_score * frac_ec) / (mean_f1_score + frac_ec), 8)
				else: ec_pred_score = 0
				annot.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12), "\t" + str(ec).ljust(12),
				"\t" + str(ec_pred_score).ljust(12), "\tuncertified"))
			
			# when the user emits the wish to build the file to draw roc curbe, the raw data saved above are written 
			# in an output file
			if args.oroc is not None:
				for ec in list(set(raw_ecs)):
					mean_f1_score = float(sum(raw_ecs[ec]["f1_score"])) / len(raw_ecs[ec]["f1_score"])
					frac_ec = raw_ecs[ec]["frac"]
					den = mean_f1_score + frac_ec
					if den > 0:	ec_pred_score = 2 *(mean_f1_score * frac_ec) / (mean_f1_score + frac_ec)
					else: ec_pred_score = 0
					roc_outpout.write("{0} {1} {2} {3}\n".format(str(protID).ljust(12),"\t"+str(ec).ljust(12),
					"\t"+str(ec_pred_score).ljust(12), "\t"))
	print "Number of news good annotations: "+ str(nb_ec_new_good_affectation)
	print "\tFINISHED"
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
#	parser.add_argument("-priam", dest = "priam",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted priam program file")
#	parser.add_argument("-e2p2", dest = "e2p2",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted e2p2 program file")
#	parser.add_argument("-b2g", dest = "b2g",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted blast2go program file")
#	parser.add_argument("-kaas", dest = "kaas",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted kaas program file")
#	parser.add_argument("-koala", dest = "koala",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted koala program file")
#	parser.add_argument("-iprscan", dest = "iprscan",
#						type = lambda arg: is_valid_file(parser, arg),
#					   	help = "The converted interproscan program file")
	parser.add_argument("-listec", dest = "listec",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The file containing the list of ec numbers of reference")
	parser.add_argument("-plud", dest = "plud", help = "Path to the limputils directory")
	parser.add_argument("-score", dest = "score", default = "0.75", help = "The f1 score value")
	parser.add_argument("-frac", dest = "frac", default = "0.75", 
						help = "The minimum number of programs repporting an ec number to be considered as FP") 
	parser.add_argument("-level", dest = "level", default = "4", help = "The level of ec classe")
	parser.add_argument("-out", dest = "out", default = "annot.out.txt", help = "The output file name")
	parser.add_argument("-oroc", dest = "oroc", default = "out_roc.txt", help = "The output for roc displaying")	
	args = parser.parse_args()
	if args.score is None or args.plud is None or args.listec is None or len(args.conv_files) < 2: 
		parser.print_help()
		exit(0)
	main(args = args)	
