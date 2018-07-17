#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
def comp_prog_and_std(pred_dic, ref_dic, level, output):
	"""
	This function produces a comparison file reporting a binary matrix of ec numbers found or not.
		
	TP (True-Positive): Enzymatic activities predicted by the program reported in SWISS-PROT. 
	FP (False-Positive): Enzymatic activities predicted by the program not reported in SWISS-PROT.
	FN (False-Negative): Enzymatic activities reported in SWISS-PROT missed by the program.
	TN (True-Negative): Enzymatic activities not reported in SWISS-PROT and also missed by the program.
	"""
	out = open(output, "w")
	out.write("{0} {1} {2} {3}\n".format("ecs".ljust(12), "\tpred".ljust(6), "\tscore".ljust(12), "\tstd".ljust(12)))
	for protID in ref_dic:
		if protID in pred_dic and pred_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_std = ref_dic[protID]
			list_ECs_prog = list(set(pred_dic[protID]))
			ECs_prog = list(set([ec[0].strip(" \n\t\r") for ec in list_ECs_prog]))
			 
			if (4 - int(level) > 0):
				list_ECs_std = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_std \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_std = list(set(list_ECs_std))
			elif (4 - int(level) == 0):
				list_ECs_std = [ec for ec in list_ECs_std if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				list_ECs_std = list(set(list_ECs_std))
			for ec in list_ECs_prog:
				if ec[0].strip(" \n\t\r") in list_ECs_std:
					out.write("{0} {1} {2} {3}\n".format(str(ec[0]).ljust(12), "\t1".ljust(6), "\t"+str(ec[1]).ljust(12),
					"\t1".ljust(12)))
				else:
					out.write("{0} {1} {2} {3}\n".format(str(ec[0]).ljust(12), "\t1".ljust(6), "\t"+str(ec[1]).ljust(12),
					"\t0".ljust(12))) 
			for ec in list_ECs_std:
				if ec not in ECs_prog:
					out.write("{0} {1} {2} {3}\n".format(str(ec).ljust(12), "\t0".ljust(6), "\t1".ljust(12), 
					"\t1".ljust(12)))
		
		elif protID in pred_dic and pred_dic[protID] != ["NA"] and ref_dic[protID] == ["NA"]:
			for ec in list(set(pred_dic[protID])):
				out.write("{0} {1} {2} {3}\n".format(str(ec[0]).ljust(12), "\t1".ljust(6), "\t" + str(ec[1]).ljust(12), 
				"\t0".ljust(12)))
		
		elif protID in pred_dic and pred_dic[protID] == ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_std = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_std = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_std \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_std = list(set(list_ECs_std))
			else: 
				list_ECs_std = [ec for ec in list_ECs_std if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			for ec in list_ECs_std:
				out.write("{0} {1} {2} {3}\n".format(str(ec).ljust(12), "\t0".ljust(6), "\t1".ljust(12), "\t1".ljust(12)))
	
		elif protID in pred_dic and pred_dic[protID] == ["NA"] and ref_dic[protID] == ["NA"]:
			pass
	out.close()


def get_prog_dic(pFile):
	"""
	Return a dictionnary of proteins IDs and their annotations: {protIDs: list of associated ECs}
	"""
	pred_dic = dict()
	f = open(pFile, "r")
	for line in f:
		if line.split("\t")[0].strip(" \n\t\r") != "protIDs":
			protID = line.split("\t")[0].strip(" \n\t\r")
			ec = line.split("\t")[1].strip(" \n\t\r")
			proba = line.split("\t")[2].strip(" \n\t\r")
			flag = line.split("\t")[3].strip(" \n\t\r")
			if protID in pred_dic:
				#if  flag == "certified":
				pred_dic[protID].append((ec, proba))
			else:
				#if  flag == "certified":
				pred_dic[protID] = [(ec, proba)]
	for protID, ECs in pred_dic.items():
		pred_dic[protID] = list(set(ECs))
	f.close()
	return pred_dic
	
if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)          	    
	parser.add_argument("-pp", dest = "pp", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The file produce by predictor.py script")
	parser.add_argument("-std", dest = "std",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted Swiss-Prot file")
	parser.add_argument("-level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("-out",
		                dest = "output",
		                type = str,
		                default = "matrix.txt",
		                help = "The output file name")
	args = parser.parse_args()
	if args.pp is None or args.std is None: 
		parser.print_help()
		exit(0)
	pred_dic = get_prog_dic(args.pp)
	ref_dic = get_std_data(args.std)
	comp_prog_and_std(pred_dic, ref_dic, args.level, args.output)
