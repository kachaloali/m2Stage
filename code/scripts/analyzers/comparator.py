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
	out.write("{0} {1} {2}\n".format("ec".ljust(12), "\tpred".ljust(12), "\tstd".ljust(12)))
	for protID in ref_dic:
		if protID in pred_dic and pred_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
			ECs_prog = pred_dic[protID]
			ECs_std = ref_dic[protID]
			
			if (4 - int(level) > 0):
				ECs_prog = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				ECs_prog = list(set(ECs_prog))     # to avoid redondancies
				ECs_std = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ECs_std \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				ECs_std = list(set(ECs_std))
			elif (4 - int(level) == 0):
				ECs_prog = [ec for ec in ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				ECs_std = [ec for ec in ECs_std 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			
			for ec in ECs_prog:
				if ec in ECs_std:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), 
					"\t" + str(1).ljust(12)))
				else:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), \
					"\t" + str(0).ljust(12))) 
			for ec in ECs_std:
				if ec not in ECs_prog:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(0).ljust(12), \
					"\t" + str(1).ljust(12)))
		elif protID in pred_dic and pred_dic[protID] != ["NA"] and ref_dic[protID]==["NA"]:
			ECs_prog = pred_dic[protID]
			if (4 - int(level) > 0):
				ECs_prog = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				ECs_prog = list(set(ECs_prog))
			else: 
				ECs_prog = [ec for ec in ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			for ec in ECs_prog:
				out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), \
				"\t" + str(0).ljust(12)))
		elif protID in pred_dic and pred_dic[protID] == ["NA"] and ref_dic[protID]!=["NA"]:
			ECs_std = ref_dic[protID]
			if (4 - int(level) > 0):
				ECs_std = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ECs_std \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				ECs_std = list(set(ECs_std))
			else: 
				ECs_std = [ec for ec in ECs_std 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			for ec in ECs_std:
				out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(0).ljust(12), \
				"\t" + str(1).ljust(12)))
		elif protID in pred_dic and pred_dic[protID] == ["NA"] and ref_dic[protID] ==["NA"]:
			pass
	out.close()
	
if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)          	    
	parser.add_argument("-pp", dest = "pp", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted program file")
	parser.add_argument("-std", dest = "std",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted Swiss-Prot file")
	parser.add_argument("-plud", dest = "plud",
						help = "Path to the limputils directory")
	parser.add_argument("-level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("-out",
		                dest = "output",
		                type = str,
		                default = "[First input filename without its extension].comp",
		                help = "The output file name")
	args = parser.parse_args()
	if args.pp is None or args.std is None or args.output is None or args.plud is None: 
		parser.print_help()
		exit(0)
	default_output = "[First input filename without its extension].comp"
	if args.output == default_output: args.output = get_output_name(args.pp)+ "."+ args.level + ".comp"
	pred_dic = get_prog_data(args.pp, args.plud)
	ref_dic = get_std_data(args.std)
	comp_prog_and_std(pred_dic, ref_dic, args.level, args.output)
