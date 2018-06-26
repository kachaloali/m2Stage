#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
def comp_prog_and_sp(prog_pred_dic, ref_dic, level, output):
	"""
	This function produces a comparison file reporting a binary matrix of ec numbers found or not.
		
	TP (True-Positive): Enzymatic activities predicted by the program reported in SWISS-PROT. 
	FP (False-Positive): Enzymatic activities predicted by the program not reported in SWISS-PROT.
	FN (False-Negative): Enzymatic activities reported in SWISS-PROT missed by the program.
	TN (True-Negative): Enzymatic activities not reported in SWISS-PROT and also missed by the program.
	"""
	out = open(output, "w")
	out.write("{0} {1} {2}\n".format("ec".ljust(12), "\tprog".ljust(12), "\tsp".ljust(12)))
	for protID in ref_dic:
		if protID in prog_pred_dic and prog_pred_dic[protID] != ["NA"] and ref_dic[protID] != ["NA"]:
			list_ECs_prog = prog_pred_dic[protID]
			list_ECs_uniprot = ref_dic[protID]
			
			if (4 - int(level) > 0):
				list_ECs_prog = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog = list(set(list_ECs_prog))     # to avoid redondancies
				list_ECs_uniprot = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot = list(set(list_ECs_uniprot)) # to avoid redondancies
			elif (4 - int(level) == 0):
				list_ECs_prog = [ec for ec in list_ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
				list_ECs_uniprot = [ec for ec in list_ECs_uniprot 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			
			for ec in list_ECs_prog:
				if ec in list_ECs_uniprot:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), 
					"\t" + str(1).ljust(12)))
				else:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), \
					"\t" + str(0).ljust(12)))
					if ec == "2.1.1.240":
						print protID, list_ECs_prog, list_ECs_uniprot
						#exit(0)
					#print protID, list_ECs_prog, list_ECs_uniprot  
			for ec in list_ECs_uniprot:
				if ec not in list_ECs_prog:
					out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(0).ljust(12), \
					"\t" + str(1).ljust(12)))
		elif protID in prog_pred_dic and prog_pred_dic[protID] != ["NA"] and ref_dic[protID]==["NA"]:
			list_ECs_prog = prog_pred_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_prog = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_prog \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_prog = list(set(list_ECs_prog))     # to avoid redondancies
			else: 
				list_ECs_prog = [ec for ec in list_ECs_prog \
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			for ec in list_ECs_prog:
				out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(1).ljust(12), \
				"\t" + str(0).ljust(12)))
		elif protID in prog_pred_dic and prog_pred_dic[protID] == ["NA"] and ref_dic[protID]!=["NA"]:
			list_ECs_uniprot = ref_dic[protID]
			if (4 - int(level) > 0):
				list_ECs_uniprot = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in list_ECs_uniprot \
				if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
				list_ECs_uniprot = list(set(list_ECs_uniprot)) # to avoid redondancies
			else: 
				list_ECs_uniprot = [ec for ec in list_ECs_uniprot 
				if len([x for x in ec.split(".") if x != "-"]) == int(level)]
			for ec in list_ECs_uniprot:
				out.write("{0} {1} {2}\n".format(ec.ljust(12), "\t" + str(0).ljust(12), \
				"\t" + str(1).ljust(12)))
		elif protID in prog_pred_dic and prog_pred_dic[protID] == ["NA"] and ref_dic[protID] ==["NA"]:
			pass
	out.close()
	
if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)          	    
	parser.add_argument("--prog", dest = "prog", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The converted program file")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted Swiss-Prot file")
	parser.add_argument("--plud", dest = "plud",
						help = "Path to the limputils directory")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("--out",
		                dest = "output",
		                type = str,
		                default = "[First input filename without its extension].comp",
		                help = "The output file name")
	args = parser.parse_args()
	if args.prog is None or args.sp is None or args.output is None or args.plud is None: 
		parser.print_help()
		exit(0)
	default_output = "[First input filename without its extension].comp"
	if args.output == default_output: args.output = get_output_name(args.prog)+ "."+ args.level + ".comp"
	prog_pred_dic = get_prog_data(args.prog, args.plud)
	ref_dic = get_sp_data(args.sp)
	comp_prog_and_sp(prog_pred_dic, ref_dic, args.level, args.output)
