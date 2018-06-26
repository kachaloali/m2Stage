# -*- coding: utf-8 -*-
"""
Created on wednesday Jan 31 12:00:47 2018

@author: Kachalo Ali
"""
from utils import *
def  convert_file(fasta, tab, output):
	"""
	By taking the fasta file and the swiss-prot .tab file as input, this function converts 
	the two files into a single output annotation file named [output].conv where output is 
	the name of output file given by the user. By default, this is the full name of the input 
	fasta file replacing simply its extention with '_std.conv' one.
	
	The output format is: entry name \t ec number. If entry name is not associated with an 
	ec number, it is replaced by NA.
	"""
	i = 0 # set i for counting progress bar status
	start_time = time.time() # set time for the beginning of converting
	all_protIDs = [line.split("|")[1] for line in open(fasta, "r") if line.startswith(">")]
	all_protIDs = list(set(all_protIDs))
	#
	#
	tab = open(tab, "r")
	next(tab)
	dic_ecs = {}
	for line in tab:
		i+= 1
		progress(i, len(all_protIDs), start_time, "Converting file, Please wait...")# progress
		protID = line.split()[0]
		if "(EC " in line:
			line = line.split("(EC ")
			if protID in all_protIDs:
				list_ECs = []
				for elt in line[1:]: list_ECs.append(elt.split(")")[0])
				if protID in dic_ecs: dic_ecs[protID]+= list_ECs
				else: dic_ecs[protID] = list_ECs
		else:
			if protID in all_protIDs:
				if protID in dic_ecs: dic_ecs[protID].append("NA")
				else: dic_ecs[protID] = ["NA"]
	tab.close()					
	max_name = max(len(name) for name in all_protIDs)
	out = open(output, "w")
	for protID, list_ECs in dic_ecs.items():
		for ec in list(set(list_ECs)):
			out.write( "{0} {1}\n".format(protID.ljust(max_name), ("\t"+ ec).ljust(15)))
	out.close()
		
if __name__ == '__main__':
	import subprocess, commands, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--fasta", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "uniprot fasta file")
	parser.add_argument("-t","--tab", dest = "tab", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "uniprot format tab file")
	parser.add_argument("-o", "--out",
	                    dest = "output",
	                    type = str,
	                    default = "[Input fasta filename without its extension]_std.conv",
	                    help = "The output file name")
	args = parser.parse_args()
	if args.fasta is None or args.tab is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input fasta filename without its extension]_std.conv"
	if args.output == default_output: args.output = get_output_name(args.fasta) + "_std.conv"
	convert_file(args.fasta, args.tab, args.output)
