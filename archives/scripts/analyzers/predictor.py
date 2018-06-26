#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
"""
	Perform annotation from multiple results from multiple annotation tools. The principle is as follows:
	1. From the score threshold of f1-score received as parameter, we build a list of potentials ec numbers.
	i.e the list of ec having a f1-score greater than or equal to the threshold.

	2. The potentials ec numbers predicted by all programs are retained as valid ecs.
	3. The partial intersections of at least 3 programs are retained as valid ecs: intersections containing 
	only potentials ec numbers.
	4. The partial intersections of a program with the reference (swissprot) are retained as valid ecs: 
	intersections containing only potentials ec numbers.
	5. An inventory of proteins is made, associating them with all the ecs assigned to it by the programs.
		For each protein and for each ec that has been associated with it:
			. If the ec exists in the list of valid ecs, then this protein is definitively associated with 
			this ec number.
			. Otherwise, the ec number will be removed.
"""

from utils import *

if __name__ == '__main__':
	import os, csv, itertools
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--prog", dest = "prog",
						type = lambda arg: are_valid_file(parser, arg),
					   	help = "List of the converted programs files separated by comma")
	parser.add_argument("--sp", dest = "sp",
						type = lambda arg: is_valid_file(parser, arg),
						help = "The converted swiss-prot file")
	parser.add_argument("--plud", dest = "plud",
						help = "Path to the limputils directory")
	parser.add_argument("--pas", dest = "pas",
						help = "Path to the analyzer.py script")
	parser.add_argument("--level", dest = "level", default = "4", help = "The level of ec classe")
	parser.add_argument("--score", dest = "score", default = "0.5", help = "The f1 score value")
	# minimal programs intersection : it is the required minimal number of programs constituting an intersection 
	# to consider it. By default its value is three (3).
	parser.add_argument("--mipci", dest = "mipci", default = "3", 
						help = "The minimum number of programs constituting an intersection to be considered")
	# minimal programs intersection : it is the required minimal number of programs constituting an intersection 
	# to consider it. By default its value is three (3).
	parser.add_argument("--miprec", dest = "miprec", default = "1", 
						help = "The minimum number of programs repporting an ec number to be considered as FP") 
	parser.add_argument("--out", dest = "out", default = "annotation.out.txt", help = "The output file name")
	
	args = parser.parse_args()
	if args.prog is None or args.sp is None or args.score is None or args.plud is None or args.pas is None: 
		parser.print_help()
		exit(0)
		
	j = 1 # set j for counting progress bar status
	start_time = time.time() # set time for the beginning
	args.prog = args.prog.split(",")
	prog_dic = {}
	progress(j, 2*len(args.prog)+2, start_time, "validating ec numbers...")# progress
	
	
	# Determine the list of ec numbers with f1-score greater than the thresold (potentials ec numbers)
	# The ec numbers predicted by the programs are first validated with the library of validation: lipmutils.
	time_stamp = time.time()
	outfile = str(args.prog[0]) + "."+ str(time_stamp) +".snp."+ args.level +".txt"
	command = ["python", args.pas, "--prog", args.prog[0], "--sp", args.sp, "--path", args.plud, "--level", args.level ,"--out", outfile]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	dic_global_ec = {}
	ec_list0, ec_list1 = list(), list()
	with open(outfile, 'rb') as csvfile:
		for row in csvfile:
			if not row.startswith("EC"):
				ec = row.split("\t")[0].strip(" \n\t\r")
				ec_list0.append(ec)
				if row.split("\t")[3] >= args.score:
					ec = row.split("\t")[0].strip(" \n\t\r")
					#score = row.split("\t")[3].strip(" \n\t\r")
					ec_list1.append(ec)
	prog_name = "p0"
	prog_dic[prog_name] = ec_list1
	for ec in set(ec_list0): dic_global_ec[ec] = 1
	os.remove(outfile)
	
	i = 1
	while i < len(args.prog):
		j+= 1
		progress(j, 2*len(args.prog)+2, start_time, "validating ec numbers...")# progress
		ec_list0, ec_list1 = list(), list()
		time_stamp = time.time()
		outfile = str(args.prog[i]) + "."+ str(time_stamp) +".snp."+ args.level +".txt"
		command = ["python", args.pas, "--prog", args.prog[i], "--sp", args.sp, "--path", args.plud, "--level", args.level ,"--out", outfile]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		with open(outfile, 'rb') as csvfile:
			for row in csvfile:
				if not row.startswith("EC"):
					ec = row.split("\t")[0].strip(" \n\t\r")
					ec_list0.append(ec)
					if row.split("\t")[3] >= args.score:
						ec = row.split("\t")[0].strip(" \n\t\r")
						#score = row.split("\t")[3].strip(" \n\t\r")
						ec_list1.append(ec)
		prog_name = "p" + str(i)
		prog_dic[prog_name] = ec_list1
		for ec in set(ec_list0): 
			if ec in dic_global_ec: dic_global_ec[ec]+= 1
			else: dic_global_ec[ec] = 1
		i+= 1
		os.remove(outfile)

	# Try to get list of ec numbers reported in swiss-prot database. The ec numbers directly from the swissprot database do not 
	# need to be validated by the lipmutils library. They are already up to date
	ec_list = [line.split("\t")[1].strip(" \n\t\r") for line in open(args.sp, "r") if "NA" not in line]
	if 4 - int(args.level) == 0:
		ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
	else:
		ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
		if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
	#prog_dic["sp"] = ec_list
	
	
	# All possible partitions are calculated between lists of ecs numbers from prediction programs. The partitions here are exactly 
	# the same as when applying a Venn diagramm to calculate the intersections between many sets.
	j+= 1
	progress(j, 2*len(args.prog)+2, start_time, "setting partitions...")# progress
	powerset1 = itertools.chain.from_iterable(itertools.combinations(prog_dic.keys(), y) for y in range(len(prog_dic) + 1))
	powerset2 = itertools.chain.from_iterable(itertools.combinations(prog_dic.values(), y) for y in range(len(prog_dic) + 1))
	titles = [part for part in powerset1 if any(map(None, part))]
	sources = [part for part in powerset2 if any(map(None, part))]
	intersect_dic = {}
	for elt in (zip(titles, sources)):
		key = "∩".join(str(x) for x in list(elt[0]))
		intersect_dic[key] = reduce(lambda x, y: set(x).intersection(y), list(elt[1])) 
		 
	# List of valid ec numbers are calculated according to the parameters mi that's means: The minimum number of programs 
	# constituting an intersection to be considered.
	valid_ecs = set()
	for key, value in intersect_dic.items():
		if len(key.split("∩")) >= int(args.mipci): valid_ecs.update(value)
		if len(key.split("∩")) == len(prog_dic): super_intersect_value = value
	for key, value in intersect_dic.items():
		if len(key.split("∩")) == 1:
			extreme_ec_list = set(value) - set(super_intersect_value)
			valid_ecs.update(extreme_ec_list)
	

	# An inventory of proteins ID is made, associating them with all the ecs assigned to it by the programs of predictions.
	# The rule is as follow:
	# For each protein and for each ec that has been associated with it:
	#	. If the ec exists in the list of valid ecs, then this protein is definitively associated with this ec number.
	#	. Otherwise, the ec number will be removed.
	merged_dic = {}
	dic = get_prog_data(args.prog[0], args.plud)
	for protID, ec_list in dic.items():
		if 4 - int(args.level) == 0:
			ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
		else:
			ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
			if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
		dic_ec = {}
		for ec in ec_list: 
			if ec in valid_ecs: dic_ec[ec] = 1
		merged_dic[protID] = dic_ec
	#nb_prot_annot = len([protID for protID in merged_dic if merged_dic[protID] != []])
	#print str(os.path.basename(args.prog[0]))+ ": " + str(nb_prot_annot) + "\n"
	i = 1
	while i < len(args.prog):
		j+= 1
		progress(j, 2*len(args.prog)+2, start_time, "proteins inventory...")# progress
		dic0, dic1 = get_prog_data(args.prog[i], args.plud), {}
		for protID, ec_list0 in dic0.items():
			if 4 - int(args.level) == 0:
				ec_list = [ec for ec in ec_list0 if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
			else:
				ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list0 \
				if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
			dic1[protID] = list(set(ec_list))
		
		for protID, ec_list in dic1.items():
			if protID in merged_dic:
				for ec in ec_list:
					if ec in merged_dic[protID]: 
						merged_dic[protID][ec]+= 1
					elif ec in valid_ecs:
						merged_dic[protID][ec] = 1
				#merged_dic[protID] = list(set(ec_list))
			else:
				dic_ec = {}
				for ec in ec_list: 
					if ec in valid_ecs: dic_ec[ec] = 1
				merged_dic[protID] = dic_ec
		#nb_prot_annot = len([protID for protID in dic1 if dic1[protID] != []])
		#print str(os.path.basename(args.prog[i]))+ ": " + str(nb_prot_annot) + "\n"
		i+= 1
	
	ref_dic, dic0 = get_sp_data(args.sp), {}
	for protID, ec_list in ref_dic.items():
		if 4 - int(args.level) == 0:
			ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(args.level)]
		else:
			ec_list = ['.'.join(ec.split('.')[:-4 + int(args.level)]) for ec in ec_list \
			if len([x for x in ec.split(".")[:-4 + int(args.level)] if x != "-"]) == int(args.level)]
		dic0[protID] = ec_list
	#nb_prot_annot = len([protID for protID in dic0 if dic0[protID] != []])
	#print str(os.path.basename(args.sp))+ ": " + str(nb_prot_annot) + "\n"
	
	
	# The last step is that of the annotation, it consists in writing the result in an output file to be visualized.
	nb_prot_new_annot = 0
	j+= 1
	progress(j, 2*len(args.prog)+2, start_time, "annotating phase...")# progress
	annot_file = open(args.out, "w")
	ref_dic = get_sp_data(args.sp)
	for protID, dic_ec in merged_dic.items():
		if protID in ref_dic:
			ref_ec_list = ref_dic[protID]
			ec_list = [ec for ec in dic_ec]
			#new_ec_list = [ec for ec in ec_list if ec in valid_ecs]
			if ec_list != []:
				for ec in ec_list:
					if dic_ec[ec] == dic_global_ec[ec]: 
						annot_file.write("{0} {1}\n".format(str(protID).ljust(12), "\t"+str(ec).ljust(12)))
					else:
						score_ec = float(dic_ec[ec]) / float(dic_global_ec[ec])
						print score_ec
						if score_ec >= args.miprec:
							print ec, score_ec
							annot_file.write("{0} {1}\n".format(str(protID).ljust(12), "\t"+str(ec).ljust(12)))
				nb_prot_new_annot+= 1
	#print "number of news annotations: "+ str(nb_prot_new_annot) + "\n"
	annot_file.close()
	j+= 1
	progress(j, 2*len(args.prog)+2, start_time, "annotating phase...")# progress
	time.sleep(1)

	
	
	
	
