#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import is_valid_file
from utils import get_output_name
if __name__ == '__main__':
	import commands, shutil, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	# Best params values for priam prediction: threshold = 0.6, min_pro_len = 20, max_overlap = 40 
	# Default parameters: -pt 0.5 -mo -1 -mp 70 -cc T -cg T -e T
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", dest = "ff", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")
	parser.add_argument("-pp", 
						dest = "pp", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The path to Priam program")
	parser.add_argument("-pr", dest = "pr", help = "The path to Priam Releases")	        
	parser.add_argument("-pt", dest = "th", help = "threshold", default = "0.5")
	parser.add_argument("-mp", dest = "mp", default = "-1",
	help = "Minimal length proportion of a profile that must be matched to consider it")
	parser.add_argument("-mo", dest = "mo",  default = "70",
	help = "Maximum overlap length between the matches of two profiles")
	
	args = parser.parse_args()
	if args.ff is None or args.pp is None or args.pr is None or args.th is None or args.mp is None: 
			parser.print_help()
			exit(0)
	
	# directory that will contain intermediates, temporary files and Results
	priam_folder = os.path.dirname(os.path.abspath(args.pp))
	pathToResults = os.path.join(priam_folder, 'resPriam')
	
	# running command to execute priam program with the chosen parameters
	if os.path.exists(pathToResults): shutil.rmtree(pathToResults)
	print commands.getoutput("java -jar "+ args.pp +" -i "+ args.ff + " -p "+ args.pr + " -od "
							+ pathToResults + " -n output -pt "+ args.th +" -mo "+ args.mp +" -mp "
							+ args.mo +" -cc F -cg T -e T"
							)
	# Retrieve fasta file folder
	fasta_dir= os.path.dirname(args.ff)
	resPriam = os.path.join(priam_folder, "resPriam", "RESULTS", "paj_output_seqsECs.txt")
	shutil.move(resPriam, fasta_dir)
