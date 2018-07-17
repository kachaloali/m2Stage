#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import is_valid_file
from utils import get_output_name
if __name__ == '__main__':
	import commands, os, time
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#	
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--fasta", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The input fasta file name")
	parser.add_argument("-p", "--path", 
						dest = "path", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The path to interproscan program")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input fasta filename without its extension].tsv",
		                help = "The output file name")  
		                					
	args = parser.parse_args()
	if args.file is None or args.path is None: 
			parser.print_help()
			exit(0)
	#
	# Running the interproscan program
	default_output = "[Input fasta filename without its extension].tsv"
	if args.output == default_output: args.output = get_output_name(args.file) + ".tsv"
	#if os.path.exists(os.path.abspath("runIprScan.sh")): os.remove(os.path.abspath("runIprScan.sh"))
	# Record time.
	time_stamp = str(time.time())	
	script = open("runIprScan."+ time_stamp +".sh", "w")
	command = args.path +" -i "+ args.file +" -f TSV -pa -o "+ args.output
	command+= " -iprlookup -goterms -dp -crid UNIDAUPIF -cpu 30\n"
	script.write(command)
	script.close()
	print commands.getoutput("chmod +x runIprScan."+ time_stamp +".sh")
	print commands.getoutput("qsub -q all.q@@bigmem -shell yes -S /bin/bash -cwd runIprScan."+ time_stamp +".sh")
	
