#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import is_valid_file
from utils import get_output_name
def restructure_fasta(fasta_file, output):
	"""
	The fasta file that the E2P2 program takes is a bit particular: Headers 
	in the FASTA file should begin with the sequence ID followed by a space.
	So this function allows to restructure the fasta file to be good for E2P2 
	program
	"""
	_fasta_file = open(fasta_file, 'r')
	new_fasta_file = open(output, 'w')
	currentLine = _fasta_file.readline()
	while True:
		if currentLine.startswith(">"):
			if "|" in currentLine:
				new_fasta_file.write(">" + currentLine.split("|")[1] +" \n")
			else:
				new_fasta_file.write(currentLine)
			seq = _fasta_file.readline()
			nextLine = _fasta_file.readline()
			while not nextLine.startswith(">"):
				seq+= nextLine
				nextLine = _fasta_file.readline()
				if nextLine == '': break
			new_fasta_file.write(seq)
		currentLine = nextLine
		if currentLine == '': break
	_fasta_file.close()
	new_fasta_file.close()
	
if __name__ == '__main__':
	import subprocess, commands, shutil, sys, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The input fasta file name")
	parser.add_argument("-p", "--path", 
						dest = "path", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The path to e2p2 program")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input fasta filename without its extension].pf",
		                help = "The output file name")  
	args = parser.parse_args()
	if args.file is None or args.path is None: 
			parser.print_help()
			exit(0)
	#
	# Retrieve fasta file directory path
	fasta_path = os.path.dirname(os.path.abspath(args.file))
	test_fasta = os.path.join(fasta_path, ".e2p2_test_file.fasta")
	
	# Restructure the fasta file to be good for E2P2 program
	restructure_fasta(args.file, test_fasta)
	
	# Running the E2P2 program
	default_output = "[Input fasta filename without its extension].pf"
	if args.output == default_output: args.output = get_output_name(args.file) + ".pf"
	command = ["python", args.path, "-i", test_fasta, "-o", args.output]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# Deleting the temporary and intermediate files
	e2p2_path = os.path.dirname(os.path.abspath(args.path))
	if os.path.exists(test_fasta): os.remove(test_fasta)
	if os.path.exists(os.path.join(e2p2_path, "run")): 
		shutil.rmtree(os.path.join(e2p2_path, "run"))	
