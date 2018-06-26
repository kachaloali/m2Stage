# -*- coding: utf-8 -*-
"""
Created on wednesday Jan 31 12:00:47 2018

@author: Kachalo Ali
"""

def  runBlast(_file, output, db):
	os.chdir(db)
	print commands.getoutput("blastp -query "+ _file +" -db "+db+" -outfmt 5 -out "+ output)

def runBlast2Go(_xml):
	exit(0)
	print commands.getoutput("blast2go -in "+ _xml +" -annot -v")		

if __name__ == '__main__':
	from Bio import SeqIO
	from Bio.Blast import NCBIWWW
	import subprocess, commands, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	def is_valid_file(parser, arg):
		"""
		Check if arg is a valid file that already exists on the file system.

		Parameters
		----------
		parser : argparse object
		arg    : str

		Returns
		-------
		arg
		"""
		arg = os.path.abspath(arg)
		if not os.path.exists(arg):
			parser.error("The file %s does not exist!" % arg)
		else: return arg

	def get_output_name(parser):
		"""
		Returns the file name by replacing its old extension with the "_blast.xml" one.
		"""
		arg = parser.parse_args()
		infile_extension = arg.fasta.split(".")[-1]
		return arg.fasta[:len(arg.fasta) - len(infile_extension) - 1] + "_blast.xml"
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--fasta", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "uniprot fasta file")
	parser.add_argument("-db", dest = "database",
						help = "The database that you are using")
	parser.add_argument("-o", "--out",
	                    dest = "output",
	                    type = str,
	                    default = get_output_name(parser),
	                    help = "The output file name")
	args = parser.parse_args()
	if args.fasta is None: 
			parser.print_help()
			exit(0)
	fasta, database, output = args.fasta, args.database, args.output 
	#
	#
	runBlast(fasta, output, database)
	runBlast2Go(output)

