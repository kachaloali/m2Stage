# -*- coding: utf-8 -*-
import commands, os
if __name__ == '__main__':

	import os
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
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--query", dest = "query", 
						type = lambda arg: is_valid_file(parser, arg), help = "The input fasta file")               	    
	parser.add_argument("--database", dest = "database", help = "database")
	parser.add_argument("--evalue", dest = "evalue", help = "evalue")
	args = parser.parse_args()
	
	if args.query is None or args.database is None or args.evalue is None: 
			parser.print_help()
			exit(0)
	
	print commands.getoutput("rpsblast -i "+ args.query +" -d "+ args.database +" -e "+ args.evalue) 



































