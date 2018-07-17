# -*- coding: utf-8 -*-
if __name__ == '__main__':
	import commands, os
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
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")               	    
	parser.add_argument("--db", dest = "database", help = "database")
	parser.add_argument("--out", dest = "output", help = "output aligments")
	parser.add_argument("--nb_it", dest = "nb_it",  default = "5", help = "number of iterations for psiblast")
	parser.add_argument("--pssm", dest = "pssm", help = "output pssm checkpoint matrix")
	parser.add_argument("--evalue", dest = "evalue", help = "evalue")
	args = parser.parse_args()
	
	if args.query is None or args.database is None or args.output is None or args.nb_it is None or args.pssm is None or args.evalue is None: 
			parser.print_help()
			exit(0)
	from Bio import SeqIO
	records = list(SeqIO.parse(args.query, "fasta"))
	for i in range(len(records)):
		file_ = open("train/train_"+ str(i+1)+".fasta", "w")
		file_.write(str(">"+ records[i].id).strip("\n") +"\n")
		file_.write(str(records[i].seq))
		file_.close()	
		print commands.getoutput("psiblast -query train/train_"+str(i+1)+".fasta -db "+ args.database +" -out out/align_"+str(i+1)+
		".fasta -num_iterations "+ args.nb_it +" -out_pssm pssm/pssm_"+str(i+1)+".chk -evalue "+ args.evalue)



































