# -*- coding: utf-8 -*-
"""
Created on wednesday Jan 31 11:04:57 2018

@author: Kachalo Ali
"""

def download(taxon, out_format):
	"""
	allows to download the taxon file of an organisms giving the  
	names of the taxon and the format to download [fasta or tab]
	from swiss-prot database
	"""
	taxon = taxon.strip().split()
	norm_taxon, output_file = "", ""
	if len(taxon) > 1:
		norm_taxon = taxon[0]
		output_file = taxon[0]
		for elt in taxon[1:]:
			norm_taxon+= "+" + elt
			output_file+= elt
	else: 
		norm_taxon = taxon[0]
		output_file = taxon[0]
	query = "'http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+taxonomy:"
	query+= norm_taxon +"&format="+ out_format +"'"
	output_file+= "."+ out_format
	try:
		subprocess.check_output(['bash', '-c' , "curl " + query + " > "+ output_file])
	except Exception as error: print error
	if os.stat(output_file).st_size == 0:
		if os.path.exists(output_file):	print commands.getoutput("rm " + output_file)
		raise Exception ('Bad taxon name')

if __name__ == '__main__':
	import subprocess, commands, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("--taxon", dest = "taxon", type = str, 
						help = "taxon name")
	parser.add_argument("--format", dest = "format", type = str, 
						help = "The output format [fasta or tab]")
	args = parser.parse_args()
	if args.taxon is None or args.format is None: 
			parser.print_help()
			exit(0)
	taxon, out_format = args.taxon, args.format 
	#
	#
	download(taxon, out_format)
