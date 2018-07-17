# -*- coding: utf-8 -*-
import ftplib, subprocess, commands, os

def download_db(nb_file):
	"""
	This function allows to download the few files of nr database and 
    to extract the contents in the current folder.
	"""
	host = "ftp.ncbi.nlm.nih.gov"
	connect = ftplib.FTP(host, 'anonymous', 'anonymous')
	connect.cwd('/blast/db/')
	if not os.path.exists("nr"): 
		print commands.getoutput("mkdir nr")
		os.chdir("nr/")
	if nb_file > 0:
		for i in range(nb_file):
			# just for formatting: 0 => '00', 1 => '01', 2 => '02',...,9 => '09'
			if i < 10: count = "{:02x}".format(i)
			else: count = i
			filename = "nr."+ count +".tar.gz"
			with open(filename, 'wb') as infile:
				connect.retrbinary('RETR %s' % filename, infile.write)
		print commands.getoutput("cat nr.*.tar.gz | tar -zxvi -f - -C .")
	else: 
		print commands.getoutput("wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'")
		print commands.getoutput("cat nr.*.tar.gz | tar -zxvi -f - -C .")

if __name__ == '__main__':
	import subprocess, commands, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("--number", dest = "number", type = int, default = "-1",
						help = "The number of the nr files you want to retrieve")
	args = parser.parse_args()
	if args.number is None: 
			parser.print_help()
			exit(0)
	download_db(args.number)             

