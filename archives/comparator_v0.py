if __name__ == '__main__':
	import subprocess, commands, shutil, glob, time, sys, os
	from multiprocessing import Process
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	def is_valid_file(parser, arg):
		"""
		Check if arg is a valid file that already exists on the file system.
		"""
		arg = os.path.abspath(arg)
		if not os.path.exists(arg):
			parser.error("The file %s does not exist!" % arg)
		else: return arg

	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--fas", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")           	    
	parser.add_argument("--tab", dest = "tab", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input Swiss-Prot tab file")
	parser.add_argument("--priam", dest = "priam", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The output priam file")
	parser.add_argument("--e2p2", dest = "e2p2", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The output e2p2 file")
	parser.add_argument("--iprscan", dest = "iprscan", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The InterproScan output file")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	args = parser.parse_args()
	
	if args.fasta is None or args.tab is None or args.level is None: 
		parser.print_help()
		exit(0)
	elif args.level == "0": 
		print("Oops! Level was no valid number. Try again with level in [1..4]\n")
		exit(0)
		
	fasta, tab, level = args.fasta, args.tab, args.level
	fasta_extension = fasta.split(".")[-1]
	fasta_without_extension = fasta[:len(fasta) - len(fasta_extension) - 1]
	#
	# Retrieve fasta file path
	fasta_path = os.path.dirname(os.path.realpath(fasta))
	#
	# Retrieve this script directory path
	script_path = os.path.dirname(os.path.abspath(__file__))
	#
	# Convert uniprot reference file
	convUniprot = os.path.join(script_path, "scripts", "converters", "convUniprot.py")
	command = ["python", convUniprot, "--fasta", fasta, "--tab", tab]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	sp_file =    os.path.join(fasta_path, fasta_without_extension + "_std.conv")
	
	# Convert priam output files
	if args.priam is not None:
		convPRIAM 	= os.path.join(script_path, "scripts", "converters", "convPRIAM.py")
		resPriam = os.path.join(fasta_path, "paj_PR_seqsECs.txt")
		command = ["python", convPRIAM, "--file", resPriam]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		resPriam_out_name = fasta_without_extension + "_pri.conv"
		shutil.move(os.path.join(fasta_path, "paj_PR_seqsECs_pri.conv"), resPriam_out_name)
		priam_file = os.path.join(fasta_path, fasta_without_extension + "_pri.conv")
		# Compare priam output file with sp reference file
		compPRIAM 	= os.path.join(script_path, "scripts", "comparators", "compPRIAM.py")
		command = ["python", compPRIAM, "--fasta", fasta, "--priam", priam_file, "--sp", sp_file]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		if os.path.exists(priam_file) : os.remove(priam_file)
		
	# Convert e2p2 output files
	if args.e2p2 is not None:
		convE2P2 	= os.path.join(script_path, "scripts", "converters", "convE2P2.py")
		e2p2_output = fasta_without_extension + ".pf.pf"
		command = ["python", convE2P2, "--file", e2p2_output]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		e2p2_file =  os.path.join(fasta_path, fasta_without_extension + ".pf_e2p2.conv")
		# Compare e2p2 output file with sp reference file
		compE2P2 	= os.path.join(script_path, "scripts", "comparators", "compE2P2.py")
		command = ["python", compE2P2, "--fasta", fasta, "--e2p2", e2p2_file, "--sp", sp_file]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		if os.path.exists(e2p2_file) : os.remove(e2p2_file)
	
	# Convert iprscan output file
	if args.iprscan is not None:
		convIPRSCAN = os.path.join(script_path, "scripts", "converters", "convIprScan.py")
		iprscan_output = fasta_without_extension + ".tsv"
		command = ["python", convIPRSCAN, "--file", iprscan_output]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		iprscan_file =  os.path.join(fasta_path, fasta_without_extension + "_iprscan.conv")
		# Compare iprscan output file with sp reference file
		compIPRSCAN 	= os.path.join(script_path, "scripts", "comparators", "compIprScan.py")
		command = ["python", compIPRSCAN, "--fasta", fasta, "--iprscan", iprscan_file, "--sp", sp_file]
		subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
		if os.path.exists(iprscan_file) : os.remove(iprscan_file)
	
	if os.path.exists(sp_file) : os.remove(sp_file)
	print "\tProcess finished\n"






    
    


