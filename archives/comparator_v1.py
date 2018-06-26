#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->

def run_process(cmd):
    # Makes the actual system call with the command passed as parameter
   	subprocess.call(cmd, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)

def create_process(cmd):
    # Launches the process by sending the command arguments to a subroutine.
    # Returns a Process object that can then be queried for process status.
    p = Process(target = run_process, args = (cmd,))
    p.start()
    return p
	
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
	# Runners programs paths
	e2p2    = os.path.join(script_path,  "scripts", "runners", "runE2P2.py")
	priam   = os.path.join(script_path, "scripts", "runners", "runPriam.py")
	iprscan = os.path.join(script_path, "scripts", "runners", "runIprScan.py")
	#
	# Programs commands
	# Best params values for priam prediction: th, mp, mo = str(0.6), str(20), str(40)
	priam_command = ["python", priam, "-f", fasta, "-pt", "0.6", "-mp", "20", "-mo", "40"]
	e2p2_command    = ["python", e2p2, "-f", fasta]
	iprscan_command = ["python", iprscan, "-f", fasta]
	#
	# programs process: run all programs in parallel at the same time
	priam_process    = create_process(priam_command)
	e2p2_process     = create_process(e2p2_command)
	iprscan_process  = create_process(iprscan_command)
	#
	# Hold until the end of the last program. 
	# is_alive() method allows to check if the program is still running
	while e2p2_process.is_alive() or priam_process.is_alive() or iprscan_process.is_alive(): 
		time.sleep(5)
	#
	# Convert uniprot reference file
	convUniprot = os.path.join(script_path, "scripts", "converters_1", "convUniprot.py")
	command = ["python", convUniprot, "--fasta", fasta, "--tab", tab]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# Convert priam output files
	convPRIAM 	= os.path.join(script_path, "scripts", "converters_1", "convPRIAM.py")
	resPriam = os.path.join(fasta_path, "paj_PR_seqsECs.txt")
	command = ["python", convPRIAM, "--file", resPriam]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	resPriam_out_name = fasta_without_extension + "_pri.conv"
	shutil.move(os.path.join(fasta_path, "paj_PR_seqsECs_pri.conv"), resPriam_out_name)
	
	# Convert e2p2 output files
	convE2P2 	= os.path.join(script_path, "scripts", "converters_1", "convE2P2.py")
	e2p2_output = fasta_without_extension + ".pf.pf"
	command = ["python", convE2P2, "--file", e2p2_output]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	# Delete files that don't interest us much among the e2p2 output files 
	os.remove(fasta_without_extension + ".pf")
	os.remove(fasta_without_extension + ".pf.long")
	os.remove(fasta_without_extension + ".pf.orxn.pf")
	
	# Convert iprscan output file
	convIPRSCAN = os.path.join(script_path, "scripts", "converters_1", "convIprScan.py")
	iprscan_output = fasta_without_extension + ".tsv"
	command = ["python", convIPRSCAN, "--file", iprscan_output]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	#
	# Output files after converting 
	sp_file =    os.path.join(fasta_path, fasta_without_extension + "_std.conv")
	priam_file = os.path.join(fasta_path, fasta_without_extension + "_pri.conv")
	e2p2_file =  os.path.join(fasta_path, fasta_without_extension + ".pf_e2p2.conv")
	iprscan_file =  os.path.join(fasta_path, fasta_without_extension + "_iprscan.conv")
	#
	# Compare priam output file with sp reference file
	compPRIAM 	= os.path.join(script_path, "scripts", "comparators", "compPRIAM.py")
	command = ["python", compPRIAM, "--fasta", fasta, "--priam", priam_file, "--sp", sp_file]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# Compare e2p2 output file with sp reference file
	compE2P2 	= os.path.join(script_path, "scripts", "comparators", "compE2P2.py")
	command = ["python", compE2P2, "--fasta", fasta, "--e2p2", e2p2_file, "--sp", sp_file]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# Compare iprscan output file with sp reference file
	compIPRSCAN 	= os.path.join(script_path, "scripts", "comparators", "compIprScan.py")
	command = ["python", compIPRSCAN, "--fasta", fasta, "--iprscan", iprscan_file, "--sp", sp_file]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	exit(0)
	
	
	fasta_extension = args.fasta.split(".")[-1]
	converted_sp = args.fasta[:len(args.fasta) - len(fasta_extension) - 1] + "_std.conv"
	converted_priam = args.fasta[:len(args.fasta) - len(fasta_extension) - 1] + "_pri.conv"

	#runPriam(fasta, tab, threshold, min_pro_len, max_overlap)
	priam_prediction_dic, ref_dic = load_and_get_data(fasta, converted_sp, converted_priam)
	comp_and_write(fasta, priam_prediction_dic, ref_dic)
	get_stats(fasta, priam_prediction_dic, ref_dic, level)






    
    


