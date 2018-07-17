#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import is_valid_file
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
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", dest = "fasta", 
						type = lambda arg: is_valid_file(parser, arg), 
					   	help = "The input fasta file")           	    
	parser.add_argument("--e2p2", dest = "e2p2",
					   	help = "switch on running E2P2 program or No. The E2P2 path is required")
	parser.add_argument("--iprscan", dest = "iprscan",
					   	help = "switch on running InterproScan or No. The interproscan path is required")
	parser.add_argument("--priam", dest = "priam",
					   	help = "switch on running Priam Program or No. The Priam path and Priam releases are required separated by comma")
	args = parser.parse_args()
	if args.fasta is None: 
		parser.print_help()
		exit(0)
	list_running_programs = []
	fasta_extension = args.fasta.split(".")[-1]
	fasta_without_extension = args.fasta[:len(args.fasta) - len(fasta_extension) - 1]
	#
	# Retrieve fasta file path
	fasta_path = os.path.dirname(os.path.realpath(args.fasta))
	#
	# Retrieve this script directory path
	script_path = os.path.dirname(os.path.abspath(__file__))
	
	if args.e2p2 is not None:
		e2p2_path = args.e2p2
		e2p2    = os.path.join(script_path,  "scripts", "runners", "runE2P2.py")
		e2p2_command    = ["python", e2p2, "-f", args.fasta, "-p", e2p2_path]
		e2p2_process    = create_process(e2p2_command)
		list_running_programs.append(e2p2_process)
	
	if args.priam is not None:
		if len(args.priam.split(",")) > 1: priam_path, priam_releases = args.priam.split(",")
		else:
			parser.print_help()
			exit(0)
		if priam_path.rstrip(" \n\r\t") == "" or priam_releases.rstrip(" \n\r\t") == "":
			parser.print_help()
			exit(0)
		priam   = os.path.join(script_path, "scripts", "runners", "runPriam.py")
		priam_command = ["python", priam, "-f", args.fasta, "-pp", priam_path, "-pr", priam_releases]
		priam_process = create_process(priam_command)
		list_running_programs.append(priam_process)
	
	if args.iprscan is not None:
		iprscan = os.path.join(script_path, "scripts", "runners", "runIprScan.py")
		iprscan_path = args.iprscan
		#iprscan_path = os.path.join("/", "usr", "local", "bioinfo", "interproscan", "interproscan.sh")
		iprscan_command = ["python", iprscan, "-f", args.fasta, "-p", iprscan_path]
		iprscan_process = create_process(iprscan_command)
		list_running_programs.append(iprscan_process)
	
	if list_running_programs == []:
		parser.print_help()
		exit(0)
	
	# Hold until the end of the last program. 
	# is_alive() method allows to check if the program is still running
	for i in range(len(list_running_programs)):
		while list_running_programs[i].is_alive():	time.sleep(5)
	print "\tProcess finished\n"
