#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
import subprocess, commands, csv
class convIprScan:
	def __init__(self, inFile, outputFile):
		"""
		Check the type of inFile and convert it in a special format

		:Param inFile     : The input file name
		:Param outputFile : The name of the output file after the Converting
		
		:TypeOf inFile    : str
		:TypeOf outputFile: str 
		
		:Seealso:: check_type_file()
		:Seealso:: convert_file()
		"""
		self.inFile	  = str(inFile)
		self.outputFile = str(outputFile)
		self.typeFile = self.check_type_file()
		
	def check_type_file(self):
		"""
		Check if arg is a valid file which is of type InterproScan output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : str
		"""
		# Try to check if input file is the output of Iprscan program
		typeFile , Error, output = False, "", []
		try:
			inputFile = open(self.inFile, "r")
			output = [l.split("\t")[9] for l in inputFile if len(l.split("\t")) > 9]
			output = list(set(output))
			inputFile.close()
		except Exception as Error:	pass
		if len(output) < 1 or len(output) > 1 or Error != "": return typeFile
		elif Error == "" and len(output) == 1 and output[0] == "T": return True		
		return typeFile
	
	def convert_file(self):
		"""
		Convert the input file in the format that we have described. 
		The format is to be checked in the output file of this function
	
		:Returns 		: no returns
		:Actions 		: creating a file
		"""	
		if not self.typeFile:
			print("\nThe file %s is not in the correct format!" % self.inFile)
			exit(0)
		
		# This part concerns the analysis of the contents of the files in order to bring out  
		# interesting informations. Thus, to optimize and  speed up the processing, 
		# bash commands have been invoked. 
		#
		if self.typeFile:
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format the display in the output file 
			iprscan = open(self.inFile, "r")
			names = list(set([line.split("\t")[0] for line in iprscan]))
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			iprscan.close()
			#
			start_time = time.time() # set time for the beginning of converting
			progress(1, 5, start_time, "Converting file, Please wait...") # progress
			#
			# Getting interpro2go and ec2go [intermediary files: will be deleted at the end]
			curl_ipr2go = os.path.join("http://www.geneontology.org","external2go","interpro2go")
			curl_ec2go = os.path.join("http://www.geneontology.org", "external2go", "ec2go")
			subprocess.check_output("wget -cq -P tmp " + curl_ipr2go, shell = True)
			subprocess.check_output("wget -cq -P tmp " + curl_ec2go, shell = True)
			#
			# Check if the downloading of intermediary files has been done properly
			if os.path.exists("tmp"):
				outputFile = open(self.outputFile, "w")
			else:
				text = "The files interpro2go and ec2go are not correctly downloaded please "
				text += "restart the program and make sure that you have a good connection"
				print ("%s\n", text)
				exit(0)
			#
			progress(2, 5, start_time, "Converting file, Please wait...") # progress
			csvfile = open(self.inFile, 'rb')
			inputFile = csv.reader(csvfile, delimiter = "\t")
			# In the input file, each protID is associated with its corresponding iprID (s)
			dic_associated_IPRs = {}
			for line in inputFile:
				protID = line[0]
				if len(line) > 11:
					iprID = line[11]
					if protID in dic_associated_IPRs: dic_associated_IPRs[protID].append(iprID)
					else: dic_associated_IPRs[protID] = [iprID]
				else: dic_associated_IPRs[protID] = ["NA"]
			#
			# In the file interpro2go, each iprID is associated with its corresponding GO Terms
			ipr2go_file = open(os.path.join("tmp", "interpro2go"), "r")
			dic_IPRs2GOs = {}
			for line in ipr2go_file:
				if "InterPro" in line:
					iprID = line.split()[0].split(":")[1]
					goTerm = line.split("; ")[-1].strip(" \n\t\n")
					if iprID in dic_IPRs2GOs: dic_IPRs2GOs[iprID].append(goTerm)
					else: dic_IPRs2GOs[iprID] = [goTerm]
			#
			# Each protID is associated with its corresponding GO Terms
			dic_associated_GOs = {}
			for protID in dic_associated_IPRs:
				IPRs = dic_associated_IPRs[protID]
				GOs = [dic_IPRs2GOs[ipr] for ipr in IPRs if ipr!="NA" and ipr in dic_IPRs2GOs]
				if GOs != []: GOs = reduce(lambda list_1, list_2: list_1 + list_2, GOs)
				if protID in dic_associated_GOs: 
					dic_associated_GOs[protID]+= GOs
				else: dic_associated_GOs[protID] = GOs
				if GOs == []: dic_associated_GOs[protID] = []
			#
			progress(3, 5, start_time, "Converting file, Please wait...") # progress
			# In the file ec2go, each GO is associated with its corresponding ECs
			ec2go = open(os.path.join("tmp", "ec2go"), "r")
			dic_GOs_ECs = {}
			for line in ec2go:
				if "EC:" in line:
					GOID = line.split(";")[-1].strip("\n").strip(' \n\t\r')
					EC = line.split(">")[0].split(":")[1].strip(" \n\t\r")
					if len(EC.split(".")) < 4:
						for i in range(4 - len(EC.split("."))):	EC+= ".-"
					if GOID in dic_GOs_ECs: dic_GOs_ECs[GOID].append(EC)
					else: dic_GOs_ECs[GOID]= [EC]
			ec2go.close()
			#
			progress(4, 5, start_time, "Converting file, Please wait...") # progress
			# Finally each protID is associated with its corresponding EC numbers
			dic_associates_ECs = {}
			for protID, GOs in dic_associated_GOs.items():
				associates_ECs = []
				for go in GOs:
					if go in dic_GOs_ECs:
						associates_ECs+= dic_GOs_ECs[go]
				dic_associates_ECs[protID] = list(set(associates_ECs))
			#
			progress(5, 5, start_time, "Converting file, Please wait...") # progress
			# writing the ProtID with its associated ECs in the output file
			outputFile = open(self.outputFile, "w")
			for protID, ECs in dic_associates_ECs.items():
				if ECs != []:
					for ec in ECs:
						outputFile.write( "{0} {1}\n".format(protID.ljust(max_name), \
						"\t"+ ec.ljust(15)))	
				else: 
					outputFile.write("{0} {1}\n".format(protID.ljust(max_name), \
					"\tNA".ljust(15)))	
			outputFile.close()
		#	
		# At the end of the program, the temporary files are deleted
		if os.path.exists("tmp"): print commands.getoutput("rm -r tmp")

if __name__ == '__main__':
	from convIprScan import convIprScan
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The interproscan output file name")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input filename without its extension]_iprscan.conv",
		                help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_iprscan.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_iprscan.conv"
	convIprScan(args.file, args.output).convert_file()
