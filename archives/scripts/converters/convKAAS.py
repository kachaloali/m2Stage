#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
import subprocess, commands
class convKAAS:
	def __init__(self, KoFile, outputFile):
		"""
		Check the type of inFile and convert it in a special format

		:Param inFile_ko     : The input file name (kaas.ko)
		:Param outputFile 	 : The name of the output file after the Converting
		
		:TypeOf inFile_ko    : str
		:TypeOf outputFile	 : str 
		
		:Seealso:: get_type_file()
		:Seealso:: convert_file()
		"""
		self.inFile_ko	  = str(KoFile)
		self.outputFile   = str(outputFile)
		self.typeFile     = self.check_type_file()
		
	def check_type_file(self):
		"""
		Check if arg is a valid file which is of type kaas output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : str
		"""
		typeFile , Error, output = False, "", []
		# Try to check if input file is the output of kaas program
		try:
			kaas_ko = open(self.inFile_ko, "r")
			output = [line.split("\t")[1][0] for line in kaas_ko if len(line.split("\t")) > 1]
			output = list(set(output))
			kaas_ko.close()
		except Exception as Error:	pass
		if Error == "" and len(output) == 1 and output[0] == "K": return True
		return typeFile
	
	def convert_file(self):
		"""
		Convert the input file in the format that we have described. 
		The format is to be checked in the output file of this function
	
		:Returns 		: no returns
		:Actions 		: creating a file
		"""	
		
		if not self.typeFile:
			print("\nThe file %s is not in the correct format!" % self.inFile_ko)
			exit(0)
		
		# This part concerns the analysis of the contents of the files in order to bring out  
		# interesting informations.
		if self.typeFile:
			# Get the longest name among all protein. Knowing the longest name will allow us 
			# to better  format displaying in the output file
			kaas_ko = open(self.inFile_ko, "r")
			names = [line.split("\t")[0] for line in kaas_ko]
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			kaas_ko.close()
			#
			start_time = time.time() # set time for the beginning of converting	
			progress(1, 4, start_time, "Converting file, Please wait...") # progress
			kaas_ko = open(self.inFile_ko, "r")
			dic_kaas_ko = {}
			for line in kaas_ko:
				if line != "\n":
					line = line.strip("\n").split()
					nameProtein = line[0]
					if "|" in nameProtein: nameProtein = nameProtein.split("|")[1]
					if len(line) > 1 :
						associate_ko = line[1]
						if nameProtein in dic_kaas_ko: 
							dic_kaas_ko[nameProtein].append(associate_ko)
						else: dic_kaas_ko[nameProtein] = [associate_ko]
					else: dic_kaas_ko[nameProtein] = []
			kaas_ko.close
		
			# Getting ko2ec [intermediary file: will be deleted at the end]
			curl_ko2ec = os.path.join("http://www.genome.jp", "kegg", "files", "ko2ec.xl")
			subprocess.check_output("wget -cq -P tmp " + curl_ko2ec, shell = True)
			progress(2, 4, start_time, "Converting file, Please wait...") # progress
			# We check if the downloading of intermediary files has been done properly
			if os.path.exists("tmp"):
				outputFile = open(self.outputFile, "w")
			else:
				text = "The files ko2ec.xl is not correctly downloaded please restart "
				text += "the program and make sure that you have a good connection"
				print ("%s\n", text)
				exit(0)
			ko2ec = open(os.path.join("tmp", "ko2ec.xl"), "r")
			dic_ko2ec = {}
			for line in ko2ec:
				if "EC:" in line:
					line = line.rstrip("\n").split()
					nameKo = line[0]
					associates_ecs = line[1:]
					associates_ecs[0] = associates_ecs[0].strip("[EC:")
					associates_ecs[-1] = associates_ecs[-1].strip("]")
					dic_ko2ec[nameKo] = associates_ecs
			ko2ec.close()
			# At the end of the program, the temporary file is deleted
			if os.path.exists("tmp"): print commands.getoutput("rm -r tmp")
			#
			progress(3, 4, start_time, "Converting file, Please wait...") # progress
			dic_associates_ECs = {}
			for nameProtein in dic_kaas_ko:
				associates_ecs = []
				for nameKo in dic_kaas_ko[nameProtein]:	
					if nameKo in dic_ko2ec:	associates_ecs+= dic_ko2ec[nameKo]
				dic_associates_ECs[nameProtein] = associates_ecs
			
			progress(4, 4, start_time, "Converting file, Please wait...") # progress
			# writing in output file
			for nameProtein, ECs in dic_associates_ECs.items():
				if ECs != []:
					for ec in ECs:
						outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name), \
						("\t"+ ec).ljust(15)))
				else: 
					outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name), \
					"\tNA".ljust(15)))
		outputFile.close()

if __name__ == '__main__':
	from convKAAS import convKAAS
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The kaas output file name")
	parser.add_argument("-o", "--out",
	                    dest = "output",
	                    type = str,
	                    default = "[Input filename without its extension]_kaas.conv",
	                    help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_kaas.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_kaas.conv"
	convKAAS(args.file, args.output).convert_file()
