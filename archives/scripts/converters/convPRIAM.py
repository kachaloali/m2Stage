#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
import subprocess, commands
class convPRIAM:
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
		Check if arg is a valid file which is of type priam output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : boolean
		"""
		# Try to check if input file is the output of Priam program
		typeFile , Error, output = False, "", []
		try:
			inputFile = open(self.inFile, "r")
			output = [l.split("\t")[3].strip(' \t\n\r')  for l in inputFile if len(l.split("\t")) > 1]
			output = list(set(output))
			inputFile.close()
		except Exception as Error:	pass
		if len(output) < 1 or len(output) > 2 or Error != "" : return typeFile
		elif Error == "" and "F" in output and "T" in output : return True
		elif Error == "" and ("F" in output or "T" in output): return True	
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
			inputFile = open(self.inFile, "r")
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format the display in the output file 
			names = [l.strip(' \t\n\r')  for l in inputFile if l.startswith(">")]
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			inputFile.close()
			#
			# In order to better understand the parsing progress of this block of code, please 
			# inspect the structure of the Priam file. The procedure is as follows: as long as we 
			# are not at the end of the file and there exist a protein ID, it is associated with 
			# its EC numbers.
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting
			outputFile = open(self.outputFile, "w")
			inputFile = open(self.inFile, "r")
			currentLine = inputFile.readline()
			while currentLine.startswith(">"):
				i+= 1
				progress(i, len(names), start_time, "Converting file, Please wait...")# progress
				#
				nameProtein = currentLine[1:].rstrip(' \t\n\r')
				if "|" in nameProtein: nameProtein = nameProtein.split("|")[1]
				nextLine = inputFile.readline()
				# To have a good text formatting in output file, the spacing for the first 
				# column is based on the length of the longest protein ID.
				ec = nextLine.split("\t")[0]
				proba = nextLine.split("\t")[2]
				kept = nextLine.split("\t")[3].strip(" \n\t\r")
				test = False
				if kept == "T":
					test = True
					outputFile.write( "{0} {1} {2}\n".format(nameProtein.rstrip().ljust(max_name), \
					("\t"+ ec.strip(' \t\n\r')).ljust(15), ("\t"+ proba.strip(' \t\n\r')).ljust(12)))
				nextLine = inputFile.readline()
				while nextLine.rstrip() != "":
					ec = nextLine.split("\t")[0]
					proba = nextLine.split("\t")[2]
					kept = nextLine.split("\t")[3].strip(" \n\t\r")
					if kept == "T":
						test = True
						outputFile.write( "{0} {1} {2}\n".format(nameProtein.rstrip().ljust(max_name), \
						("\t"+ ec.strip(' \t\n\r')).ljust(15), ("\t"+ proba.strip(' \t\n\r')).ljust(12)))
					nextLine = inputFile.readline()
				if not test:
					outputFile.write( "{0} {1} {2}\n".format(nameProtein.rstrip().ljust(max_name), \
					("\tNA").ljust(15), ("\tNone").ljust(12)))
				currentLine = inputFile.readline()
			inputFile.close()	
			outputFile.close()

if __name__ == '__main__':
	from convPRIAM import convPRIAM
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The Priam output file name")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input filename without its extension]_pri.conv",
		                help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_pri.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_pri.conv"
	convPRIAM(args.file, args.output).convert_file()
