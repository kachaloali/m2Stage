#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
import subprocess, commands
class convDETECT:
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
			output  = inputFile.readline()
			output = output.split("\t")
			output = list(set([elt.strip(" \n\t\r") for elt in output]))
			inputFile.close()
		except Exception as Error:	pass
		if len(output) != 5 or Error != "" : return typeFile
		elif Error == "" and "ID" in output and "EC" in output and "probability" in output\
		and "positive_hits" in output and "negative_hits" in output: return True
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
		# interesting informations.
		if self.typeFile:
			inputFile = open(self.inFile, "r")
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format the display in the output file 
			names = [l.split("\t")[0].strip(' \t\n\r')  for l in inputFile 
			if not l.startswith("ID")]
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			inputFile.close()
			outputFile = open(self.outputFile, "w")
			inputFile = open(self.inFile, "r")
			for line in inputFile:
				if not line.startswith("ID"):
					nameProtein = line.split("\t")[0].strip(' \t\n\r')
					if "|" in nameProtein: nameProtein = nameProtein.split("|")[1]
					ec = line.split("\t")[1].strip(' \t\n\r')
					outputFile.write( "{0} {1}\n".format(nameProtein.ljust(max_name), \
					("\t"+ ec.strip(' \t\n\r')).ljust(15)))
			inputFile.close()	
			outputFile.close()

if __name__ == '__main__':
	from convDETECT import convDETECT
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The DETECT output file name")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input filename without its extension]_detect.conv",
		                help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_detect.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_detect.conv"
	convDETECT(args.file, args.output).convert_file()
