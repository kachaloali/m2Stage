#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
class convE2P2:
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
		:TypeOf typeFile : str
		"""
		typeFile, output = False, []
		# Try to check if input file is the output of e2p2 program
		e2p2 = open(self.inFile, "r")
		output = []
		for line in e2p2:
			if line.startswith("ID"): output.append(line.split()[0])
			elif line.startswith("NAME"): output.append(line.split()[0])
			elif line.startswith("PRODUCT-TYPE"): output.append(line.split()[0])
			elif line.startswith("EC"): output.append(line.split()[0])
		output = list(set(output))
		e2p2.close()
		
		if "ID" in output and "NAME" in output and "PRODUCT-TYPE" in output and "EC" in output: 
			return True
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
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format displaying in the output file
			inputFile = open(self.inFile, "r")
			names = [line.split("\t")[1].strip("\n") for line in inputFile if line.startswith("ID")]
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			inputFile.close()
		
			# In order to better understand the parsing progress of this block of code, please 
			# inspect the structure of the E2P2 file. The procedure is as follows: as long as we 
			# are not at the end of the file and there exist a protein ID, we recover this ID and 
			# then extract the EC numbers that have been associated with it.
			inputFile = open(self.inFile, "r")
			outputFile = open(self.outputFile, "w")
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting
			currentLine = inputFile.readline()
			while True:
				if currentLine.startswith("ID"):
					i+= 1
					progress(i, len(names), start_time, "Converting file, Please wait...")# progress
					#
					nameProtein = currentLine.split('\t')[1].rstrip('\n')
					nextLine = inputFile.readline()
					canPass = True
					while not nextLine.startswith("EC"):
						nextLine = inputFile.readline()
						if nextLine.startswith("//"): 
							canPass = False
							break
					# To have a good text formatting in output file, the spacing for the first 
					# column is based on the length of the longest protein ID.
					if canPass:
						list_ = nextLine.split("\t")
						outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name), \
						("\t"+ list_[1].rstrip('\n')).ljust(15)))
						nextLine = inputFile.readline()
						while nextLine.startswith("EC"):
							list_ = nextLine.split("\t")
							outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name), \
							("\t"+ list_[1].rstrip('\n')).ljust(15)))
							nextLine = inputFile.readline()
							if nextLine.startswith("//"): break
					else:
						outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name), \
						"\tNA".ljust(15)))
				currentLine = inputFile.readline()
				if currentLine == '': break # At the end of the file, python return an empty line
		inputFile.close()	
		outputFile.close()

if __name__ == '__main__':
	from convE2P2 import convE2P2
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg),
						help = "The E2P2 output file name")
	parser.add_argument("-o", "--out", 
						dest = "output", 
						type = str, 
						default = "[Input filename without its extension]_e2p2.conv",
						help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_e2p2.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_e2p2.conv"
	convE2P2(args.file, args.output).convert_file()
