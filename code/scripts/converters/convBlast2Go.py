#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from utils import *
class convBlast2Go:
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
		Check if arg is a valid file which is of type blast2go output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : boolean
		"""
		# Try to check if input file is the output of Blast2go program
		typeFile , Error, output = False, "", []
		try:
			inputFile = open(self.inFile, "r")
			output = [l.split("\t")[1][:2] for l in inputFile if len(l.split("\t")) > 1]
			output = list(set(output))
			inputFile.close()
		except Exception as Error:	pass
		if len(output) < 1 or len(output) > 2 or Error != "": return typeFile
		elif "EC" in output and "GO" in output: return True
		elif "EC" in output or "GO" in output:  return True
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
			# will allow us to better format the display in the output file
			blast2go = open(self.inFile, "r")
			names = list(set([line.split("\t")[0] for line in blast2go]))
			names = [line.split("|")[1] if "|" in line else line for line in names]
			max_name = max(len(name) for name in list(set(names)))
			blast2go.close()
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting
			i+= 1
			progress(i, len(names)+1, start_time, "Converting file, Please wait...")
			#
			# In order to better understand the parsing progress of this block of 
			# code, please inspect the structure of the Blast2go file. The procedure  
			# is as  follows: For each protein, we retrieve all lines where it appears 
			# associated  with an EC number and stock it in a dictionnary. Then,  we  
			# proceed to formatting the output text for great vision.
			inputFile = open(self.inFile, "r") 
			protIDs_with_ECs = {}
			for line in inputFile:
				line = line.split("\t")
				protID = line[0]
				if "|" in protID: protID = protID.split("|")[1]
				if "EC:" in line[1]:
					ecNumber = line[1].split("\t")[0].strip(' \n\t\r')
					if protID in protIDs_with_ECs: protIDs_with_ECs[protID].append(ecNumber)
					else: protIDs_with_ECs[protID] = [ecNumber]
			inputFile.close()
			
			outputFile = open(self.outputFile, "w")
			for protID in names:
				# progress bar
				i+= 1
				progress(i, len(names)+1, start_time, "Converting file, Please wait...")
				if protID in protIDs_with_ECs:
					for ec in protIDs_with_ECs[protID]:
						outputFile.write("{0} {1}\n".format(protID.strip().ljust(max_name), \
						"\t"+ ec.strip("EC:").rstrip().ljust(15)))	
				else: 
					outputFile.write( "{0} {1}\n".format(protID.strip().ljust(max_name), \
					"\tNA".ljust(15)))
		outputFile.close()

if __name__ == '__main__':
	from convBlast2Go import convBlast2Go
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "The Blast2go output file name")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = "[Input filename without its extension]_b2g.conv",
		                help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	default_output = "[Input filename without its extension]_b2g.conv"
	if args.output == default_output: args.output = get_output_name(args.file) + "_b2g.conv"
	convBlast2Go(args.file, args.output).convert_file()
