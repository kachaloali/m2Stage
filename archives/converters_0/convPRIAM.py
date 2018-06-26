#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
import subprocess, commands, sys, os, time

class convPRIAM:
	def __init__(self, inFile, outputFile):

		"""
		Check the type of inFile and convert it in a special format

		:Param inFile     : The input file name
		:Param outputFile : The name of the output file after the Converting
		:param typeFile   : The type of input file (output of priam)
		
		:TypeOf inFile    : str
		:TypeOf outputFile: str
		:TypeOf typeFile  : str 
		
		:Seealso:: get_type_file(inFile)
		:Seealso:: convert_file(inFile, typeFile)
		"""
		self.inFile	  = str(inFile)
		self.outputFile = str(outputFile)
		self.typeFile = self.get_type_file()
		
	
	def get_type_file(self):
		"""
		Check if arg is a valid file which is of type priam output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : str
		"""
		typeFile , Error, output = "Incorrect", "", []
		# Try to check if input file is the output of Priam program
		try:
			test = subprocess.check_output("grep -m 1 '>' "+ self.inFile, shell = True)
			Error = ""
		except Exception as Error:	pass
		if Error == "":
			output = subprocess.check_output("head -n 50 "+ self.inFile +
					 " | grep -v '>' | cut -d ' ' -f8 | sort | uniq", shell = True)
			output = list(filter(None, output.splitlines()))
		if len(output) == 2: return "PRIAM"
		else: Error = ""
		return typeFile
	
	def convert_file(self):
		"""
		Convert the input file in the format that we have described. 
		The format is to be checked in the output file of this function
	
		:Returns 		: no returns
		:Actions 		: creating a file
		"""	
		if self.typeFile == "Incorrect":
			print("\nThe file %s is not in the correct format!" % self.inFile)
			exit(0)
		outputFile = open(self.outputFile, "w")
		inputFile = open(self.inFile, "r")
		
		# This part concerns the analysis of the contents of the files in order to bring out  
		# interesting informations. Thus, to optimize and  speed up the processing, 
		# bash commands have been invoked. 
		#
		if self.typeFile == "PRIAM":
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format the display in the output file 
			names = subprocess.check_output("awk -F '\t' '{if ($1 ~ /^>/) print $1}' "+ self.inFile,
										   shell = True)
			names = list(set(filter(None, names.splitlines())))
			max_name = max(len(name) for name in names)
			#
			# In order to better understand the parsing progress of this block of code, please 
			# inspect the structure of the Priam file. The procedure is as follows: as long as we 
			# are not at the end of the file and there exist a protein ID, it is associated with 
			# its EC numbers.
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting
			currentLine = inputFile.readline()
			while currentLine.startswith(">"):
				# progress bar
				i+= 1
				self.progress(i, len(names), start_time, "Converting file, Please wait...")
				#
				nameProtein = currentLine[1:].rstrip()
				nextLine = inputFile.readline()
				# To have a good text formatting in output file, the spacing for the first 
				# column is based on the length of the longest protein ID.
				list_ = nextLine.split("\t")
				outputFile.write( "{0} {1} {2}\n".format(nameProtein.ljust(max_name + 1),
													    ("EC:"+ list_[0]).ljust(15), 
													    (list_[1]).ljust(12)))
				nextLine = inputFile.readline()
				while nextLine.rstrip() != "":
					list_ = nextLine.split("\t")
					outputFile.write( "{0} {1} {2}\n".format("".ljust(max_name + 1),
													   		("EC:"+ list_[0]).ljust(15), 
													   		(list_[1]).ljust(12)))
					nextLine = inputFile.readline()
				currentLine = inputFile.readline()
		inputFile.close()	
		outputFile.close()
	
	def progress(self, count, total, start_time, status = 'Running'):
		"""
		This function allows us to track the progress of program execution in the console
		The execution bar, the percentage and the elapsed time are visible
		"""
		bar_len = 60
		filled_len = int(round(bar_len * count / float(total)))

		percents = round(100.0 * count / float(total), 1)
		if count == total: bar = '=' * filled_len +'-' * (bar_len - filled_len)
		else: bar = '=' * filled_len + '>' +'-' * (bar_len - filled_len - 1)
		m, s = divmod(time.time() - start_time, 60)
		h, m = divmod(m, 60)

		if count == total:
			time.sleep(0.5)
			sys.stdout.write("\033[K")
			sys.stdout.write('|%s|%s%s|%d:%02d:%02d| %s\r\n' % 
							(bar, percents, '%', h, m, s, status +" OK"))
		else:
			sys.stdout.write('|%s|%s%s|%d:%02d:%02d| %s\r' % 
							(bar, percents, '%', h, m, s, status))
		sys.stdout.flush()



if __name__ == '__main__':

	import os
	from convPRIAM import convPRIAM
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	
	def is_valid_file(parser, arg):
		"""
		Check if arg is a valid file that already exists on the file system.

		Parameters
		----------
		parser : argparse object
		arg    : str

		Returns
		-------
		arg
		"""
		arg = os.path.abspath(arg)
		if not os.path.exists(arg):
			parser.error("The file %s does not exist!" % arg)
		else: return arg

	def get_output_name(parser):
		"""
		Returns the file name by replacing its old extension with the ".conv" one.
		"""
		arg = parser.parse_args()
		infile_extension = arg.file.split(".")[-1]
		return arg.file[:len(arg.file) - len(infile_extension) - 1] + ".conv"
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--file", 
						dest = "file", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "File name")
	parser.add_argument("-o", "--out",
		                dest = "output",
		                type = str,
		                default = get_output_name(parser),
		                help = "The output file name")
	args = parser.parse_args()
	if args.file is None: 
			parser.print_help()
			exit(0)
	#
	#
	conv = convPRIAM(args.file, args.output)
	print "\n{0} :{1}".format("Type Of File".ljust(18), conv.typeFile.ljust(10))
	conv.convert_file()
	print "{0} :{1} \n".format("The output file".ljust(18), conv.outputFile.ljust(10))
