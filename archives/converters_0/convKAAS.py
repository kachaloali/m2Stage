#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
import subprocess, commands, sys, os, time

class convKAAS:
	def __init__(self, KoFile, Ko2ecFile, outputFile):

		"""
		Check the type of inFile and convert it in a special format

		:Param inFile_ko     : The input file name (kaas.ko)
		:Param inFile_ko2ec  : The input file name (ko2ec)
		:Param outputFile 	 : The name of the output file after the Converting
		:param typeFile   	 : The type of input file (output of kaas)
		
		:TypeOf inFile_ko    : str
		:TypeOf inFile_ko2ec : str
		:TypeOf outputFile	 : str
		:TypeOf typeFile     : str 
		
		:Seealso:: get_type_file(inFile)
		:Seealso:: convert_file(inFile, typeFile)
		"""
		self.inFile_ko	  = str(KoFile)
		self.inFile_ko2ec = str(Ko2ecFile)
		self.outputFile   = str(outputFile)
		self.typeFile     = self.get_type_file()
		
	
	def get_type_file(self):
		"""
		Check if arg is a valid file which is of type kaas output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : str
		"""
		typeFile , Error, output = "Incorrect", "", []
		# Try to check if input file is the output of kaas program
		try:
			output = subprocess.check_output("awk '{print $2}' "+ self.inFile_ko +
					 "| grep -m 50 '^K' | cut -c 1 | sort | uniq", shell = True)
			output = list(filter(None, output.splitlines()))
		except Exception as Error:	pass
		if Error == "":
			if len(output) == 1 and output[0] == "K": return "KAAS"
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
			print("\nThe file %s is not in the correct format!" % self.inFile_ko)
			exit(0)
		outputFile = open(self.outputFile, "w")
		kaas_ko = open(self.inFile_ko, "r")
		# This part concerns the analysis of the contents of the files in order to bring out interesting 
		# informations. Thus, to optimize and  speed up the processing, bash commands have been invoked. 
		#
		if self.typeFile == "KAAS":
			# Get the longest name among all protein. Knowing the longest name will allow us to better 
			# format displaying in the output file 
			names = subprocess.check_output("awk -F '\t' '{print $1}' "+ self.inFile_ko, shell = True)
			names = list(set(filter(None, names.splitlines())))
			max_name = max(len(name) for name in names)
			#
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting	
			for line in kaas_ko:
				# progress bar
				i+= 1
				self.progress(i, len(names), start_time, "Converting file, Please wait...")
				#
				line = line.rstrip("\n").split("\t")
				if len(line) > 1:
					
					Error = ""
					nameProtein = line[0]
					ko_associated = line[1]
					try:
						list_associated_ECs = subprocess.check_output("grep '"+ ko_associated +"' "+  
											  self.inFile_ko2ec, shell = True)
						list_associated_ECs = list(set(filter(None, list_associated_ECs.splitlines()))) 
					except Exception as Error: pass
					if Error == "":
						list_ECs = list_associated_ECs[0].split("\t")[1].strip("[").strip("]").split(" ")
						outputFile.write("{0} {1}\n".format(nameProtein.ljust(max_name + 1), 
														(list_ECs[0]).ljust(15)))
						for ec in list_ECs[1:]:
							outputFile.write("{0} {1}\n".format(" ".ljust(max_name + 1), 
															   ("EC:" + ec).ljust(15)))
				else: 
					outputFile.write("{0}\n".format(nameProtein.ljust(max_name + 1)))
		kaas_ko.close()	
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
	from convKAAS import convKAAS
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
		infile_extension = arg.koFile.split(".")[-1]
		return arg.koFile[:len(arg.koFile) - len(infile_extension) - 1] + ".conv"
	
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
	                    formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-k", "--ko", 
						dest = "koFile", 
						type = lambda arg: is_valid_file(parser, arg), 
						help = "ko file name")
	parser.add_argument("-e", "--ko2ec",
	                    dest = "ko2ecFile",
	                    type = str,
	                    help = "ko2ec file name")
	parser.add_argument("-o", "--out",
	                    dest = "output",
	                    type = str,
	                    default = get_output_name(parser),
	                    help = "The output file name")
	args = parser.parse_args()
	if args.koFile is None or args.ko2ecFile is None: 
			parser.print_help()
			exit(0)
	#
	#
	conv = convKAAS(args.koFile, args.ko2ecFile, args.output)
	print "\n{0} :{1}".format("Type Of File".ljust(18), conv.typeFile.ljust(10))
	conv.convert_file()
	print "{0} :{1} \n".format("The output file".ljust(18), conv.outputFile.ljust(10))
