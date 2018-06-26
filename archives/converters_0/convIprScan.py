#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
import subprocess, commands, sys, os, time

class convIprScan:
	def __init__(self, inFile, outputFile):

		"""
		Check the type of inFile and convert it in a special format

		:Param inFile     : The input file name
		:Param outputFile : The name of the output file after the Converting
		:param typeFile   : The type of input file (output of InterproScan)
		
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
		Check if arg is a valid file which is of type InterproScan output.

		:Returns 	     : typeFile 
		:TypeOf typeFile : str
		"""
		typeFile , Error, output = "Incorrect", "", []
		# Try to check if input file is the output of Iprscan program
		try:
			output = subprocess.check_output("head -n 50 "+ self.inFile +" | cut -d'	' -f10"
					 +" | sort | uniq", shell = True)
			output = list(filter(None, output.splitlines()))
		except Exception as Error:	pass
		
		if Error == "" and len(output) == 1 and output[0] == "T":
			return "IPRSCAN"		
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
		
		# This part concerns the analysis of the contents of the files in order to bring out  
		# interesting informations. Thus, to optimize and  speed up the processing, 
		# bash commands have been invoked. 
		#
		if self.typeFile == "IPRSCAN":
			# Get the longest name among all protein. Knowing the longest name  
			# will allow us to better format the display in the output file 
			names = subprocess.check_output("awk '{print $1}' "+ self.inFile, shell = True)
			names = list(set(filter(None, names.splitlines())))
			max_name = max(len(name) for name in names)
			
			# Getting interpro2go and ec2go [intermediary files: will be deleted at the end]
			curl_ipr2go = "http://www.geneontology.org/external2go/interpro2go"
			curl_ec2go = "http://www.geneontology.org/external2go/ec2go"
			curl_list = [curl_ec2go, curl_ipr2go]
			start_time = time.time()
			for i in range(len(curl_list)):
				# progress bar
				i+= 1
				self.progress(i, len(curl_list), start_time, "Downloading intermediary files...")
				#
				subprocess.check_output("wget -cq -P ./tmp " + curl_list[i-1], shell = True)
			#
			# We check if the downloading of intermediary files has been done properly
			if os.path.exists("./tmp"):
				outputFile = open(self.outputFile, "w")
			else:
				text = "The files interpro2go and ec2go are not well downloaded please "
				text += "restart the program and make sure that you have a good connection"
				print ("%s\n", text)
				exit(0)
			
			i = 0 # set i for counting progress bar status
			start_time = time.time() # set time for the beginning of converting
			for ProtID in names:
				# progress bar
				i+= 1
				self.progress(i, len(names), start_time, "Converting file, Please wait...")
				#
				#In the input file, ProtID is associated with its corresponding 
				# InterproScan identifiers and scores.
				associated_IPRs = subprocess.check_output("grep '"+ ProtID + "' "+ self.inFile
					 				 	 +" | cut -d'	' -f 9,12 --output-delimiter='|'", 
					 				 	 shell = True)
				associated_IPRs = list(set(filter(None, associated_IPRs.splitlines())))
				
				# Each IPRID is associated with its score
				dic_associated_IPRs = {
					IPRID.split("|")[1]: IPRID.split("|")[0] for IPRID in associated_IPRs 
					if len(IPRID.split("|")) > 1 and IPRID.split("|")[0].strip() != "-"
				}
#					rest_associated_IPRs = [elt.split("|")[0] for elt in associated_IPRs 
#						if elt not in dic_associated_IPRs]
				#
				# In the file interpro2go, each IPRID is associated with its corresponding GOs
				dic_associated_GOs = {}
				for IPR in dic_associated_IPRs:
					ipr2go = subprocess.check_output("grep '"+ IPR + "' ./tmp/interpro2go"
					 				 	 +" | cut -d';' -f2", shell = True)
					ipr2go = list(set(filter(None, ipr2go.splitlines())))
					dic_ipr2go = {go:dic_associated_IPRs[IPR] for go in ipr2go}
					dic_associated_GOs.update(dic_ipr2go)
#					rest_associated_GOs = []
#					for IPR in rest_associated_IPRs:
#						ipr2go = subprocess.check_output("grep '"+ IPR + "' ./tmp/interpro2go"
#						 				 	 +" | cut -d';' -f2", shell = True)
#						ipr2go = list(set(filter(None, ipr2go.splitlines())))
#						rest_associated_GOs+= ipr2go
				#
				# In the file ec2go, each GO is associated with its corresponding ECs
				dic_associated_ECs = {}
				for GO in dic_associated_GOs:
					go2ec = subprocess.check_output("grep '"+ GO + "' ./tmp/ec2go"
					 				 	 +" | cut -d'>' -f1", shell = True)
					go2ec = list(set(filter(None, go2ec.splitlines())))
					dic_go2ec = {ec:dic_associated_GOs[GO] for ec in  go2ec}
					dic_associated_ECs.update(dic_go2ec)
				tuple_associated_ECs = [(k, v) for k, v in dic_associated_ECs.iteritems()]
#					rest_associated_ECs = []
#					for GO in rest_associated_GOs:
#						go2ec = subprocess.check_output("grep '"+ GO + "' ./tmp/ec2go"
#						 				 	 +" | cut -d'>' -f1", shell = True)
#						go2ec = list(set(filter(None, go2ec.splitlines())))
#						rest_associated_ECs+= go2ec
				
				# writing the ProtID with its associated ECs in the output file
				if len(tuple_associated_ECs) >= 1:
					outputFile.write( "{0} {1} {2}\n".format(ProtID.ljust(max_name + 1),
				 					tuple_associated_ECs[0][0].rstrip("\n").ljust(15), 
				 					tuple_associated_ECs[0][1].ljust(12)))
				else:
					outputFile.write( "{0}\n".format(ProtID.ljust(max_name + 1)))
				if len(tuple_associated_ECs) > 1:
				 	for (ec, score) in tuple_associated_ECs[1:]:
						outputFile.write( "{0} {1} {2}\n".format(" ".ljust(max_name + 1),
				 									 ec.rstrip("\n").ljust(15), 
				 									 score.ljust(12)))
#					if len(rest_associated_ECs) > 0:
#						for ec in rest_associated_ECs:
#							tmp_file.write( "{0} {1}\n".format(" ".ljust(max_name + 1),
#					 									 ec.rstrip("\n").ljust(12), 
#					 									 " ".ljust(12)))
#					print name, ":\n", tuple_associated_ECs, rest_associated_ECs					
			
		outputFile.close()
		# At the end of the program, the temporary files are deleted
		if os.path.exists("./tmp"):	print commands.getoutput("rm -r ./tmp")

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
	from convIprScan import convIprScan
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
	conv = convIprScan(args.file, args.output)
	print "\n{0} :{1}".format("Type Of File".ljust(18), conv.typeFile.ljust(10))
	conv.convert_file()
	print "{0} :{1} \n".format("The output file".ljust(18), conv.outputFile.ljust(10))
