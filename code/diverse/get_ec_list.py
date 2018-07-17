def write_list_ecNumbers(prog_file, lipmutils, level, outfile):
	"""
	Extract the list of ec number from a converted file of a program and write this list in a file.
	"""
	# try to transfer ec numbers with lipmutils program
	lipmutils_program = os.path.join(lipmutils, "bin", "lipm_valid_ec_numbers.pl")
	command = [lipmutils_program, "--infile", prog_file, "--outfile", ".lipm_out"]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	
	# try to write the ec numbers list in a file after being transferred
	ec_list0 = [x.split("\t")[1].strip(" \n\t\r") for x in open(".lipm_out", "r") if "NA" not in x]
	ec_list = list()
	for ec in ec_list0: 
		if "," in ec: 
			_list = ec.split(",")
			ec_list+= _list
		else: ec_list+= [ec]
	if os.path.exists(".lipm_out"): os.remove(".lipm_out")
	if os.path.exists(".lipm_out.err"): os.remove(".lipm_out.err")
	if 4 - int(level) == 0:
		ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == int(level)]
	else:
		ec_list = ['.'.join(ec.split('.')[:-4 + int(level)]) for ec in ec_list \
		if len([x for x in ec.split(".")[:-4 + int(level)] if x != "-"]) == int(level)]
	file = open(outfile, "w")
	for ec in set(ec_list): file.writelines(list("%s\n" % item for item in ec.split(",")))
	file.close()
	  	
if __name__ == '__main__':
	import warnings, os, subprocess
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	def is_valid_file(parser, arg):
		arg = os.path.abspath(arg)
		if not os.path.exists(arg):
			parser.error("The file %s does not exist!" % arg)
		else: return arg
	#
	# Checking arguments
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)	    
	parser.add_argument("--prog", dest = "prog",
						type = lambda arg: is_valid_file(parser, arg),
					   	help = "The converted program file")
	parser.add_argument("--path", dest = "path",
						help = "Path to the limputils directory")
	parser.add_argument("--level", dest = "level", default = "4", help = "level of ec classes")
	parser.add_argument("--out", dest = "output", default = "ec_numbers.list",
						help = "The output file name")
	args = parser.parse_args()
	if args.prog is None or args.level is None or args.path is None: 
		parser.print_help()
		exit(0)
	write_list_ecNumbers(args.prog, args.path, args.level, args.output)
