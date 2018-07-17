import sys, os, time, subprocess

def get_std_data(sp_file):
	uniprot_ref_file = open(sp_file, "r")
	ref_dic = {}
	for line in uniprot_ref_file:
		protID = line.split()[0].strip(" \n\t\r")
		ec = line.split()[1].strip(" \n\t\r")
		if protID in ref_dic: ref_dic[protID].append(ec)
		else: ref_dic[protID] = [ec]
	uniprot_ref_file.close()
	return ref_dic
	
def get_prog_data(prog_file, lipmutils):
	# try to transfer ec numbers with lipmutils program
	lipmutils_program = os.path.join(lipmutils, "bin", "lipm_valid_ec_numbers.pl")
	time_stamp = time.time()
	outfile = str(prog_file) + "."+ str(time_stamp) +".transf"
	command = [lipmutils_program, "--infile", prog_file, "--outfile", outfile]
	subprocess.call(command, stdout = open('/dev/null', 'w'), stderr = subprocess.STDOUT)
	prog_pred_file = open(outfile, "r")
	transferred_dic = {}
	for line in prog_pred_file:
		protID = line.split("\t")[0].strip(" \n\t\r")
		if "|" in protID: protID = line.split("|")[1].strip(" \n\t\r")
		ec_list = line.split("\t")[1].strip(" \n\t\r").split(",")
		if protID in transferred_dic:
			transferred_dic[protID]+= ec_list
		else: transferred_dic[protID] = ec_list
	prog_pred_file.close()
	os.remove(outfile)
	os.remove(outfile + ".err")
	
	prog_pred_file = open(prog_file, "r")
	prog_pred_dic = {}
	for line in prog_pred_file:
		protID = line.split("\t")[0].strip(" \n\t\r")
		if "|" in protID: protID = line.split("|")[1].strip(" \n\t\r")
		if protID in transferred_dic:
			if protID in prog_pred_dic:
				ec_list = prog_pred_dic[protID]
				ec_list+= transferred_dic[protID]
				prog_pred_dic[protID] = list(set(ec_list))
			else: prog_pred_dic[protID] = transferred_dic[protID]
		else: prog_pred_dic[protID] = ['NA']
	prog_pred_file.close()
	return prog_pred_dic
	
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

def are_valid_file(parser, arg):
	"""
	Check if arg are valid files that already exists on the files system.

	Parameters
	----------
	parser : argparse object
	arg    : str

	Returns
	-------
	arg
	"""
	Type, i = "", 0
	while i < len(arg.split(",")):
		if i+1 == len(arg.split(",")): Type += is_valid_file(parser, arg.split(",")[i])
		else: Type += is_valid_file(parser, arg.split(",")[i]) + ","
		i+= 1
	return Type

def get_output_name(filename):
		"""
		Returns the file name without its extension.
		"""
		infile_extension = filename.split(".")[-1]
		return filename[:len(filename) - len(infile_extension) - 1]

def progress(count, total, start_time, status = 'Running'):
	"""
	This function allows us to track the progress of program execution in 
	the console
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
