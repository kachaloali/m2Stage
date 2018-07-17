#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from multiprocessing import Process
from Bio import ExPASy
from Bio import SeqIO
import csv, subprocess, time
def main(cdhit, output_dir, database):
	"""
	This script is the complete procedure to retrieve the dataset constructed in the article
	DEEPre: sequence-based enzyme EC number prediction by deep learning.
	"""
	if not os.path.exists(os.path.join("database", "uniprot", "sp.tab")):
		url = os.path.join("http://www.uniprot.org","uniprot","?query=reviewed:yes&format=tab")
		subprocess.check_output("wget -cq -P database/uniprot '"+ url + "'", shell = True)
		subprocess.check_output("mv database/uniprot/*=tab database/uniprot/sp.tab",shell=True)
	if not os.path.exists(os.path.join("database", "uniprot", "sp.tab")): 
		print ("%s\n", "Missing uniprot database!")
		exit(0)
	csvfile = os.path.join("database", "uniprot", "sp.tab")
	readCSV = csv.reader(csvfile, delimiter = '\t')
	dic_sp = dict()
	for row in readCSV: 
		if row[0] != "Entry":
			seqID = row[0]
			seqName = row[3]
			seqLength = row[6]
			dic_sp[seqID] = dict()
			dic_sp[seqID]['name'] = seqName
			dic_sp[seqID]['length'] = seqLength
	sp_protIDs = dic_sp.keys()
	sp_protIDs = [elt.strip(" \n\t\r") for elt in sp_protIDs if elt != ""]
	sp_protIDs = [elt for elt in sp_protIDs if int(dic_sp[elt]['length']) >= 50]


	# i
	def is_annotated(description):
		if 	 "EC 1." in description: return True
		elif "EC 2." in description: return True
		elif "EC 3." in description: return True
		elif "EC 4." in description: return True
		elif "EC 5." in description: return True
		elif "EC 6." in description: return True
		else: return False
	sp_enzyme_protIDs = [seq for seq in sp_protIDs if is_annotated(dic_sp[seq]['name'])]
	sp_non_enzyme_protIDs=[s for s in sp_protIDs if not is_annotated(dic_sp[s]['name'])]


	# ii
	def get_complete_ecs(description):
		ec_list = [element.split(")")[0] for element in description.split("(EC ")[1:]]
		ec_list = [ec for ec in ec_list if len([x for x in ec.split(".") if x != "-"]) == 4]
		return [ec for ec in ec_list if not "n" in ec.split(".")[3]]
	dic_enzy_sp = {iD: get_complete_ecs(dic_sp[iD]['name']) for iD in sp_enzyme_protIDs}
	dic_enzy_sp = {iD: value for iD, value in dic_enzy_sp.items() if len(value) == 1 }


	# iii
	dic_enzy_sp = {k:v for k, v in dic_enzy_sp.items() if not "Fragment" in dic_sp[k]['name']}
	dic_enzy_sp = {iD:val for iD,val in dic_enzy_sp.items() if int(dic_sp[iD]['length']) >=50}
	dic_enzy_sp = {iD: v for iD,v in dic_enzy_sp.items() if int(dic_sp[iD]['length']) <= 5000}
	sp_non_enzyme_protIDs=[k for k in sp_non_enzyme_protIDs if int(dic_sp[k]['length']) >= 50]
	sp_non_enzyme_protIDs=[k for k in sp_non_enzyme_protIDs if int(dic_sp[k]['length'])<=5000]


	# Downlod and constructing fasta file of the dataset
	def run_process(list_seqs, filename, output_dir, database, cdhit):
		# @nested function
		# Downlod and constructing fasta file of suclasses
		file = open(os.path.join(output_dir, filename + ".ids.list"), 'w')
		for seqID in list_seqs: file.write("%s\n" % seqID)
	  	file.close()
		fasta = os.path.join(output_dir, filename + ".faa")
		batch = os.path.join(output_dir, filename +".ids.list")
		print commands.getoutput("blastdbcmd -db "+ database +" -entry_batch "+ batch +" > "+ fasta)
		os.remove(batch)
		# run cdhit program
		cdhitout = os.path.join(output_dir, filename +".cdhit.faa")
		cdhitverbose = os.path.join(output_dir, filename +".out")
		print commands.getoutput(cdhit +" -i "+ fasta +" -d 0 -o " + cdhitout
			+" -c 0.4 -n 2  -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000 > " + cdhitverbose)  	
	def create_process(list_seqs, filename, output_dir, database, cdhit):
		# @nested function:
		p = Process(target = run_process, args = (list_seqs, filename, output_dir, database, cdhit,))
		p.start()
		return p
	# iv
	create_process(dic_enzy_sp, "new_enzyme", output_dir, database, cdhit)
	
	# v
	create_process(sp_non_enzyme_protIDs, "new_non_enzyme", output_dir, database, cdhit)

if __name__ == '__main__':
	import commands, shutil, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-p", dest = "p", help = "The path to the cdhit program")	        
	parser.add_argument("-d", dest = "d", help = "The directory of the output files")
	parser.add_argument("-db", dest = "db", help = "The indexed database")
	args = parser.parse_args()
	if args.p is None or args.d is None or args.db is None: 
			parser.print_help()
			exit(0)
	main(args.p, args.d, args.db)
