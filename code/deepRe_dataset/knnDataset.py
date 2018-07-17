#!/usr/bin/python2
# <!-- -*- coding: utf-8 -*- -->
from multiprocessing import Process
import subprocess, os, pandas as pd
from Bio.ExPASy import Enzyme
from Bio import SwissProt
from Bio import ExPASy
from Bio import SeqIO
import numpy as np, sys
import requests, commands,csv,time
import xml.etree.ElementTree as ET

def get_knndataset(cdhit, output_dir, database): 
	#Reads a Expasy Enzyme .dat file and writes a numpy data frame where the first column is 
	#EC number, the second column is the reaction description, the third column is the associated 
	#sequenceID ids separated by '|', and the fourth column indicates whether the reactions described 
	#by this EC have been transferred to other EC numbers.
	if not os.path.exists(os.path.join("database", "enzyme", "enzyme.dat")):
		curl_enzyme = os.path.join("ftp://ftp.expasy.org", "databases","enzyme", "enzyme.dat")
		subprocess.check_output("wget -cq -P database/enzyme " + curl_enzyme, shell = True)
	if not os.path.exists(os.path.join("database", "enzyme", "enzyme.dat")): 
		print ("%s\n", "Missing enzyme database!")
		exit(0)
	input_name = os.path.join("database", "enzyme", "enzyme.dat")
	output_name = os.path.join("database", "enzyme", "enzyme.tsv")
	records = Enzyme.parse(open(input_name))
	out = dict() # dict of dicts, first key: EC number, second key: field
	transferred = dict() #dict of lists
	for record in records:
		if 'Transferred entry:' in record['DE']:
			record['DE'] = record['DE'].rstrip('.')
			record['DE'] = record['DE'].replace('Transferred entry:',' ')
			record['DE'] = record['DE'].replace(',',' ')
			record['DE'] = record['DE'].replace('and',' ')
			point_to = record['DE'].split()
			transferred[record['ID']] = point_to
		else:
			out[record['ID']] = dict()
			out[record['ID']]['sequenceID'] = '|'.join([x[0] for x in record['DR']])
			out[record['ID']]['description'] = record['DE']
			out[record['ID']]['transferred'] = False
	for id in transferred:
		out[id] = dict()
		out[id]['sequenceID'] = '|'.join([out[x]['sequenceID'] for x in transferred[id]])
		out[id]['description'] = 'Transferred entry: ' + ' '.join(transferred[id])
		out[id]['transferred'] = True
	df = pd.DataFrame.from_dict(out, orient = 'index')
	df.index.name = 'EC'
	# write all data in a enzyme.csv file
	df.to_csv(output_name, sep = '\t')
	# ignore EC numbers with no sequenceID ids associated
	df.dropna(subset = ['sequenceID'], inplace = True)
	# ignore EC numbers that are obsolete due to transfer 
	df = df[df.transferred == False]
	
	
	# The numpy data frame is converted to a python dictionnary 
	mydic = df.to_dict()
	enzyme_protIDS = [mydic["sequenceID"][ec].split("|") for ec in mydic["sequenceID"].keys()]
	enzyme_protIDS = list(set(reduce(lambda x, y: x + y, enzyme_protIDS)))
	enzyme_protIDS = [elt.strip(" \n\t\r") for elt in enzyme_protIDS if elt != ""]
	dic_ecs = dict()
	dic_ecs["1"] = set()
	dic_ecs["2"] = set()
	dic_ecs["3"] = set()
	dic_ecs["4"] = set()
	dic_ecs["5"] = set()
	dic_ecs["6"] = set()
	if not os.path.exists(os.path.join("database", "uniprot", "sp.tab")):
		url = os.path.join("http://www.uniprot.org","uniprot","?query=reviewed:yes&format=tab")
		subprocess.check_output("wget -cq -P database/uniprot '"+ url + "'", shell = True)
		subprocess.check_output("mv database/uniprot/*=tab database/uniprot/sp.tab",shell=True)
	if not os.path.exists(os.path.join("database", "uniprot", "sp.tab")): 
		print ("%s\n", "Missing uniprot database!")
		exit(0)
	csvfile = os.path.join("database", "uniprot", "sp.tab")
	readCSV = csv.reader(csvfile, delimiter = '\t')
	non_valids_enzyme = set()
	dic_sp = dict()
	for row in readCSV: 
		if row[0] != "Entry":
			seqID = row[0]
			seqName = row[3]
			seqLength = row[6]
			dic_sp[seqID] = dict()
			dic_sp[seqID]['name'] = seqName
			dic_sp[seqID]['length'] = seqLength

	#===================================================o========================================
	# Selection rules for the Main functional classes							 
	#===================================================o========================================
	# step 1
	# those enzymes whose sequences were annotated with ‘‘fragment’’ were excluded
	# those enzymes whose sequences had less than 50 amino acids were excluded
	for ec in mydic["description"].keys():
		sequenceID_iDs = mydic["sequenceID"][ec]
		protIDs = sequenceID_iDs.strip(" \n\t\r").split("|")
		protIDs = [elt for elt in protIDs if elt != ""]
		frag_seqs = list()
		short_seqs = list()
		for seqID in protIDs:
			if "Fragment" in dic_sp[seqID]['name']:
				 frag_seqs.append(seqID)
			if int(dic_sp[seqID]['length']) < 50:
				short_seqs.append(seqID)
		protIDs=[e for e in protIDs if not e in frag_seqs and not e in short_seqs]
		if ec.startswith("1"):
			dic_ecs["1"].update(protIDs)
		elif ec.startswith("2"):
			dic_ecs["2"].update(protIDs)
		elif ec.startswith("3"):
			dic_ecs["3"].update(protIDs)
		elif ec.startswith("4"):
			dic_ecs["4"].update(protIDs)
		elif ec.startswith("5"):
			dic_ecs["5"].update(protIDs)
		elif ec.startswith("6"):
			dic_ecs["6"].update(protIDs)
		non_valids_enzyme.update(frag_seqs)
		non_valids_enzyme.update(short_seqs)
	
	# step 2
	# for the uniqueness, those enzymes that occur in two or more classes were excluded
	for ec in ["2", "3", "4", "5", "6"]:
		dic_ecs["1"] = dic_ecs["1"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["1"].intersection(dic_ecs[ec]))
	for ec in ["1", "3", "4", "5", "6"]:
		dic_ecs["2"] = dic_ecs["2"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["2"].intersection(dic_ecs[ec]))
	for ec in ["2", "1", "4", "5", "6"]:
		dic_ecs["3"] = dic_ecs["3"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["3"].intersection(dic_ecs[ec]))
	for ec in ["2", "3", "1", "5", "6"]:
		dic_ecs["4"] = dic_ecs["4"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["4"].intersection(dic_ecs[ec]))
	for ec in ["2", "3", "4", "1", "6"]:
		dic_ecs["5"] = dic_ecs["5"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["5"].intersection(dic_ecs[ec]))
	for ec in ["2", "3", "4", "5", "1"]:
		dic_ecs["6"] = dic_ecs["6"].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["6"].intersection(dic_ecs[ec]))
	
	# these following two functions are internal and allow to create processes to parallel the fasta 
	# files downloading and their passage to the cd-hit program
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
	
	# step 3: 
	# to reduce the homology bias, a redundancy cutoff was operated by cd-hit program to winnow
	# those sequences which have >=40% sequence identity to any other in a same functional class
	# making fasta files for the six main classes
	for ec in dic_ecs:
		create_process(dic_ecs[ec], str(ec), output_dir, database, cdhit)

	
	#===================================================o===========================================
	# Selection rules for the subclasses: same screening procedures	than the Main functional classes 
	#===================================================o===========================================
	# step 1
	# those enzymes whose sequences were annotated with 'fragment' were excluded
	# those enzymes whose sequences had less than 50 amino acids were excluded
	dic_subclasses = dict()
	for ec in mydic["description"].keys():
		sequenceID_iDs = mydic["sequenceID"][ec]
		protIDs = sequenceID_iDs.strip(" \n\t\r").split("|")
		protIDs = [elt for elt in protIDs if elt != ""]
		frag_seqs = list()
		short_seqs = list()
		for seqID in protIDs:
			if "Fragment" in dic_sp[seqID]['name']:
				 frag_seqs.append(seqID)
			if int(dic_sp[seqID]['length']) < 50:
				short_seqs.append(seqID)
		protIDs=[e for e in protIDs if not e in frag_seqs and not e in short_seqs]
		list_ec_digits = [x for x in ec.split(".") if x != "-"]
		if len(list_ec_digits) >= 2:
			ec_on_l2 = '.'.join(list_ec_digits[:2])
			if ec_on_l2 in dic_subclasses: dic_subclasses[ec_on_l2].update(protIDs)
			else: dic_subclasses[ec_on_l2] = set(protIDs)
	
	# step 2
	# for the uniqueness, those enzymes that occur in two or more classes were excluded
	for ec1 in dic_subclasses.keys():
		for ec2 in dic_subclasses.keys():
			if ec1 != ec2: dic_subclasses[ec1] = dic_subclasses[ec1].difference(dic_subclasses[ec2])
	#print(len(dic_subclasses))
	excluded_ecs = list()
	for ec in dic_subclasses:
		if len(dic_subclasses[ec]) < 10: excluded_ecs.append(ec)
	dic_subclasses = {k: v for k, v in dic_subclasses.items() if k not in excluded_ecs}
	
	# step 3: 
	# to reduce the homology bias, a redundancy cutoff was operated by cd-hit program to winnow
	# those sequences which have >=40% sequence identity to any other in a same functional class
	for ec in dic_subclasses:
		# making fasta files for the subclasses: after retrieving associated fasta file and 
		# reducing redundancy with cd-hit program
		create_process(dic_subclasses[ec], str(ec), output_dir, database, cdhit)
	#
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
	get_knndataset(args.p, args.d, args.db)
