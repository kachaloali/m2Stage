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
	'''
	Reads a Expasy Enzyme .dat file and writes a tab separated file where the first column is 
	EC number, the second column is the reaction description, the third column is the associated 
	uniprot ids separated by '|', and the fourth column indicates whether the reactions described 
	by this EC have been transferred to other EC numbers.
	'''
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
			out[record['ID']]['uniprot'] = '|'.join([x[0] for x in record['DR']])
			out[record['ID']]['description'] = record['DE']
			out[record['ID']]['transferred'] = False
	for id in transferred:
		out[id] = dict()
		out[id]['uniprot'] = '|'.join([out[x]['uniprot'] for x in transferred[id]])
		out[id]['description'] = 'Transferred entry: ' + ' '.join(transferred[id])
		out[id]['transferred'] = True
	df = pd.DataFrame.from_dict(out, orient = 'index')
	df.index.name = 'EC'
	
	# write all data in a enzyme.csv file
	df.to_csv(output_name, sep = '\t')
	
	# ignore EC numbers with no uniprot ids associated
	df.dropna(subset = ['uniprot'], inplace = True)
	
	# ignore EC numbers that are obsolete due to transfer 
	df = df[df.transferred == False]
	
	# construct a dictionnary from dataframe
	mydic = df.to_dict()
	
	enzyme_protIDS = [mydic["uniprot"][ec].split("|") for ec in mydic["uniprot"].keys()]
	enzyme_protIDS = list(set(reduce(lambda x, y: x + y, enzyme_protIDS)))
	enzyme_protIDS = [elt.strip(" \n\t\r") for elt in enzyme_protIDS if elt != ""]
	dic_ecs = dict()
	dic_ecs["1."] = set()
	dic_ecs["2."] = set()
	dic_ecs["3."] = set()
	dic_ecs["4."] = set()
	dic_ecs["5."] = set()
	dic_ecs["6."] = set()

	csvfile = open('uniprot-reviewed%3Ayes.tab', 'r')
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
		uniprot_iDs = mydic["uniprot"][ec]
		protIDs = uniprot_iDs.strip(" \n\t\r").split("|")
		protIDs = [elt for elt in protIDs if elt != ""]
		frag_seqs = list()
		short_seqs = list()
		for seqID in protIDs:
			if "Fragment" in dic_sp[seqID]['name']:
				 frag_seqs.append(seqID)
			if int(dic_sp[seqID]['length']) < 50:
				short_seqs.append(seqID)
		protIDs=[e for e in protIDs if not e in frag_seqs and not e in short_seqs]
		if ec.startswith("1."):
			dic_ecs["1."].update(protIDs)
		elif ec.startswith("2."):
			dic_ecs["2."].update(protIDs)
		elif ec.startswith("3."):
			dic_ecs["3."].update(protIDs)
		elif ec.startswith("4."):
			dic_ecs["4."].update(protIDs)
		elif ec.startswith("5."):
			dic_ecs["5."].update(protIDs)
		elif ec.startswith("6."):
			dic_ecs["6."].update(protIDs)
		non_valids_enzyme.update(frag_seqs)
		non_valids_enzyme.update(short_seqs)
	
	# step 2
	# for the uniqueness, those enzymes that occur in two or more classes were excluded
	for ec in ["2.", "3.", "4.", "5.", "6."]:
		dic_ecs["1."] = dic_ecs["1."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["1."].intersection(dic_ecs[ec]))
	for ec in ["1.", "3.", "4.", "5.", "6."]:
		dic_ecs["2."] = dic_ecs["2."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["2."].intersection(dic_ecs[ec]))
	for ec in ["2.", "1.", "4.", "5.", "6."]:
		dic_ecs["3."] = dic_ecs["3."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["3."].intersection(dic_ecs[ec]))
	for ec in ["2.", "3.", "1.", "5.", "6."]:
		dic_ecs["4."] = dic_ecs["4."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["4."].intersection(dic_ecs[ec]))
	for ec in ["2.", "3.", "4.", "1.", "6."]:
		dic_ecs["5."] = dic_ecs["5."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["5."].intersection(dic_ecs[ec]))
	for ec in ["2.", "3.", "4.", "5.", "1."]:
		dic_ecs["6."] = dic_ecs["6."].difference(dic_ecs[ec])
		non_valids_enzyme.update(dic_ecs["6."].intersection(dic_ecs[ec]))
	
	# step 3: 
	# to reduce the homology bias, a redundancy cutoff was operated by cd-hit program to winnow
	# those sequences which have >=40% sequence identity to any other in a same functional class

	#
	# Downlod and constructing fasta file of the Main functional classes
	def split_sequence(seq, l):
		new_seq = ""
		if len(seq) > l:
			new_seq = seq[:l]
			k = l
			while k + l < len(seq):
				new_seq+= "\n"+str(seq[k:k+l])
				k+= l
			new_seq+= "\n" + str(seq[k:])
			return new_seq + "\n"
		else: return seq + "\n"
	def run_process(list_seqs, filename):
		# @nested function
		session = requests.Session()
		outfile = open(filename, "a")
		for seqID in list_seqs:
			#handle = ExPASy.get_sprot_raw(seqID.strip(" \n\r\t"))
			#record = SeqIO.read(handle, "swiss")
			#SeqIO.write(record, outfile, "fasta")
			req = "http://wwwdev.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession="+str(seqID)
			res = session.get(req, headers = {'User-Agent' : 'application/XML Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11',
			"content-type":"text"})
			# parse the returned XML
			uniprot = ET.fromstring(res.text)
			for isoform in uniprot.getchildren():
				# get the sequence
				iso_sequence = isoform.find('{http://uniprot.org/uniprot}sequence')
				# get the accession number
				iso_accession = isoform.find('{http://uniprot.org/uniprot}accession')
				outfile.write(">"+str(iso_accession.text)+"\n")
				outfile.write(split_sequence(str(iso_sequence.text), 60))		   		
		outfile.close()
	def create_process(list_seqs, filename):
		# @nested function:
		p = Process(target = run_process, args = (list_seqs, filename,))
		p.start()
		return p
	
	#ec1 = create_process(dic_ecs["1."], "knnDataset/ec_1.*.faa")
	#ec2 = create_process(dic_ecs["2."], "knnDataset/ec_2.*.faa")
	#ec3 = create_process(dic_ecs["3."], "knnDataset/ec_3.*.faa")
	#ec4 = create_process(dic_ecs["4."], "knnDataset/ec_4.*.faa")
	#ec5 = create_process(dic_ecs["5."], "knnDataset/ec_5.*.faa")
	#ec6 = create_process(dic_ecs["6."], "knnDataset/ec_6.*.faa")
	
	
	#===================================================o===========================================
	# Selection rules for the subclasses: same screening procedures	than the Main functional classes 
	#===================================================o===========================================
	# step 1
	# those enzymes whose sequences were annotated with 'fragment' were excluded
	# those enzymes whose sequences had less than 50 amino acids were excluded
	dic_subclasses = dict()
	for ec in mydic["description"].keys():
		uniprot_iDs = mydic["uniprot"][ec]
		protIDs = uniprot_iDs.strip(" \n\t\r").split("|")
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
	
	# making fasta files
#	list_process = list()
#	for ec in dic_subclasses:
#		process = create_process(dic_subclasses[ec], os.path.join(output_dir, str(ec)+".faa"))
#		list_process.append(process)
#	for i in range(len(list_process)):
#		while list_process[i].is_alive(): time.sleep(60)
	for ec in dic_subclasses:
		file = open(os.path.join(output_dir, str(ec)+".ids.list"), 'w')
		for seqID in dic_subclasses[ec]: file.write("%s\n" % seqID)
	  	file.close()
	for ec in dic_subclasses:
		batch = os.path.join(output_dir, str(ec)+".ids.list")
		fasta  = os.path.join(output_dir, str(ec)+".faa")
		print commands.getoutput("blastdbcmd -db "+ database +" -entry_batch "+ batch +" > "+ fasta)
		#outfile = open(os.path.join(output_dir, str(ec)+".faa"), "a")
		
		#for seqID in dic_subclasses[ec]:
#			handle = ExPASy.get_sprot_raw(seqID.strip(" \n\r\t"))
#			record = SeqIO.read(handle, "swiss")
#			SeqIO.write(record, outfile, "fasta") 
#			req = "http://wwwdev.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession="+str(seqID)
#			#res = requests.get(req, headers = {'User-Agent' : 'application/XML Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11'})
#			print commands.getoutput("wget -cq -P "+ output_dir +" '" + req + "'")
#			tree = ET.parse(os.path.join(output_dir, os.path.basename(req)))
#			uniprot = tree.getroot()
#			# parse the returned XML
#			#uniprot = ET.fromstring(res.text)
#			for isoform in uniprot.getchildren():
#				# get the sequence
#				iso_sequence = isoform.find('{http://uniprot.org/uniprot}sequence')
#				# get the accession number
#				iso_accession = isoform.find('{http://uniprot.org/uniprot}accession')
#				outfile.write(">"+str(iso_accession.text)+"\n")
#				outfile.write(split_sequence(str(iso_sequence.text), 60))
#			os.remove(os.path.join(output_dir, os.path.basename(req)))		
		#outfile.close()
		#process = create_process(dic_subclasses[ec], os.path.join(output_dir, str(ec)+".faa"))
		#while process.is_alive():
		#	time.sleep(60)
	
	# step 3: 
	# to reduce the homology bias, a redundancy cutoff was operated by cd-hit program to winnow
	# those sequences which have >=40% sequence identity to any other in a same functional class
#	for ec in dic_subclasses:
#		print commands.getoutput(cdhit +" -i "+os.path.join(output_dir, str(ec)+".faa")
#			+" -d 0 -o "+ os.path.join(output_dir, str(ec) +".cdhit.faa")
#			+" -c 0.4 -n 2  -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000 > "
#			+ os.path.join(output_dir, str(ec) +".out"))

	print "\tFINISHED"
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
	
##	handle = ExPASy.get_sprot_raw("P84071")
##	seq_record = SeqIO.read(handle, "swiss")
##	print(seq_record.id)
##	print(seq_record.name)
##	print(seq_record.description)
##	print(repr(seq_record.seq))
##	print("Length %d" % len(seq_record))
##	print(seq_record.annotations["keywords"])
##	print(df)
