
#!/usr/bin/env python
# <!-- -*- coding: utf-8 -*- -->
from multiprocessing import Process
import subprocess, os, pandas as pd
from Bio.ExPASy import Enzyme
import requests, sys
import warnings
warnings.filterwarnings("ignore")
# -*- coding: utf-8 -*-

try: 
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

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
import pprint
#pprint.pprint(mydic["sequenceID"])
#exit(0)
    
f = open("ecs_to_pick_l4.tab", "r").readlines()[1:]
progs = ["PRIAM", "E2P2", "INTERPRO", "KAAS", "KOALA", "Blast2Go"]
i = 5
seuilSens = 0.75
seuilPrec = 0.75
fw1 = open(progs[i]+".1", "w")
fw2 = open(progs[i]+".2", "w")
fw1.write("{0} {1} {2} {3} {4} {5} {6}\n".format("ECs".ljust(12),"\tSENS".ljust(12),"\tPREC".ljust(12),"\tF1-Score".ljust(12), "\tnbSEQS".ljust(12), "\tProfils".ljust(12), "\tKEEG-Hiostory".ljust(12)))
fw2.write("{0} {1} {2} {3} {4} {5} {6}\n".format("ECs".ljust(12),"\tSENS".ljust(12),"\tPREC".ljust(12),"\tF1-Score".ljust(12), "\tnbSEQS".ljust(12), "\tProfils".ljust(12), "\tKEEG-Hiostory".ljust(12)))
for l in f:
	l = l.split("\t")
	ec = l[0].strip(" \n\t\r")
	prog = l[1].strip(" \n\t\r")
	sens = float(l[5].strip(" \n\t\r"))
	prec = float(l[6].strip(" \n\t\r"))
	f1score = float(l[7].strip(" \n\t\r"))
	
	if prog == progs[i]:
		if sens < seuilSens and prec < seuilPrec:
			#print ec, sens, prec, f1score
			requestURL = "http://www.genome.jp/dbget-bin/www_bget?ec:"+ str(ec)
			r = requests.get(requestURL, headers={ "Accept" : "application/HTML"})
			soup = BeautifulSoup(r.content, 'lxml')
			res = soup.find_all('td', class_ = "td21")
			indice = 0
			if res is not None:
				for j in range(len(res)):
					if ">EC " in str(res[j]):
						indice = j
						break
			if indice != 0:
				history = res[indice].find('div').text.strip("\n")
			else:
				history = "None"
			
			######################
			reqPriam = "http://priam.prabi.fr/cgi-bin/PRIAM_profiles_CurrentRelease.pl?EC="+str(ec)
			r = requests.get(reqPriam, headers={ "Accept" : "application/HTML"})
			soup = BeautifulSoup(r.content, 'lxml')
			res = soup.find_all('a')
			profils = [elt.text for elt in res if "psg" in str(elt)]
			#print profils
			AA = len(mydic["sequenceID"][str(ec)].split("|"))
			lP = len(profils)
			fw1.write("{0} {1} {2} {3} {4} {5} {6}\n".format(str(ec).ljust(12),"\t"+str(sens).ljust(12),"\t"+str(prec).ljust(12),"\t"+str(f1score).ljust(12), "\t"+str(AA).ljust(12), "\t"+str(lP).ljust(12), "\t"+str(history).ljust(12)))
			print ec, ":\t", sens, ":\t", prec, ":\t", f1score, ":\t", len(mydic["sequenceID"][str(ec)]) , ":\t", len(profils), ":\t", history
		else:
			#print ec, sens, prec, f1score
			requestURL = "http://www.genome.jp/dbget-bin/www_bget?ec:"+ str(ec)
			r = requests.get(requestURL, headers={ "Accept" : "application/HTML"})
			soup = BeautifulSoup(r.content, 'lxml')
			res = soup.find_all('td', class_ = "td21")
			indice = 0
			if res is not None:
				for j in range(len(res)):
					if ">EC " in str(res[j]):
						indice = j
						break
			if indice != 0:
				history = res[indice].find('div').text.strip("\n")
			else:
				history = "None"
			
			######################
			reqPriam = "http://priam.prabi.fr/cgi-bin/PRIAM_profiles_CurrentRelease.pl?EC="+str(ec)
			r = requests.get(reqPriam, headers={ "Accept" : "application/HTML"})
			soup = BeautifulSoup(r.content, 'lxml')
			res = soup.find_all('a')
			profils = [elt.text for elt in res if "psg" in str(elt)]
			#print profils
			AA = len(mydic["sequenceID"][str(ec)].split("|"))
			lP = len(profils)
			fw2.write("{0} {1} {2} {3} {4} {5} {6}\n".format(str(ec).ljust(12),"\t"+str(sens).ljust(12),"\t"+str(prec).ljust(12),"\t"+str(f1score).ljust(12), "\t"+str(AA).ljust(12), "\t"+str(lP).ljust(12), "\t"+str(history).ljust(12)))






