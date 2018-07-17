# 1. get_taxon.py

The script __get_taxon.py__ allows to download the taxon files of a giving organism taxon names from swiss-prot database.
To use it you must give in parameters the taxon correct name and the format to be download. The format type is [fasta or tab]

	usage: get_taxon.py [-h] [--taxon TAXON] [--format FORMAT]

	  -h, --help       show this help message and exit
	  --taxon TAXON    taxon name (default: None)						: (is requiered)
	  --format FORMAT  The output format [fasta or tab] (default: None)	: (is requiered)
	
Example:
Test to download tracheophyta taxonomy. So use the __get_taxon.py__ like this:
>>```bash
>> python get_taxon.py --taxon "tracheophyta" --format "fasta"
>> python get_taxon.py --taxon "tracheophyta" --format "tab"
```
You must have the files tracheophyta.fasta and tracheophyta.tab

# 2. split_fasta_file.py
The script __split_fasta_file.py__ is used to split large fasta file into slightly smaller files.

	usage: split_fasta_file.py [-h] [-f FILE] [-n NUMBER]

	  -h, --help            show this help message and exit
	  -f FILE, --fasta FILE The input fasta file name (default: None)					: (is requiered)
	  -n NUMBER             The number of sequences from which to cut (default: None)	: (is requiered)
exmple:
>>```bash
>> python split_fasta_file.py -f inputFastaFile.fasta -n 3000
```	
This command produce several spliting files from __inputFastaFile.fasta__

# 3. get_ec_list.py
The script __get_ec_list.py__ is used to build the list of ECs numbers of a converted results file of a programme.

	usage: get_ec_list.py [-h] [--prog PROG] [--path PATH] [--level LEVEL]
                      [--out OUTPUT]

		optional arguments:
		  -h, --help     show this help message and exit
		  --prog PROG    The converted program file (default: None)		  
		  --path PATH    Path to the limputils directory (default: None)		  
		  --level LEVEL  level of ec classes (default: 4)		  
		  --out OUTPUT   The output file name (default: ec_numbers.list)

As you can see it, the script requires the lipmutils library. EC numbers are converted before making the list.

# 4. Installing PRIAM locally on linux system
[![Dependencies Blast](https://img.shields.io/badge/Blast-2-green.svg)](https://ubuntu.pkgs.org/16.04/ubuntu-universe-amd64/blast2_2.2.26.20120620-10_amd64.deb.html)
[![Dependencies Releases](https://img.shields.io/badge/Releases-v15-blue.svg)](http://priam.prabi.fr/REL_MAR15/index_mar15.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

##### Step 1: Install Blast suite
>> Download the package __blast2_2.2.26.20120620-10_amd64.deb__ and install it with the command line:
```bash
>> sudo gdebi blast2_2.2.26.20120620-10_amd64.deb
```
>> Or just run the following command to install it from the remote apt repository:
```bash
>> sudo apt-get install blast2
```

##### Step 2: Download complete releases of PRIAM
>> Download only releases of PRIAM constructed after August 2010.

##### Step 3: Download PRIAM program PRIAM_search.jar
>> Download PRIAM program at [PRIAM_Search utility](http://priam.prabi.fr/REL_MAR15/index_mar15.html)

##### Requirements
>> java: is the only one program requiered to execute __PRIAM_search.jar__ after the dependencies listed above.

##### Test example on sunflower proteome : **HanXRQr1.0-20151230-EGN-r1.2.prot_without_annotation.faa**
Once you have established all the things listed above, create a directory that you will name testPriam. 

Firstly decompress the archive that contain releases of PRIAM and move it in testPriam directory that you have just created.
Move also the package __PRIAM_search.jar__ into the folder testPriam. 

Secondly, create a folder that you will name TOURNESOL and put there the sunflower proteome fasta file (Be careful to give the 
correct file name for the __pathToFastaFile__ variable in the following python script).

Finally, still in testPriam directory, create a file that you will name testPriam.py and insert the following python script.
```python
	# -*- coding: utf-8 -*-
	import os
	import commands

	pathToPriamReleases = "PRIAM_MAR15" # directory that contain PRIAM releases files
	pathToFastaFile = "TOURNESOL/HanXRQr1.0-20151230-EGN-r1.2.prot_without_annotation.faa"

	# directory that will contain intermediates, temporary files and Results
	pathToResults = "RESULTS" 


	if os.path.exists(pathToResults): print commands.getoutput("rm -r "+ pathToResults)
	print commands.getoutput(
		"java -jar PRIAM_search.jar -i "+ pathToFastaFile + " -p "+pathToPriamReleases + 
		" -od "+ pathToResults + " -n testPriam -pt 0.5 -mo -1 -mp 70 -cc T -cg T -e T")

```

If all is well, run the testPriam.py script using the following command and check your results in the __RESULTS/RESULTS__ directory.
>>```bash
>> python testPriam.py
```

# 5. How to create profiles database

First of all, we will need the PSI-BLAST tool from [HERE](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) 
>> Or just run the following command to install it from the remote apt repository:
```bash
>> sudo apt-get install ncbi-blast+
```

We will also need a database for exemple the [nr database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) dataset: So download for exemple the numbered files from __nr.00.tar.gz__ to __nr.05.tar.gz__ and extract all of them somewhere into a same folder called __nr/__  or create this __nr/__ folder and move in to it then run the __get_nr_db.py__ script into the __nr/__ folder like this:
```bash
>> python get_nr_db.py --number 5
```
Remarque: it take long time ! 

Let's say we have a fasta file called __tracheophyta.Fasta__ containing some protein sequences and  located in the examples directory

PSI-BLAST needs each sequence as separate input files, so we need to create a fasta file for each one. We will do this at the same time as running PSI-BLAST. We will assume we extracted the [nr database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) dataset into a directory called __nr/__.

Now to run PSI-BLAST, We will move in to the __nr/__ directory. This ensures PSI-BLAST can find the dataset. Edit the file called __nr.pal__ by removing all unsused databases names in the __DBLIST__. In our case, beacause we are using 5 datasets of nr, this list must be: DBLIST "nr.00" "nr.01""nr.02" "nr.03" "nr.04" "nr.05"

From the __nr/__ directory run the __get_pssm.py__ python script using these parameters below:

>>```bash
	python get_pssm.py --query examples/tracheophyta.Fasta --db nr --nb_it 2 --evalue 0.0001
```

The script above outputs each sequence of __tracheophyta.Fasta__ into a file called train/train_N.fasta, where N is a number that is incremented for each protein. we will have to make sure the directory __train__ where we are putting these files in exists. There will also need to be directories called __out__ and __pssm__. These will be filled with PSI-BLAST outputs and the pssm files respectively.

Once you have successfully completed your profile directory called __pssm__, you can now build your database containing only the list of the profiles that you just built. To do so, you would start by making a simple text file called __pssm.pn__ in to __pssm__ directory. This file contain the names of the associated PSSM profiles (the order doesn't matter): In our case this file must contain:
```text
	pssm_1.chk
	pssm_2.chk
	pssm_3.chk
	pssm_4.chk
	pssm_5.chk
```
Then run the __makeprofiledb__ tool to build a database with the following command:
>>```bash
>> makeprofiledb -in pssm.pn -title pssm_db -out pssm_db -dbtype 'rps'
```

The above assumed you are in the same directory as the __pssm.pn__ file and all the __*.chk__ files. After running, this does creates the nine files: __pssm.aux__, __pssm.freq__, __pssm.loo__, __pssm.phr__, __pssm.pin__, __pssm.psd__, __pssm.psi__, __pssm.psq__ and __pssm.rps__ which together make up the database. 

# 6. How to run RPS-BLAST

As PSI-BLAST, RPS-BLAST needs also each sequence as separate input files so we are going to prepare a fasta file called __tracheophyta1.Fasta__ which is in examples directory.
Then run the the following command (Make sure that you are in the __pssm__ directory!)
```bash
	rpsblast -i examples/tracheophyta1.Fasta -d pssm -e 0.00001 -o rpsblast_result.txt
```
You should find your results in the __rpsblast_result.txt__ file. It remains to parser these RPS-BLAST results.
