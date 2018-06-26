## How To Run converters programs:
[![Dependencies Python2](https://img.shields.io/badge/Python-2.7-green.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

	I explain for a single program "convPRIAM.py", but it is the same approach for all the following converters 
	programs:
		convBlast2Go.py
		convIprScan.py
		convKAAS.py
		convE2P2.py
	The only file includes in each of Theses above programs is : utils.py
	
	Run "convPRIAM.py" using input file ("-f") and output file name ("-o"). If no output file name is specified, 
	the default output is the input file name without its extension followed by this signature "_pri.conv". Below 
	is the general command:
    
	python convPRIAM.py -f inputfile -o outputfile

	If no output file name is specified for others converters programs, theirs signatures ares:
		"_b2g.conv" 	for convBlast2Go.py
		"_iprscan.conv"	for convIprScan.py
		"_kaas.conv"	for convKAAS.py
		"_e2p2.conv"	for convE2P2.py
	
	I will give an example of using "convPRIAM.py" program. Take for example the file tracheophyta.txt which is the 
	result of annotation from the priam program. To convert this file, use the following command line:
	
	python convPRIAM.py -f examples/tracheophyta.text -o examples/tracheophyta_out1.text
	
	You should find tracheophyta_out1.txt file in the examples directory which is the converted file produce by "convPRIAM.py"
	
	The converter "convUniprot.py" is a bit special because it takes as input two files. it also includes utils.py
	Below is the general command to run "convUniprot.py" program:
	
	python convUniprot.py -f inputfile.fasta -t inputfile.tab -o outputfile
	
	The fasta and tab files are both from the swiss-prot database
	
	Bellow, an example of using "convUniprot.py" program:
	consider the tracheophyta.Fasta and tracheophyta.tab files in the examples directory. To test this program, 
	use the following command line:
	
	python convUniprot.py -f examples/tracheophyta.Fasta -t examples/tracheophyta.tab -o examples/tracheophyta_out2.text
	
	You should find tracheophyta_out2.text file in the examples directory which is the converted file produce by "convUniprot.py"
	
The onverters programs are used to convert annotation results from multiple annotation tools into one format that 
will allow a simple manipulation.
