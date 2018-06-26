INTRODUCTION
------------

Runners scripts are used to run automated enzymatic annotation programs with each of its specific parameters.
Each program has its own runner associate with the best parameters to produce more good results.
The main role of these runners is to be able to automate the execution of annotation programs at the same time.

	runE2P2.py:
		usage: runE2P2.py [-h] [-f FILE] [-p PATH] [-o OUTPUT]
	
		  -h, --help            show this help message and exit
		  -f FILE, --file FILE  The input fasta file name (default: None)	: (is required)
		  -p PATH, --path PATH  The path to e2p2 program (default: None)	: (is required)
		  -o OUTPUT, --out OUTPUT
				                The output file name (default: [Input fasta filename without its extension].pf)
	
	runIprScan.py:
		The interproscan is installed on LIPM nod cluster. To execute this program directly on the command line 
		we use the command below:
		/usr/local/bioinfo/interproscan/interproscan.sh --pathways --outfile OUTPUTFILE --iprlookup --input PROTEOME --goterms --formats TSV --disable-precalc --clusterrunid UNIDAUPIF --mode cluster
		However, there is a particularity when we want to submit the work via qsub command. We must create a bash 
		script containing the command above then execute the created bash script via qsub otherwise interproscan will not work.
		Suppose we have the script runIprScan.sh containing the command so we can run this script via qsub with the 
		appropriate arguments as follows: qsub -shell yes -S /bin/bash -q all.q@@bigmem -cwd runIprScan.sh
		
		The description above is exactly what is in the file runIprScan.py
			
			usage: runIprScan.py [-h] [-f FILE] [-p PATH] [-o OUTPUT]
				  -h, --help            show this help message and exit
				  -f FILE, --fasta FILE
								        The input fasta file name (default: None)			: (is required)
				  -p PATH, --path PATH  The path to interproscan program (default: None)	: (is required)
				  -o OUTPUT, --out OUTPUT
								        The output file name (default: [Input fasta filename without its extension].tsv)
	
	runPriam.py:
		usage: runPriam.py [-h] [-f FF] [-pp PP] [-pr PR] [-pt TH] [-mp MP] [-mo MO]

		  -h, --help  show this help message and exit
		  -f FF       The input fasta file (default: None)			: (is required)
		  -pp PP      The path to Priam program (default: None)		: (is required)
		  -pr PR      The path to Priam Releases (default: None)	: (is required)
		  -pt TH      threshold (default: 0.6)
		  -mp MP      Minimal length proportion of a profile that must be matched to
				      consider it (default: 20)
		  -mo MO      Maximum overlap length between the matches of two profiles
				      (default: 40)
		
		The default parameters such as pt, mp, mo are the best parameters that I determined for the annotation
		
