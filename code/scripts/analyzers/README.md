INTRODUCTION
------------

The analyzer.py program makes it possible to perform an SNP analysis on the ec numbers predicted by the automatic annotation programs.
For each ec number, it is calculated a sensitivity, a precision and so on, of its detection by the automatic annotation programs.
The output results looks like something like:

EC      |     	PROG     |   	TP  |  	FP |   	FN  |  	SENS      |  	PREC   |     	F1SCORE    | 	N   |  	ZSCORE     
--------|----------------|----------|------|--------|-------------|------------|-------------------|--------|--------------
3.1.22.4|     	Blast2Go |    	1   |   4  |    0   |   1.0       |   	0.2    |      	0.33333333 |  	1   |   0.0        
4.6.1.16|     	E2P2     |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.58461282 
4.6.1.16|     	PRIAM    |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.58461282 
4.6.1.16|     	KOALA    |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.58461282 
4.6.1.16|     	KAAS     |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.58461282 
4.6.1.16|     	INTERPRO |    	1   |   0  |    2   |   0.33333333|   	1.0    |      	0.5        |  	6   |   -2.14358035 
4.6.1.16|     	Blast2Go |    	3   |   1  |    0   |   1.0       |   	0.75   |      	0.85714286 |  	6   |   -0.19487093 
4.6.1.17|     	Blast2Go |    	1   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	3   |   0.0        
4.6.1.17|     	KOALA    |    	1   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	3   |   0.0        
4.6.1.17|     	KAAS     |    	1   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	3   |   0.0        
4.6.1.12|     	E2P2     |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.4472136  
4.6.1.12|     	PRIAM    |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.4472136  
4.6.1.12|     	KOALA    |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.4472136  
4.6.1.12|     	KAAS     |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.4472136  
4.6.1.12|     	INTERPRO |    	1   |   0  |    2   |   0.33333333|   	1.0    |      	0.5        |  	6   |   -2.23606798
4.6.1.12|     	Blast2Go |    	3   |   0  |    0   |   1.0       |   	1.0    |      	1.0        |  	6   |   0.4472136  

sensitivity: measures the ability of a program to detect all ec numbers that appear true for an enzyme in the reference database.

precision: measures the ability of a program to only detect all ec numbers that appear true for an enzyme in the reference database.

DEPENDANCIES
------------

The analyzer.py script depends on lipmutils that is a LIPM perl library for sequence analysis. The used script in this library is one which
allows to retrieve the last ec numbers(__lipm_valid_ec_numbers.pl__). To install this library, follow the steps below:
1. Move yourself in the folder where the analyzer.py script is located
2. Recover the lipmutils library by running the command below: svn checkout --depth=immediates http://lipm-svn.toulouse.inra.fr/svn/lipmutils/trunk lipmutils
3. Move yourself in the lipmutils directory: cd lipmutils
4. svn update --set-depth=infinity bin corelib doc library
5. svn update --set-depth=exclude archive bdbh cfg demo javascript t web

Once this library is correctly installed, we can now test the analyzer.py script.


Command Line Arguments
----------------------
usage: analyzer.py [-h] [-std STD] [-plud PLUD] [-level LEVEL] [-out OUTPUT]
                   [conv_files [conv_files ...]]

positional arguments:

args		|helps
---------------|--------------------------------------------------------
  conv\_files  |The list of the converted programs files. Each file must be preceded by the name of the corresponded program (ex: E2P2 e2p2\_file.conv PRIAM priam\_file.conv) (default: None)

optional arguments:
       
args			|helps
--------------- |--------------------------------------------------------
-std STD        |  The converted standard file (default: None)
-path PATH      |  Path to the limputils directory (default: None)
-level LEVEL    |  level of ec classes (default: 4)
-out OUTPUT     |  The output file name (default: ecs_to_pick.tab)

-level: (optional) is a number from 1 to 4 specifying the number of digit on which the SNP analysis will be carried out. 
By default level is 4


--------------------
# predictor.py (PLEAP)

Perform annotation from multiple results from multiple annotation tools. The principle is as follows:
1. From the score threshold of f1-score received as parameter, we build a list of potentials ec numbers.
i.e the list of ec having a f1-score greater than or equal to the threshold.
2. An inventory of proteins is made, associating them with all the ecs assigned to it by the programs.
	
	For each protein and for each ec that has been associated with it:
		
		. If the ec exists in the list of potentials ecs  and the fraction between the number of programs that predicted 
		this ec  number and the number of programs wich are able of predicting it must be greater than or equal to the fraction 
		of threshold passed as a parameter (with option: -frac), then this protein is definitively associated with 
		this ec number.
		
		. Otherwise, the ec number will be removed. 
		
	The programs that have predicted an ec number and are not in the list of programs capable of detecting the this 
	ec number, and also they are not programs eliminated because of the insufficiency of their F1score, then their 
	prediction is retained. on the other hand for this precise case, the annotation is followed by a flag 
	indicating that the assignment is putative (the flag is: uncertified). The annotation having for flag certified, 
	means that all the conditions of the different thresholds are reached.
		
		
Command Line Arguments
----------------------
usage: predictor.py [-h] [-listec LISTEC] [-plud PLUD] [-score SCORE]
                    [-frac FRAC] [-level LEVEL] [-out OUT] [-oroc OROC]
                    [conv_files [conv_files ...]]

positional arguments:

args		|helps
---------------|--------------------------------------------------------
  conv\_files  |The list of the converted programs files. Each file must be preceded by the name of the corresponded program (ex: E2P2 e2p2\_file.conv PRIAM priam\_file.conv) (default: None)

optional arguments:

args			|helps
--------------- |---------------
-listec LISTEC  |The file containing the list of ec numbers of reference (default: None)
-plud PLUD      |to the limputils directory (default: None)
-score SCORE    |The f1 score value (default: 0.75)
-frac FRAC      |The minimum number of programs repporting an ec number to be considered as FP (default: 0.75)
-level LEVEL    |The level of ec classe (default: 4)
-out OUT        |The output file name (default: annot.out.txt)
-oroc OROC      |The output for roc displaying (default: out_roc.txt)

-oroc: Is optional argument: Use only when you want the file that will be used to build the Roc curve of the reported predictions. 


--------------------
# comparator.py

The __comparator.py__ script is used to produce files that will be used to build the barplot of comparison of pmrogrammes. 

The files produced with this script are the direct entries of the __comparator.R__ script.

__comparator.R__: is a script that can be found in the __Rscripts__ folder

Command Line Arguments
----------------------

usage: comparator.py [-h] [-pp PP] [-std STD] [-plud PLUD] [-level LEVEL]
                     [-out OUTPUT]
                     
args			|helps
--------------- |----------------------------------------------
-pp PP    		|The converted program file (default: None)
-std STD      	|The converted Swiss-Prot file (default: None)
-plud PLUD      |to the limputils directory (default: None)
-level LEVEL    |The level of ec classe (default: 4)
-out OUT        |The output file name (default: [First input filename without its extension].comp)

For each program in [PRIAM, E2P2, Blast2Go, ...], we produce a comparision file with this script and these files are used to build the image of barplot using __comparator.R__.


--------------------
# make_matrix.py
The script __make_matrix.py__ allows it to produce the file for computing the ROC curve of the predictions reported by PLEAP program.

Command Line Arguments
----------------------
usage: make_matrix.py [-h] [-pp PP] [-std STD] [-level LEVEL] [-out OUTPUT]

args			|helps
--------------- |----------------------------------------------
-pp PP    		|The converted program file (default: None)
-std STD      	|The converted Swiss-Prot file (default: None)
-level LEVEL    |The level of ec classe (default: 4)
-out OUT        |The output file name (default: matrix.txt)



