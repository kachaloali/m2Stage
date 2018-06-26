INTRODUCTION
------------

The analyzer.py program makes it possible to perform an SNP analysis on the ec numbers predicted by the automatic annotation programs.
For each ec number, it is calculated a sensitivity, a specificity and an accuracy of its detection by the automatic annotation programs.
The output results looks like something like:

EC              | sens            | spec            | prec
--------------- | --------------- | --------------- | ------------
1.17.4.1        | 1.0             | 0.99984114      | 0.83333333
3.1.4.46        | 1.0             | 1.0             | 1.0         
1.8.7.2         | 1.0             | 1.0             | 1.0         
1.8.7.1         | 1.0             | 1.0             | 1.0         
6.3.2.17        | 1.0             | 0.99984119      | 0.75        
1.4.1.14        | 1.0             | 1.0             | 1.0         
1.1.1.234       | 0.38461538      | 1.0             | 1.0         
1.1.1.237       | 1.0             | 1.0             | 1.0         
1.1.1.236       | 1.0             | 0.99952358      | 0.5         
6.3.2.19        | 0.0             | 0.99222222      | 0.0         

Sensitivity: measures the ability of a program to detect all ec numbers that appear true for an enzyme in the reference database.

Acurancy: measures the ability of a program to only detect all ec numbers that appear true for an enzyme in the reference database.

Specificity: measures the ability of a program to not detect the ec numbers that do not appear true for an enzyme in the reference database.

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

python analyzer.py --prog inputfile1 --sp inputfile2 --level level

--prog | --sp: (required) Specifies the input files: __inputfile1__ is the converted program output file and __inputfile2__ is the output 
file produce by __convUniprot.py__ script.

--level: (optional) is a number from 1 to 4 specifying the number of digit on which the SNP analysis will be carried out. 
By default level is 4

Example:

	python analyzer.py --prog examples/tracheophyta_pri.text --sp examples/tracheophyta_sp.text
	
Check the results in examples directory:

	tracheophyta_pri.ec.list: is the list of the ec numbers predicted by the program
	
	tracheophyta_pri.snp.4.txt: is the output SNP analysis file that is descrived above
	
	
-------------------------
The predictor.py program
-------------------------

Perform annotation from multiple results from multiple annotation tools. The principle is as follows:
1. From the score threshold of f1-score received as parameter, we build a list of potentials ec numbers.
i.e the list of ec having a f1-score greater than or equal to the threshold.

2. The potentials ec numbers predicted by all programs are retained as valid ecs.
3. The partial intersections of at least 3 programs are retained as valid ecs: intersections containing 
only potentials ec numbers.
4. The partial intersections of a program with the reference (swissprot) are retained as valid ecs: 
intersections containing only potentials ec numbers.
5. An inventory of proteins is made, associating them with all the ecs assigned to it by the programs.
	For each protein and for each ec that has been associated with it:
		. If the ec exists in the list of valid ecs, then this protein is definitively associated with 
		this ec number.
		. Otherwise, the ec number will be removed.
		
		
Command Line Arguments
----------------------

python predictor.py [--prog PROG] [--sp SP] [--plud PLUD] [--pas PAS] [--level LEVEL] [--score SCORE] [--mi MI] [--out OUT]

args			|helps
--------------- |---------------
--prog PROG     |List of the converted programs files separated by comma (default: None)
--sp SP         |The converted swiss-prot file (default: None)
--plud PLUD     | Path to the limputils directory (default: None)
--pas PAS       |Path to the analyzer.py script (default: None)
--level LEVEL   |The level of ec classe (default: 4)
--score SCORE   |The f1 score value (default: 0.5)
--mi MI         |The minimum number of programs constituting an intersection to be considered (default: 3)
--out OUT       |The output file name (default: annotation.out.txt)

