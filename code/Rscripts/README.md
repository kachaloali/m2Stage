##### comparator.R

This Rscript create a comparision graphic (barplots of sensitivities precision and f1-score) between the automatic enzyme annotation programs. The Rscript work in command line. If command line invoked, pass external arguments to R when calling this Rscript which needs, as an input, the comparision files names produce by __comparator.py__ script (__comparator.py__ is in analyzer directory). The result provides by this program is a file containing the comparison image.

Usage: comparator.R [options]

Options:
	--progs=CHARACTER
		The programs files name (Files produce with comparator.py script)

	--out=CHARACTER
		output file name [default= myIMG.png]

	-h, --help
		Show this help message and exit

comparator.R --progs Prog1\ file1\ Prog2\ file2 --out myImage

As you see it, you give the name of the program, you escape the space, and then you give the name of the program's comprison file. Everything is in couple (that's to say: each name of the program is accompanied with its comparision file). After that, you give the name of the output image without any extension.

##### vennDiagram.R
This Rscript compute a venn diagram graphic from a list of ec predicted by the automated annotation programs 
The Rscript work as an interactive mode. If command line invoked, pass external arguments to R when calling 
this Rscript which needs, as an input, the list of ec files names produce by __get_ec_list.py__ (__get_ec_list.py__ is in the directory named __diverse__) script. 
The result provides by this program is a file containing the venn diagram image.

Usage: vennDiagram.R [options]

Options:

	--priam=CHARACTER
		Priam file name

	--e2p2=CHARACTER
		E2P2 file name

	--b2g=CHARACTER
		Blast2go file name

	--kaas=CHARACTER
		KAAS file name

	--koala=CHARACTER
		KOALA file name

	--iprscan=CHARACTER
		InterproScan file name

	--std=CHARACTER
		The swiss prot standard file name

	--level=CHARACTER
		The level of ec number classe

	--pleap=CHARACTER
		The new pipline file name

	--out=CHARACTER
		output file name [default= myIMG.png]

	-h, --help
		Show this help message and exit


##### heatmap.R
It is not executed on command line this script. Open it with Rstudio, learn about the required files and execute it.
This script builds a Heatmap of the 3 essentials informations containing in the evaluation file.


##### rocCurve.R
It is not also executed on command line. It allows to build ROC curves individually or several arranged in the same image.
It is executed using the Rstudio environment. It will be necessary to indicate the files required for its execution (description of the files in the self script)

##### multi_comp.py
Allows to build several barPlots and arrange them in the same graphic. This script is not executed on the command line. You have to look in the script to see the files required for its execution.
