## knnDataset.py: 
C'est un script établit pour télécharger le 1er jeu de données de l'article DEEPre: sequence-based enzyme EC number prediction by deep learning

Le script s'exécute sur le serveur nod ou symbiose(A cause de la base de données uniprot qui est indexée sur nod et symbiose). Il faut donc le copier à destination de nod ou symbiose.

Le script depend du programme cdhit.

En étant sur le serveur nod par exemple, télécharger et compiler le programme cd-hit:
 
	1. git clone https://github.com/weizhongli/cdhit.git
	
	2. cd cdhit
	
	3. make
Le programme cdhit se trouve alors dans le dossier: cdhit

usage: knnDataset.py [-h] [-p P] [-d D] [-db DB]

optional arguments:
  -h, --help  show this help message and exit
  
  -p P        The path to the cdhit program (default: None)
  
  -d D        The directory of the output files (default: None)
  
  -db DB      The indexed database (default: None)
 
Pour construire le jeu de donner, on exécute le script de cette manière: 

python knnDataset.py -p cdhit/cd-hit -d knnDataset -db /db/biomaj/blastdb/uniprot

Avec __cdhit__ le dossier où se trouve le programme cd-hit, __knnDataset__ un dossier préalablement crée pour sauvegarder les données de sortie et __/db/biomaj/blastdb/uniprot__ le chemin vers la base de données indexée pour extraire les fichiers fasta.

Sur nod et symbiose, le chemin vers la base uniprot est le même: __/db/biomaj/blastdb/uniprot__


## newDataset.py
Il est établit pour télécharger le 2e jeu de données de l'article DEEPre: sequence-based enzyme EC number prediction by deep learning

Il s'exécute aussi sur le serveur nod ou symbiose(A cause de la base de données uniprot qui est indexée sur nod et symbiose).

usage: newDataset.py [-h] [-p P] [-d D] [-db DB]

optional arguments:

  -h, --help  show this help message and exit
  
  -p P        The path to the cdhit program (default: None)
  
  -d D        The directory of the output files (default: None)
  
  -db DB      The indexed database (default: None)

Pour construire le jeu de donner, on exécute le script de cette manière: 

python newDataset.py -p cdhit/cd-hit -d newDataset -db /db/biomaj/blastdb/uniprot

Avec __cdhit__ le dossier où se trouve le programme cd-hit, __newDataset__ un dossier préalablement crée pour sauvegarder les données de sortie et __/db/biomaj/blastdb/uniprot__ le chemin vers la base de données indexée pour extraire les fichiers fasta.

		
