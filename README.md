# Calcul des OTU 

Vous trouverez la description complète du TP [ici](https://docs.google.com/document/d/1qWNqPZ9Ecd-yZ5Hpl6n2zd7ZGtHPjf3yaW1ulKRdWnk/edit?usp=sharing).

## Introduction

L’objectif de ce TP sera de calculer les OTU obtenues à partir d’un séquençage “mock”. Nous n’avons amplifié que les bactéries (et non les champignons). 8 espèces sont ainsi attendues.

Vous devrez développer un programme effectuant une dé-duplication en séquence complète (“dereplication full length”), une recherche des séquences chimériques et un regroupement basé sur un algorithme glouton (“Abundance Greedy Clustering”).  


## Installation des dépendances

Vous utiliserez les librairies nwalign3, pytest et pylint de Python:
```
pip3 install --user nwalign3 pytest pylint pytest-cov
```

## Utilisation

Vous devrez développer un programme python3 effectuant une dé-duplication en séquence complète (“dereplication full length”), une recherche des séquences chimériques et un regroupement basé sur un algorithme glouton (“Abundance Greedy Clustering”). Il prendra pour arguments:

 -i, -amplicon_file fichier contenant des séquences au format FASTA
 -s, -minseqlen Longueur minimum des séquences (optionnel - valeur par défaut 400)
 -m, -mincount Comptage minimum des séquences (optionnel - valeur par défaut 10)
 -c, -chunk_size Taille des partitions de séquence (optionnel - valeur par défaut 100)
 -k, -kmer_size Longueur des “kmer” (optionnel - valeur par défaut 8)
 -o, -output_file fichier de sortie avec les OTU au format FASTA

 ## Tests

Vous testerez vos fonctions à l’aide de la commande pytest --cov=agc à exécuter dans le dossier agc-tp/. En raison de cette contrainte, les noms des fonctions ne seront pas libre. Il sera donc impératif de respecter le nom des fonctions “imposées”, de même que leur caractéristique et paramètres. 
Vous vérifierez également la qualité syntaxique de votre programme en exécutant la commande: pylint agc.py

## Contact

En cas de questions, vous pouvez me contacter par email: amine.ghozlane[at]pasteur.fr

