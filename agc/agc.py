#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Pierre Cote"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Pierre Cote"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Pierre Cote"
__email__ = "cotepierre@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

import re
def read_fasta(amplicon_file, minseqlen):
    """
    Generator that yield full length sequences from a fasta file
    """
    with gzip.open(amplicon_file, 'rb') as fasta:
        content = fasta.read().decode('utf-8')
    p = re.compile(r">\w+\n([ATCG\n]+)")
    completes = p.findall(content)
    completes = list(map(lambda x: x.replace('\n', ''), completes))
    for seq in completes:
        if len(seq) >= minseqlen:
            yield seq

from collections import defaultdict
def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Generator that yield sequence present more than mincount times
    """
    dict = defaultdict(int)
    for seq in read_fasta(amplicon_file, minseqlen):
        dict[seq] += 1
    for key, value in sorted(dict.items(), key = lambda x: x[1], reverse=True):
        if value >= mincount:
            yield key, value


def get_chunks(sequence, chunk_size):
    """
    Return sequence cut in 4 chunks of size chunk_size.
    Raise ValueError if length of sequence is smaller than 4*chunk_size.
    """
    l = len(sequence)
    if 4*chunk_size > l:
        raise ValueError
    return [sequence[i*chunk_size:i*chunk_size+chunk_size] for i in range(0, 4)]

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    Generator of kmer of size kmer_size from sequence
    """
    for i in range(len(sequence) - kmer_size +1):
        yield sequence[i:i+kmer_size]

def detect_chimera(perc_identity_matrix):
    """
    Si l’écart type moyen des pourcentages d’identité est supérieur à 5 et
    que 2 segments minimum de notre séquence montrent une similarité différente
    à un des deux parents, nous identifierons cette séquence comme chimérique
    """
    stddev = 0
    sim1 = {}
    sim2 = {}
    for identities in perc_identity_matrix:
        stddev += statistics.stdev(identities)
        sim1.add(identities[0])
        sim2.add(identities[1])
    if len(sim1) >= 2 or len(sim2) >= 2:
        stddev_mean = stddev/len(perc_identity_matrix)
        if stddev_mean > 5:
            return True
    return False


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    Populate and return kmer_dict iteratively.
    Kmer_dict is a dictionary of key (= kmer of size kmer_size from sequence) and value (=list of sequence id that contain the kmer key).
    """
    for i in cut_kmer(sequence, kmer_size):
        if i not in kmer_dict:
            kmer_dict[i] = list()
        kmer_dict[i].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, seq, kmer_size):
    """
    Return a list of the 8 most similar sequence id from kmer_dict, compared with target sequence seq.
    The most similar sequences compared to seq have the most kmer in common.
    """
    res=""
    for kmer in cut_kmer(seq, kmer_size):
        if kmer in kmer_dict:
            res+="".join(map(str, kmer_dict[kmer]))
    count = Counter(res).most_common(8)
    return list(map(lambda x: int(x[0]), count))


def get_identity(alignment_list):
    """
    Compare and score alignment of 2 sequences, score is nb_matching_amino_acids / length_aligned_sequences.
    """
    a, b = alignment_list
    count = match = 0
    for i in range(len(a)):
        if a[i] == b[i]:
            match += 1
        count += 1
    return match/count*100

def calcul_identity_matrix(chunks, parents, chunk_size, non_chimeres):
    """
    Compute identity matrix between current sequences (represented in chunks) and parent non chimere sequences (list of 2 raw sequences)
    """
    perc_id_matrix = [[] for i in range(len(chunks))]
    for parent in parents:
        parent_chunks = get_chunks(non_chimeres[parent], chunk_size)
        for i, chunk in enumerate(chunks):
            alignement = nw.global_align(chunk, parent_chunks[i])
            identity = get_identity(alignement)
            perc_id_matrix[i].append(identity)
    return perc_id_matrix

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Generator that detect chimera in fasta amplicon_file, yields non chimera sequences with their presence counts
    Iteratively detect chimera over dereplicated sequences from fasta amplicon file
    For every sequences present more than mincount times in amplicon file:
        cut sequences in chunk
        search mates in kmer_dict for each chunk
        find parents from chunk mates list (parents = mates that appears in every chunk)
        if we detect more than 2 parents, maybe their is a chimera
        detect it and continue or populate kmer_dict and yield non chimera sequence with its count
    """
    kmer_dict = {}
    non_chimeres = []
    id_seq = 0
    for seq, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        chunks = get_chunks(seq, chunk_size)
        mates = [search_mates(kmer_dict, chunk, kmer_size) for chunk in chunks]
        parents = mates[0]
        for i, mate in enumerate(mates):
            parents = common(parents, mate)

        print("nb parents", len(parents))
        perc_id_matrix = []
        if len(parents) >= 2:
            perc_id_matrix = calcul_identity_matrix(chunks, parents[:2], chunk_size, non_chimeres)
        print('here is perc id matrix')
        print(perc_id_matrix)
        if not detect_chimera(perc_id_matrix):
            kmer_dict = get_unique_kmer(kmer_dict, seq, id_seq, kmer_size)
            non_chimeres.append(seq)
            id_seq += 1
            yield seq, count

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Apply chimera_removal clustering for amplicon file
    """
    return [rm_chimer for rm_chimer in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)]

def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """
    write OTUs to output_file with specific format
    """
    with open(output_file, "w") as file:
        for i, _ in enumerate(OTU_list):
            file.write(">OTU_" + str(i + 1) + " occurrence:"+ str(OTU_list[i][1]) + "\n")
            file.write(fill(str(OTU_list[i][0])))
            file.write("\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file = args.amplicon_file
    if isfile(file):
        otu_cluster = abundance_greedy_clustering(file, args.minseqlen,
                  args.mincount, args.chunk_size, args.kmer_size)
        write_OTU(otu_cluster, args.output_file)


if __name__ == '__main__':
    main()
