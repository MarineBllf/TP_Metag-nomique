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
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import numpy as np 
np.int = int


__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file,'rt') as fasta_dezip:   
        seq = ""
        for line in fasta_dezip:
            line_clear = line.strip()
            if line_clear.startswith('>'):
                if len(seq)  >= minseqlen : 
                    yield seq
                    seq = ""
                else : 
                    seq = ""
            else:
                seq = seq + line_clear
        if len(seq)  >= minseqlen :
            yield seq
    pass


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    generateur_amplicon = read_fasta(amplicon_file,minseqlen)
    dico_amplicon = {}
    for sequence in generateur_amplicon :
        if sequence in dico_amplicon:
            dico_amplicon[sequence] += 1
        else:
            dico_amplicon[sequence] = 1
    
    for key in sorted(dico_amplicon, key=dico_amplicon.get,reverse=True):
        #print(key)
        if dico_amplicon[key] >= mincount : 
            yield key, dico_amplicon[key]
    pass



def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences
    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq1 = str( alignment_list[0])
    seq2 = str(alignment_list[1])
    longueur = len(seq1)
    identity = 0
    for i in range(0,longueur):
        if seq1[i] == seq2[i] : 
            identity+=1

    pc_identity = (identity/longueur)*100
    return pc_identity
    pass

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int,chunk_size: int=0, kmer_size: int=0) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    OTU_liste = []

    seq_N_occ = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    
    for seq, count in seq_N_occ:
        is_otu = True
        for existing_otu, existing_count in OTU_liste:
            alignement = nw.global_align(seq, existing_otu, gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
            align_id = get_identity(alignement)
            if align_id > 97 and count <= existing_count:
                is_otu = False
                break
        if is_otu:
            OTU_liste.append([seq, count])

    return OTU_liste
    # OTU_liste = []
    
    # seq_N_occ = list(dereplication_fulllength(amplicon_file,minseqlen,mincount))
    # occ = seq_N_occ[0][1]
    # for i in range (0,len(seq_N_occ)) : 
    #     if i == 0 : 
    #         OTU_liste.append(seq_N_occ[0])
    #         print(OTU_liste)
    #         #print(occ)
    #         occ = seq_N_occ[0][1]
    #         #print(len(seq_N_occ))
    #     else:
    #         if seq_N_occ[i][1] == occ :
    #             OTU_liste.append(seq_N_occ[i])
    #             print(occ)
    #             print(OTU_liste)
    #         elif seq_N_occ[i][1] <= occ : 
    #             flag_occ = True
    #             for j in range(len(OTU_liste)):
    #                 alignement = nw.global_align(seq_N_occ[j][0], OTU_liste[j][0], gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
    #                 align_id = get_identity(alignement)
    #                 if align_id > 97 : 
    #                     flag_occ = False
    #             if flag_occ : 
    #                  OTU_liste.append(seq_N_occ[i])
    # return (OTU_liste)
    pass


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as fasta_file:
        for i, (seq, occ) in enumerate(OTU_list, start=1):
            header = f">OTU_{i} occurrence:{occ}\n"
            wrapped_sequence = textwrap.fill(seq, width=80)
            fasta_file.write(header + wrapped_sequence + "\n")

    pass


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    list_OTUs = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount) 
    write_OTU(list_OTUs, args.output_file)


if __name__ == '__main__':
    main()
