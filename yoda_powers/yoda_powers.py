#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################
# Modules
##################################################
# Python modules
import os
import subprocess
import glob
import re
import warnings
from pathlib import Path

# BIO Python modules
from Bio import SeqIO

# ParseGFF modules
from collections import namedtuple
import gzip
import urllib

##################################################
# Variables Globales
__version__ = "1.0.0"

__doc__ = """
:author: Sebastien Ravel
:contact: sebastien.ravel@cirad.fr
:date: 22-07-2020
:version: """ + __version__ + """

Use it to import very handy functions.


Install
-------
::

    # not working yet, module in development
    pip3 install yoda_powers

Install developing version
--------------------------

If you want to use an **unofficial version** of the ``yoda_powers`` module, you need to work from
a clone of this ``git`` repository.
following actions:

1. Clone from github ::

    $ git clone https://github.com/sravel/yoda-powers.git

2. Go in the cloned directory ::

    $ cd yoda-powers

3. Install ``yoda-powers`` in editable mode ::

    $ sudo pip3 install -e .

Required Module install
-----------------------

This module run with Python 3.x and not Python 2.x

** Modules install with :class:`yoda_powers`:**
    - BioPython
    - pyvcf
    - pyfaidx

more info at https://yoda-powers.readthedocs.io/en/latest/

Exemple of usage
----------------

Example:

>>> from yoda_powers import dict2txt
>>> dico = {'key1':'value1','key2':'value2','key3':'value3'}
>>> dict2txt(dico)
key1	value1
key2	value2
key3	value3

Table of contents
-----------------

"""


##################################################
# Functions

def compare_list(list1, list2):
    """
    Function to compare two list and return common, uniq1 and uniq2

    Arguments:
        list1 (list): the first python :class:`list`
        list2 (list): the second python :class:`list`

    Returns:
        list: common, u1, u2
        common: the common elements of the 2 list,
        u1: uniq to list1,
        u2: uniq to list2

    Notes:
        ens1 = set([1, 2, 3, 4, 5, 6])\n
        ens2 = set([2, 3, 4])\n
        ens3 = set([6, 7, 8, 9])\n
        print(ens1 & ens2) set([2, 3, 4]) car ce sont les seuls à être en même temps dans ens1 et ens2\n
        print(ens1 | ens3) set([1, 2, 3, 4, 5, 6, 7, 8, 9]), les deux réunis\n
        print(ens1 & ens3) set([6]), même raison que deux lignes au dessus\n
        print(ens1 ^ ens3) set([1, 2, 3, 4, 5, 7, 8, 9]), l'union moins les éléments communs\n
        print(ens1 - ens2) set([1, 5, 6]), on enlève les éléments de ens2

    Examples:
        >>> l1 = [1, 2, 3, 4, 5, 6]
        >>> l2 = [6, 7, 8, 9]
        >>> com, u1, u2 = compare_list(l1, l2)
        >>> print(com)
        [6]
        >>> print(u1)
        [1, 2, 3, 4, 5]
        >>> print(u2)
        [7, 8, 9]

    """

    ens1 = set(list1)
    ens2 = set(list2)
    common = list(ens1 & ens2)
    uniq1 = list(ens1 - ens2)
    uniq2 = list(ens2 - ens1)
    return sorted(common, key=sort_human), sorted(uniq1, key=sort_human), sorted(uniq2, key=sort_human)


def concat_fasta_files(path_directory):
    """
    Return a fasta dictionnary of concatenation fasta file's find in directory ("fasta", "fa", "fas")

    Warning:
        Sequence on fasta must have the same name

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        path_directory (str): a path to fasta file directory

    Returns:
        :class:`dict`: python dict with the concatenation of fasta filename in path_directory ( file with extention "fa", "fasta", "fas" )

    Raises:
         ValueError: If `path_directory` does not exist.
         ValueError: If `path_directory` is not a valid directory.

    Example:
        >>> dico_concat = concat_fasta_files('path/to/directory/')
        >>> print(dico_concat)
            {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
    """
    from pathlib import Path
    from Bio import SeqIO

    if not Path(path_directory).exists():
        raise ValueError(f'ERROR: directory "{path_directory}" does not exist')
    elif not Path(path_directory).is_dir():
        raise ValueError(f'ERROR: "{path_directory} " is not a valid directory')

    path_directory = Path(path_directory).resolve()
    fa_files = path_directory.glob("*.fa")
    fasta_files = path_directory.glob("*.fasta")
    fas_files = path_directory.glob("*.fas")

    output_dico_seqs = {}

    for fasta_file in (*fa_files, *fasta_files, *fas_files):
            # print(fasta_file)
        with open(fasta_file.as_posix(), "rU") as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        for seq_name in sorted(record_dict.keys()):
             if seq_name not in output_dico_seqs.keys():
                name = record_dict[seq_name].id
                seq = record_dict[seq_name].seq
                output_dico_seqs[name] = seq
             else:
                name = record_dict[seq_name].id
                seq = record_dict[seq_name].seq
                output_dico_seqs[name] += seq

    return output_dico_seqs

def convert_fasta_2_nexus(path_directory, path_directory_out):
    """
    Return the number of fasta file's convert  find in directory ("fasta", "fa", "fas") where are converted

    Warning:
        Sequence on fasta must align and have the same length

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        path_directory (str): a path to fasta file directory
        path_directory_out (str): a directory path to write nexus file

    Returns:
        :class:`int`: the number of file converted

    Raises:
         ValueError: If `path_directory` does not exist.
         ValueError: If `path_directory` is not a valid directory.
         ValueError: If fasta is not align.

    Example:
        >>> nb_file = convert_fasta_2_nexus('path/to/directory/',' path/to/directory/')
        >>> print(nb_file)
            "4172"
    """
    from pathlib import Path
    from Bio import AlignIO
    from Bio.Nexus import Nexus

    if not Path(path_directory).exists():
        raise ValueError(f'ERROR: directory "{path_directory}" does not exist')
    elif not Path(path_directory).is_dir():
        raise ValueError(f'ERROR: "{path_directory} " is not a valid directory')

    path_directory = Path(path_directory).resolve()
    path_directory_out = Path(path_directory_out).resolve()
    if not path_directory_out.exists():
        path_directory_out.mkdir()

    fa_files = path_directory.glob("*.fa")
    fasta_files = path_directory.glob("*.fasta")
    fas_files = path_directory.glob("*.fas")

    # general variables:
    minimal_record = f'#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=dna; end;'
    count_convert = 0

    for fasta_file in (*fa_files, *fasta_files, *fas_files):
        try:
            count_convert += 1
            print(Path(fasta_file))
            basename = Path(fasta_file).stem
            alignment = AlignIO.read(fasta_file.as_posix(), format='fasta')
            n = Nexus.Nexus(minimal_record)
            n.alphabet = alignment._alphabet
            for record in alignment:
                n.add_sequence(record.id.replace("-", "_"), str(record.seq))
            n.write_nexus_data(f"{path_directory_out.as_posix()}/{basename}.nex", interleave=False)
        except ValueError as e:
            raise ValueError(f"ERROR on file {fasta_file}, with message: {e}, please check")

    return count_convert


def dict_2_txt(dico, sep="\t"):
    """
    Function that takes a dictionary and returns a string with separator::

    Arguments:
        dico (:obj:`dict`): the python dict to translate to formated string.
        sep (:obj:`str`): the separator for join . Defaults to '\\t'.

    Returns:
        str: formated dict to string

    Example:
        >>> dico = {"key1":"value1","key2":"value2","key3":"value3"}
        >>> dict_2_txt(dico)
        key1	value1
        key2	value2
        key3	value3
        >>> dict_2_txt(dico, sep=";")
        key1;value1
        key2;value2
        key3;value3

    Warning: if the value of the dict is list or dictionary, the display will be plain and without formatting data.
    """

    return "".join([f"{key}{sep}{dico[key]}\n" for key in sorted(dico.keys(), key=sort_human)])

def dict_2_fasta(dico, fasta_out):
    """
    Function that takes a dictionary where key are ID and value Seq, and write a fasta file.

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        dico (dict): python dict with ID in key and Seq on values
        fasta_out (str): a directory path to write nexus file

    Returns:
        :class:`str`: the output fasta file name

    Example:
        >>> dico = {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
        >>> dict_2_fasta(dico)
        >Seq1
        ATGCTGCAGTAG
        >Seq2
        ATGCCGATCGATG
        >Seq3
        ATGCTCAGTCAGTAG
    """
    from pathlib import Path
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    with open(Path(fasta_out).resolve().as_posix(), "w") as output_handle:
        for seq_name in sorted(dico.keys()):
            seq_obj = Seq(dico[seq_name])
            record = SeqRecord(seq_obj, id=seq_name, name=seq_name, description="")
            SeqIO.write(record, output_handle, "fasta")
    return output_handle.name


def max_key_dict(dico):
    """
    Function return the key of max value in dico values()

    Arguments:
        dico (:obj:`dict`): a python :class:`dict`

    Returns:
        str: key of the dict

    Example:
        >>> dico = {"A":0.5, "C":0.7, "T":0.01, "G":0.9}
        >>> key_max = max_key_dict(dico)
        >>> print(key_max)
        G
    """
    return max(dico, key=dico.get)


def sort_human(list, _nsre=re.compile('([0-9]+)')):
    """
    Sort a :class:`list` with alpha/digit on the way that humans expect,\n
    use list.sort(key=sort_human) or\n
    sorted(list, key=sort_human)).

    :param list: a python :class:`list`
    :type list: :class:`list`
    :param _nsre: re expression use for compare , defaults re.compile('([0-9]+)'
    :type _nsre: re.compile, optional
    :return: :class:`list` sorted with human sort number
    :rtype: :class:`list`

    :Example:

    >>> list_to_sorted = ["something1","something32","something17","something2","something29","something24"]
    >>> print(list_to_sorted.sort(key=sort_human))
    ['something1', 'something17', 'something2', 'something25', 'something29', 'something32']
    >>> print(sorted(list_to_sorted, key=sort_human))
    ['something1', 'something17', 'something2', 'something25', 'something29', 'something32']

    """
    try:
        return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, list)]
    except TypeError:
        if not isinstance(list, int):
            warnings.warn(f"WARNNING Yoda_powers::sort_human : List {list} value not understand so don't sort \n",
                          SyntaxWarning, stacklevel=2)
        return list


def existant_file(path):
    """
    'Type' for argparse - checks that file exists but does not open by default.
    return the absolute path as PosixPath() with pathlib

    :param path: a file path
    :type path: :class:`str`
    :return: Absolute path of file
    :rtype: :class:`PosixPath`
    :raise: :class:`argparse.ArgumentTypeError` if file path does not exist or not a file

    """
    from argparse import ArgumentTypeError
    if not Path(path).exists():
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise ArgumentTypeError(f'ERROR: file "{path}" does not exist')
    elif not Path(path).is_file():
        raise ArgumentTypeError(f'ERROR: "{path} " is not a valid file')

    return Path(path).resolve()




def dict_list2txt(dico, sep="\t"):
    """
    Function that takes a dictionary of list and returns a tabular string with::

        "key\\tvalue\\n".

    :param dico: a python dictionary
    :type dico: dict()
    :param sep: separator character default tab
    :type sep: str()
    :rtype: str()
    :return: string with "key\t\tvalue1\tvalue1\n

    Example:
        >>> dico = {"key1":["value1","value1"], "key2":["value2","value2"],"key3":["value3","value3"]}
        >>> dict_list2txt(dico)
        key1	value1	value1
        key2	value2	value2
        key3	value3	value3

    """

    txt_output = ""
    for key in sorted(dico.keys(), key=sort_human):
        txt_output += f"{key}{sep}{sep.join(sorted(dico[key], key=sort_human))}\n"
    return txt_output


def dict_dict2txt(dico, first="Info", sep="\t"):
    """
    Function that takes a dictionary and returns a tabular string with::

        "key\\tvalue\\n".

    :param dico: a python dictionary
    :type dico: dict()
    :param first: string for name first column
    :type first: str()
    :rtype: str()
    :param: sep: separator character default tab
    :type sep: str()
    :return: string with "key\\tvalue\\n

    Example:
        >>> dico = {"Souche1":{"NUM":"171","MIN":"2042","MAX":"3133578","N50 BP":"938544","N50 NUM":"11"},
                    "Souche2":{"NUM":"182","MIN":"5004","MAX":"74254","N50 BP":"45245"}}
        >>> dict_dict2txt(dico,"souches")
        souches	NUM	MIN	MAX	N50 BP	N50 NUM
        Souche1	171	2042	3133578	938544	11
        Souche2	182	5004	74254	45245	None

    """

    txt_output = f'{first}{sep}{sep.join(list(dico.values())[0].keys())}\n'
    for row_name in sorted(dico.keys(), key=sort_human):
        value = sep.join([str(dico[row_name][key2]) if key2 in dico[row_name].keys() else "None" for key2 in
                          list(dico.values())[0].keys()])
        txt_output += f"{row_name}{sep}{value}\n"
    return txt_output





def fasta2dict(filename):
    """
    Function that take a file name (fasta), and return a dictionnary of sequence

    :param filename: a fasta file
    :type filename: file in fasta format
    :rtype: record_dict()
    :return: dict() - dictionnary with keys are Id and value SeqRecord() fields
    :requires: this function require BIO Python modules:
    :requires: from Bio import SeqIO
    :requires: from Bio.SeqRecord import SeqRecord
    :requires: from Bio.Seq import Seq
    :requires: from Bio.Alphabet import SingleLetterAlphabet

    Example:
        >>> filename = "sequence.fasta"
        >>> fasta2dict(filename)
        {">Seq1":"SeqRecord()"}
    """

    # chargement du fasta des MGG en mémoire
    with open(filename, "rU") as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def lenSeq2dict(filename):
    """
    Function that take a file name (fasta), and return a dictionnary with length of sequence

    :param filename: a fasta file
    :type filename: file in fasta format
    :rtype: record_dict()
    :return: dict() - contain length of sequences
    :requires: this function require fasta2dict(filename)

    Example:
        >>> filename = "sequence.fasta"
        >>> lenSeq2dict(filename)
        {">Seq1":20154}
    """

    dicoLenMGG = {}
    record_dict = fasta2dict(filename)
    for gene in sorted(record_dict.keys(), key=sort_human):
        if record_dict[gene].id not in dicoLenMGG:
            lenseq = len(record_dict[gene].seq)
            dicoLenMGG[gene] = int(lenseq)
    return dicoLenMGG


def nbSeqInFile2dict(pathDirectory):
    """
    Function  that take a Path Directory and returna dictionnary with number of sequences in fasta file's

    :param pathDirectory: a directory Path
    :type pathDirectory: Path
    :rtype: dict1(), dict2()
    :return: - contient le nombre de sequences dans les fichiers (key = nom de fichier value =  nombre de sequences)\n
             - contient le nombre de fichier qui ont x sequences (key = nombre de sequence =  nombre de fichier)
    :raise print: print("ERROR: Sequence: "+nameFichier+" allready read") with nameFichier is the current file read.

    Example:
        >>> dico1,dico2 = nbSeqInFile2dict(path/to/directory/)
        >>> print(dict2txt(dico1))
        ./out/gemo10_4497_ortho_rename_add.fasta	58
        ./out/gemo10_6825_ortho_rename_add.fasta	59
        ./out/gemo10_3497_ortho_rename_add.fasta	59
        ./out/gemo10_6254_ortho_rename_add.fasta	59
        >>> print(dict2txt(dico2))
        58	1
        59	3
    """

    dicoNbSeqInFiles = {}
    dicoNbFilesNbSouche = {}
    # récupération des fichiers du repertoire
    if pathDirectory[-1] != "*":
        pathDirectory += "*"
    directoryFiles = glob.glob(pathDirectory)
    # Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
    for fichier in directoryFiles:
        try:
            nameFichier = fichier.split("/")[-1].split(".")[0]
            extentionFichier = fichier.split("/")[-1].split(".")[1]
        except:
            extentionFichier = "directory"
        if extentionFichier in ["fasta", "fa", "fas"]:
            # Ouverture des sequences fasta et chargement dans dictionnaire
            handle = open(fichier, "r")
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            handle.close()
            nbseq = len(record_dict.keys())
            # print(nameFichier+"\t"+str(nbseq))
            if nameFichier not in dicoNbSeqInFiles.keys():
                dicoNbSeqInFiles[nameFichier] = nbseq
            else:
                print("ERROR: Sequence: " + nameFichier + " allready read")

            dicoNbFilesNbSouche[nbseq] = (dicoNbFilesNbSouche.get(nbseq, 1)) + 1

    # if nbseq not in dicoNbFilesNbSouche.keys():
    # dicoNbFilesNbSouche[nbseq] = 1
    # else:
    # dicoNbFilesNbSouche[nbseq] += 1

    return dicoNbSeqInFiles, dicoNbFilesNbSouche


def readable_dir(prospective_dir):
    """
    Check if directory exist and if is readable

    :param prospective_dir: a directory Path
    :type prospective_dir: Path

    :raise error: raise argparse.ArgumentTypeError(" :{0} is not a valid path".format(prospective_dir))
    :raise error: raise argparse.ArgumentTypeError(" :{0} is not a readable dir".format(prospective_dir))

    Example:
        >>> parser = argparse.ArgumentParser(prog='make_structure_dir.py', description='''This Programme make arborescence of rep of programme structure''')
        >>> paths = parser.add_argument_group('Input PATH for running')
        >>> paths.add_argument('-p', '--path', metavar="<path/to/>", type = readable_dir, required=True, dest = 'pathParam', help = 'Path to ')
    """

    if not os.path.isdir(prospective_dir):
        raise argparse.ArgumentTypeError(" :{0} is not a valid path".format(prospective_dir))
    if os.access(prospective_dir, os.R_OK) == False:
        raise argparse.ArgumentTypeError(" :{0} is not a readable dir".format(prospective_dir))


def replace_all(repls, str):
    """
    Function that take a dictionnary and text variable and return text variable with replace 'Key' from dictionnary with 'Value'.

    :param repls: a python dictionary
    :type repls: dict()
    :param str: a string where remplace some words
    :type str: str()
    :rtype: str()
    :return: - txt with replace 'Key' of dictionnary with 'Value' in the input txt

    Example:
        >>> text =  "i like apples, but pears scare me"
        >>> print(replace_all({"apple": "pear", "pear": "apple"}, text))
        i like pears, but apples scare me
    """

    return re.sub('|'.join(re.escape(key) for key in repls.keys()), lambda k: repls[k.group(0)], str)


def loadInList(filename):
    """
    Load file in list() and then remove \\n at end of line

    :param filename: a file
    :type filename: file
    :rtype: list()
    :return: - list of row's file without \\n
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> rows = loadInList(filename)
        >>> rows
        ["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
    """
    with open(filename, "r") as fileIn:
        return [line.rstrip() for line in fileIn.readlines()]


def loadInListCol(filename, col):
    """
    Load a column of file in list() and remove \\n at end of line

    :param filename: a file
    :type filename: file
    :param col: a int of keep column
    :type col: int
    :rtype: list()
    :return: - list of row's file from column without \\n if end column
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> rows = loadInListCol(filename, 0)
        >>> rows
        ["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
    """

    with open(filename, "r") as fileIn:
        return [line.rstrip().split("\t")[col] for line in fileIn.readlines()]


def loadInListWithHeader(filename):
    """
    Load file in two list(): return header list and rows list

    :param filename: a file
    :type filename: file
    :rtype: list(), list()
    :return: - header liste\n
             - list of list of row's file
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> header, rows = loadInListWithHeader(filename)
        >>> header
        "['head1','head2','head3']
        >>> rows
        [["line1col1","line1col2","line1col3"],["line2col1","line2col2","line2col3"]]
    """

    with open(filename, "r") as fileIn:
        list = fileIn.readlines()
    header = list[0].rstrip().split("\t")
    listgood = [line.rstrip() for line in list[1:]]
    return header, listgood


def loadInDict(filename):
    """
    Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valueare other column.

    :param filename: a file
    :type filename: file
    :rtype: dict()
    :return: - dict of row's file without \\n with key is first column and value list of other column
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> dico = loadInDict(filename)
        >>> dico
        {
        "col1",["col2","col3"],
        "indiv1",["valeurcol2","valeurcol3"],
        "indiv2",["valeurcol2","valeurcol3"]
        }
    """

    dicoOut = {}
    with open(filename) as filein:
        for line in filein:
            tabLine = line.rstrip().split("\t")
            # print(tabLine[0], tabLine[1])
            if tabLine[0] not in dicoOut.keys():
                dicoOut[tabLine[0]] = [] + tabLine[1:]
    # else:
    # dicoOut[tabLine[0]].append(tabLine[1])
    return dicoOut


def loadInDictList(filename):
    """
    Load file in Dict() and remove \\n at end of line, then add first column in key of dict and value are list of column 2.

    :param filename: a file
    :type filename: file
    :rtype: dict()
    :return: - dict of row's file without \\n with key is first column and value list of other column
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> dico = loadInDictList(filename)
        >>> dico
        {
        "scaffold1",["1000","2000"],
        "scaffold12",["2000","5000"]
        }
    """

    dicoOut = {}
    with open(filename) as filein:
        for line in filein:
            tabLine = line.rstrip().split("\t")

            if tabLine[0] not in dicoOut.keys():
                dicoOut[tabLine[0]] = [] + tabLine[1:]
            else:
                dicoOut[tabLine[0]].append(tabLine[1])
    return dicoOut


def loadInDictCol(filename, columnkey, columnvalue):
    """
    Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valu specify column.

    :param filename: a file
    :type filename: file
    :param columnkey: int of column
    :type columnkey: int
    :param columnvalue: int of column
    :type columnvalue: int
    :rtype: dict()
    :return: - dict of row's file without \\n with key is first column and value column number pass
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> dico = loadInDict(filename,columnkey=1,columnvalue=3 )
        >>> dico
        {
        "col1","col3",
        "indiv1","valeurcol3",
        "indiv2","valeurcol3"
        }
    """

    dicoOut = {}
    with open(filename) as filein:
        for line in filein:
            tabLine = line.rstrip().split("\t")
            # print(tabLine[0], tabLine[1])
            if tabLine[columnkey] not in dicoOut.keys():
                dicoOut[tabLine[columnkey]] = tabLine[columnvalue]
    return dicoOut


def loadInDictLine(filename):
    """
    Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valueare other column.

    :param filename: a file
    :type filename: file
    :rtype: dict()
    :return: - dict of row's file without \\n with key is first column and value list of other column
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> dico = loadInDictLine(filename)
        >>> dico
        {
        "col1",[line1],
        "indiv1",[line2],
        "indiv2",[line3]
        }
    """

    dicoOut = {}
    with open(filename) as filein:
        for line in filein:
            tabLine = line.rstrip().split("\t")
            # print(tabLine[0], tabLine[1])
            if tabLine[0] not in dicoOut.keys():
                dicoOut[tabLine[0]] = line
    return dicoOut


def loadInDictDict(filename):
    """
    Load a file with header in dictDict().

    :param filename: a file
    :type filename: file
    :rtype: dict()
    :return: - dict of dict
    :warn: Use this function with small file !!! except more RAM are use and crash systeme.

    Example:
        >>> dico = loadInDictDict(filename)
        >>> dico
        {
        "indiv1",{"headerCol2":"toto","headerCol3":"tata"},
        "indiv2",{"headerCol2":"tutu","headerCol3":"titi"},
        "indiv3",{"headerCol2":"tete","headerCol3":"tata"},
        }
    """

    dicoOut = {}
    with open(filename) as filein:
        header = filein.readline().rstrip().split("\t")
        for line in filein:
            tabLine = line.rstrip().split("\t")
            if tabLine[0] not in dicoOut.keys():
                dicoOut[tabLine[0]] = {}
                i = 1
                for head in header[1:]:
                    dicoOut[tabLine[0]][head] = tabLine[i]
                    i += 1
            else:
                print("ERROR key %s already load exit" % tabLine[0])
    return dicoOut


def extractListFromFasta(sequenceFile, FileList):
    """
    Function who use 2 files:\\n
            - fasta : file with all Sequences
            - list : file with name of sequence to extract from fasta

    then return dict() with only Sequence in listed file.

    :param sequenceFile: a file with all Sequences
    :type sequenceFile: file
    :param FileList: a file with name of sequence to extract from fasta
    :type FileList: file
    :rtype: dict() and int()
    :return: - dict with only Sequence in listed file and nbtotal seq in file
    :requires: this function require fasta2dict(filename) and loadInList(filename)

    Example:
        >>> dictSequences = extractListFromFasta(sequenceFile, FileList)
        >>> dictSequences
        {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
    """

    dicoOutput = {}
    # Ouverture des sequences fasta MGG et chargement dans dictionnaire
    dictSequences = fasta2dict(sequenceFile)
    # ouverture des identifiants a garder
    listKeep = loadInList(FileList)
    keep = 0
    noKeep = 0
    noKeepID = []
    # for ID, record in dictSequences.items():
    for ID in sorted(dictSequences.keys(), key=sort_human):
        record = dictSequences[ID]
        if ID in listKeep:
            keep += 1
            dicoOutput[ID] = record
        else:
            noKeepID.append(ID)
            noKeep += 1
    # print("seq keep:"+str(keep))
    # print("seq nokeep:"+str(noKeep))
    # print("ID nokeep:"+str(noKeepID))
    total = noKeep + keep
    return dicoOutput, total


def extractInverseListFromFasta(sequenceFile, FileList):
    """
    Function who use 2 files:\\n
            - fasta : file with all Sequences
            - list : file with name of sequence to extract from fasta

    then return dict() with only Sequence not in listed file.

    :param sequenceFile: a file with all Sequences
    :type sequenceFile: file
    :param FileList: a file with name of sequence to not extract from fasta
    :type FileList: file
    :rtype: dict() and int()
    :return: - dict with only Sequence in listed file and nbtotal seq in file
    :requires: this function require fasta2dict(filename) and loadInList(filename)

    Example:
        >>> dictSequences = extractListFromFasta(sequenceFile, FileList)
        >>> dictSequences
        {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
    """

    dicoOutput = {}
    # Ouverture des sequences fasta MGG et chargement dans dictionnaire
    dictSequences = fasta2dict(sequenceFile)
    # ouverture des identifiants a garder
    listnotKeep = loadInList(FileList)
    keep = 0
    noKeep = 0
    noKeepID = []
    # for ID, record in dictSequences.items():
    for ID in sorted(dictSequences.keys(), key=sort_human):
        record = dictSequences[ID]
        if ID not in listnotKeep:
            keep += 1
            dicoOutput[ID] = record
        else:
            noKeepID.append(ID)
            noKeep += 1
    total = noKeep + keep
    # print("seq keep:"+str(keep))
    # print("seq nokeep:"+str(noKeep))
    # print("ID nokeep:"+str(noKeepID))
    return dicoOutput, total


def lsDirToList(pathDirectory):
    """
    Return a list of file and directory find in directory

    :param pathDirectory: a directory Path
    :type pathDirectory: Path
    :rtype: list()
    :return: list of filename in pathDirectory

    Example:
        >>> lsDirectory = lsDirToList(path/to/directory/)
        >>> print(lsDirectory)
        ["./out/gemo10_4497_ortho_rename_add.fasta", "./out/gemo10_6825_ortho_rename_add.fasta", "./out/gemo10_3497_ortho_rename_add.fasta", "./out/rename/"]
    """

    if pathDirectory[-1] != "/":
        pathDirectory += "/"
    if pathDirectory[-1] != "*":
        pathDirectory += "*"
    lsFiles = glob.glob(pathDirectory)
    return lsFiles


def lsFastaInDirToList(pathDirectory):
    """
    Return a list of fasta file's find in directory ("fasta", "fa", "fas")

    :param pathDirectory: a directory Path
    :type pathDirectory: Path
    :rtype: list()
    :return: list of fasta filename in pathDirectory ( file with extention "fa", "fasta", "fas" )

    Example:
        >>> lsDirectory = lsFastaInDirToList(path/to/directory/)
        >>> print(lsDirectory)
        ["./out/gemo10_4497_ortho_rename_add.fasta", "./out/gemo10_6825_ortho_rename_add.fasta", "./out/gemo10_3497_ortho_rename_add.fasta"]
    """

    lsFilesFasta = []
    if pathDirectory[-1] != "/":
        pathDirectory += "/"
    if pathDirectory[-1] != "*":
        pathDirectory += "*"
    directoryFiles = glob.glob(pathDirectory)
    # Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
    for fichier in directoryFiles:
        try:
            nameFichier = fichier.split("/")[-1].split(".")[0]
            extentionFichier = fichier.split("/")[-1].split(".")[-1]
        except:
            extentionFichier = "directory"
        if extentionFichier in ["fasta", "fa", "fas"]:
            lsFilesFasta.append(fichier)

    return sorted(lsFilesFasta)


def lsExtInDirToList(pathDirectory, extentionFichierKeep):
    """
    Return a list of 'ext' file's find in directory (exemple ext = "txt" or ["txt","py"])

    :param pathDirectory: a directory Path
    :type pathDirectory: Path
    :param extentionFichierKeep: a list or string with extention
    :type extentionFichierKeep: list or string
    :rtype: list()
    :return: list of 'ext' filename in pathDirectory ( file with extention find in param extentionFichierKeep )

    Example:
        >>> lsDirectory = lsExtInDirToList(path/to/directory/,"txt")
        >>> print(lsDirectory)
        ["./out/gemo10_4497_ortho_rename_add.txt", "./out/gemo10_6825_ortho_rename_add.txt", "./out/gemo10_3497_ortho_rename_add.txt"]
    """

    lsFilesFasta = []
    if pathDirectory[-1] != "/":
        pathDirectory += "/"
    if pathDirectory[-1] != "*":
        pathDirectory += "*"
    directoryFiles = glob.glob(pathDirectory)

    # Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
    for fichier in directoryFiles:
        try:
            if "." in fichier.split("/")[-1]:
                nameFichier = fichier.split("/")[-1].split(".")[0]
                extentionFichier = fichier.split("/")[-1].split(".")[-1]
            else:
                nameFichier = fichier.split("/")[-1]
                extentionFichier = ""
        except:
            extentionFichier = "directory"
        if extentionFichier == extentionFichierKeep:
            lsFilesFasta.append(fichier)

    return sorted(lsFilesFasta)

def printcolor(txt, color, noprint=1):
    """	Return the printed color txt format

    :param txt: a string
    :type txt: string
    :param color: a color value
    :type color: string
    :type noprint: int 0=noprint 1=print (default)
    :rtype: string()
    :return: string with acci color for printed
    :warn: List of avail color: reset, hicolor, underline, inverse, fblack, fred, fgreen, fyellow, fblue, fmagenta, fcyan, fwhite, bblack, bred, bgreen, byellow, bblue, bmagenta, bcyan, bwhite

    Example:
        >>> printcolor("il fait beau aujourd'hui","bgreen")
            "\\033[36mil fait beau aujourd'hui"
        >>> txtcolor = printcolor("il fait beau aujourd'hui","bgreen", 0)

    """
    dicoColor = {
            "reset" : "\033[0m", "hicolor": "\033[1m", "underline": "\033[4m", "inverse": "\033[7m",
            "fblack": "\033[30m", "fred": "\033[31m", "fgreen": "\033[32m", "fyellow": "\033[1;33m",
            "fblue" : "\033[34m", "fmagenta": "\033[35m", "fcyan": "\033[36m", "fwhite": "\033[37m",
            "bblack": "\033[40m", "bred": "\033[41m", "bgreen": "\033[42m", "byellow": "\033[43m",
            "bblue" : "\033[44m", "bmagenta": "\033[45m", "bcyan": "\033[46m", "bwhite": "\033[47m",
    }
    if color in dicoColor.keys():
        txtout = dicoColor[color] + txt
        if noprint == 0:
            return txtout
        else:
            print(txtout)
    else:
        txtout = "Error, color value non exist, please check other color\n\n" + txt


def relativeToAbsolutePath(relative):
    """	Return the absolutPath

    :param relative: a string path
    :type relative: string
    :rtype: string()
    :return: absolutePath
    :warn: need subprocess::check_output

    Example:
    >>>print(relative)
    ../test
    >>> pathDirectory = relativeToAbsolutePath(relative)
    >>>print(pathDirectory)
    /home/sebastien/test

    """
    from subprocess import check_output
    if relative[0] != "/":  # The relative path is a relative path, ie do not starts with /
        command = "readlink -m " + relative
        absolutePath = subprocess.check_output(command, shell=True).decode("utf-8").rstrip()
        return absolutePath
    else:  # Relative is in fact an absolute path, send a warning
        absolutePath = relative;
        return absolutePath


#################################################
# CLASS
#################################################

class parseGFF():
    """
    Parser of GFF3 file write in python.
    return an object iterable containt GFFRecord()

    line in GFF3 return:

    Example:
        >>> scaffold_44     prediction      gene    46      6942    0       +       .       ID=gene_1;Name=jgi.p|Mycfi2|180833;portal_id=Mycfi2;proteinId=180833;transcriptId=180833
        >>> GFFRecord(seqid='scaffold_44', source='prediction', type='gene', start=46, end=6942, score=0.0, strand='+', phase=None,
        >>> 		attributes={'portal_id': 'Mycfi2', 'transcriptId': '180833', 'proteinId': '180833', 'Name': 'jgi.p|Mycfi2|180833', 'ID': 'gene_1'}, seq=None, len=6896)

    GFFRecord has attributes can acces with record.value (ex: record.seqid):

    ==========  ===========================================================
    attribut    infos
    ==========  ===========================================================
    seqid       first column of gff3
    source      second column of gff3
    type        third column of gff3 contain type
    start       start position
    end         end position
    score       score
    strand      DNA brin
    phase       phase
    attributes  dict() with key corresponding to GFFAttributes
    seq         if fasta load can add sequence but by default = None
    len         size of sequence
    ==========  ===========================================================

    Example:
        >>> objGFF = parseGFF(gffFile)
        >>> for record in objGFF.parseGFF3():
        >>> 	print(record.seqid)
        >>> 	if record.type == "mRNA" :
        >>> 		transcriptID = record.attributes["transcriptId"]

    """

    def __init__(self, filename):
        # Initialized GeneInfo named tuple. Note: namedtuple is immutable
        self.filename = filename
        self.gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes",
                              "seq", "len"]
        self.GFFRecord = namedtuple("GFFRecord", self.gffInfoFields)

    def parseGFFAttributes(self, attributeString):
        """Parse the GFF3 attribute column and return a dict"""
        if attributeString == ".": return {}
        ret = {}
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
        return ret

    def parseGFF3(self):
        """
        A minimalistic GFF3 format parser.
        Yields objects that contain info about a single GFF3 feature.

        Supports transparent gzip decompression.
        """
        # Parse with transparent decompression
        openFunc = gzip.open if self.filename.endswith(".gz") else open
        with openFunc(self.filename) as infile:
            for line in infile:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                # If this fails, the file format is not standard-compatible
                assert len(parts) == len(self.gffInfoFields) - 2
                # Normalize data
                normalizedInfo = {
                        "seqid"     : None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                        "source"    : None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                        "type"      : None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                        "start"     : None if parts[3] == "." else int(parts[3]),
                        "end"       : None if parts[4] == "." else int(parts[4]),
                        "len"       : None if parts[4] == "." and parts[3] == "." else int(parts[4]) - int(parts[3]),
                        "score"     : None if parts[5] == "." else float(parts[5]),
                        "strand"    : None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                        "phase"     : None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                        "seq"       : None,
                        "attributes": self.parseGFFAttributes(parts[8])
                }
                # Alternatively, you can emit the dictionary here, if you need mutability:
                #	yield normalizedInfo
                yield self.GFFRecord(**normalizedInfo)


class printCol():
    """
    Classe qui ajoute des méthodes à print pour afficher de la couleur

    Example:

    >>> printCol.red("j'affiche en rouge")
    j'affiche en rouge

    """

    __RED = '\033[91m'
    __GREEN = '\033[92m'
    __YELLOW = '\033[93m'
    __LIGHT_PURPLE = '\033[94m'
    __PURPLE = '\033[95m'
    __END = '\033[0m'

    @classmethod
    def red(cls, s):
        print(cls.__RED + str(s) + cls.__END)

    @classmethod
    def green(cls, s):
        print(cls.__GREEN + str(s) + cls.__END)

    @classmethod
    def yellow(cls, s):
        print(cls.__YELLOW + str(s) + cls.__END)

    @classmethod
    def lightPurple(cls, s):
        print(cls.__LIGHT_PURPLE + str(s) + cls.__END)

    @classmethod
    def purple(cls, s):
        print(cls.__PURPLE + str(s) + cls.__END)


class AutoVivification(dict):
    """
    Implementation of perl's autovivification feature.

    Example:

    >>> a = AutoVivification()
    >>> a[1][2][3] = 4
    >>> a[1][3][3] = 5
    >>> a[1][2]['test'] = 6
    >>> print a
    >>> {1: {2: {'test': 6, 3: 4}, 3: {3: 5}}}

    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


# *********************************************** Classe directory *******************
class directory(str):
    """
    Class which derives from string.
    Checks that the string is and path to valid directory and not empty

    Returns object which able to return basename and list file (with an extention, else consider as folder)

    Example:

    >>> inDirectory=directory("/home/sravel/Documents")
    >>> print(inDirectory.pathDirectory())
    >>> /home/sravel/Documents/

    >>> print(inDirectory.listFiles())
    >>> ["File1.txt","File2.pl","file.toto"]
    """

    def __init__(self, pathDirectory=None):
        """
            Initialise variable
        """
        self.listPath = []  # all in the path
        self.listDir = []  # only directory in path
        self.listFiles = []  # only files in path

        self.current_dir = os.path.dirname(os.path.abspath(__file__))

        # Change relative path to absolute path
        self.pathDirectory = relativeToAbsolutePath(pathDirectory)

        # appel les fonctions
        self.testDirExist()
        self.lsInDir()
        self.splitFilesDir()

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        """Fonction qui permet de formater le text de sortie lors du print du dictionnaire"""
        return """
pathDirectory=%s\n
listPath=%s\n
listDir=%s\n
listFiles=%s\n
""" % (self.pathDirectory, str(self.listPath), str(self.listDir), str(self.listFiles))

    def testDirExist(self):
        """Test l'existance du répertoire"""
        if os.path.isdir(self.pathDirectory) != True:
            raise ValueError("ERROR Yoda_powers::Class-directory : path '%s' is not valid path" % self.pathDirectory)

    def lsInDir(self):
        """List all in directory"""
        self.testDirExist()
        if self.pathDirectory[-1] != "/":
            self.pathDirectory += "/"
        if self.pathDirectory[-1] != "*":
            pathDirectoryList = self.pathDirectory + "*"
        self.listPath = glob.glob(pathDirectoryList)

    def lsExtInDirToList(self, ext):
        """List specific extention file in directory"""
        lsFilesFasta = []
        for fichier in self.listPath:
            try:
                if "." in fichier.split("/")[-1]:
                    nameFichier = fichier.split("/")[-1].split(".")[0]
                    extentionFichier = fichier.split("/")[-1].split(".")[-1]
                else:
                    nameFichier = fichier.split("/")[-1]
                    extentionFichier = ""
            except:
                extentionFichier = "directory"
            if extentionFichier in ext:
                lsFilesFasta.append(fichier)
        return sorted(lsFilesFasta, key=sort_human)

    def splitFilesDir(self):
        """list files and list directory"""
        self.lsInDir()
        self.listDir = []  # only directory in path
        self.listFiles = []  # only files in path
        # Ouverture des fichiers du repertoire
        for fichier in self.listPath:
            try:
                if os.path.isdir(fichier) == True:  # C'est un dossier
                    self.listDir.append(str(fichier))
                elif os.path.exists(fichier) == True:  # C'est un fichier
                    self.listFiles.append(str(fichier))
            except:
                raise ValueError(
                        "ERROR Yoda_powers::Class-directory : path '%s' is not valide path contain other type (not files or directory)" % self.pathDirectory)


# *********************************************** Classe StrList *******************
class StrList(list):
    """
    Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien dans la liste "expectedValue"

    Example:

    >>> x=StrList(expectedValue=['true'])
    >>> x.append("true")
    >>> x
    ['true']
    >>> x.append("false")
    Traceback (most recent call last):
        ...
    ValueError: Input value is not in the following list: ['true']
    """

    def __init__(self, expectedValue=None, defaultvalue=None):
        """
        StrList(expectedValue=[])

        :rtype: list()
        :return: list check if not value of past list
        """
        self.expectedValue = expectedValue
        self.defaultValue = defaultvalue
        super(StrList, self).__init__()

    def append(self, value):
        """Redefine append() function to do testItem before append in list"""
        self.testItem(value)
        super(StrList, self).append(value)

    def insert(self, at, value):
        """Redefine insert() function to do testItem before insert in list"""
        self.testItem(value)
        return list.insert(self, at, value)

    def testItem(self, value):
        """Méthode qui test la valeur a ajouter est dans la list des expected"""
        if (value != None) and (value != "") and (self.expectedValue != None):
            if value in self.expectedValue:
                pass
            else:
                raise ValueError("Input value is not in the following list: %s" % self.expectedValue)


# *********************************************** Classe NumberList *******************
class NumberList(list):
    """Classe qui dérive de list et qui vérifie que le nombre ajouté dans l'objet (list) est bien du type défini dans "dataType"
    Returns list

    Example:

    >>> x=NumberList(datatype = int, minValue = 2, maxValue = 100)
    >>> x.append("32")
    >>> x
    ['32']
    >>> x.append("0.5")
    Traceback (most recent call last):
        ...
    ValueError: "0.5" is not type "int"
    >>> x.append("1")
    Traceback (most recent call last):
        ...
    ValueError: "1" is less than "2" ( 2 < value < 100 )
    >>> x.append("1000")
    Traceback (most recent call last):
        ...
    ValueError: "1000" is greater than "100" ( 2 < value < 100 )
    """

    def __init__(self, datatype=None, minValue=None, maxValue=None, defaultvalue=None):
        """
        :rtype: list()
        :return: list check if not value of past list
        """
        self.dataType = datatype
        self.min = minValue
        self.max = maxValue
        self.defaultValue = defaultvalue
        super(NumberList, self).__init__()

    def append(self, value):
        """L.append(object) -- append object to end after do testItem"""
        self.testItem(value)
        super(NumberList, self).append(value)

    def insert(self, at, value):
        self.testItem(value)
        return list.insert(self, at, value)

    def testItem(self, value):
        """Méthode qui test la valeur a ajouter"""

        if (value != None) and (value != ""):
            try:
                self.dataType(value)
            except:
                raise ValueError("\"%s\" is not type \"%s\"" % (value, self.dataType.__name__))
            try:
                if self.min < self.dataType(value):
                    pass
                else:
                    raise ValueError(
                            "\"%s\" is less than \"%s\" ( %s < value < %s )" % (value, self.min, self.min, self.max))
                if self.dataType(value) < self.max:
                    pass
                else:
                    raise ValueError(
                            "\"%s\" is greater than \"%s\" ( %s < value < %s )" % (value, self.max, self.min, self.max))
            except ValueError as infos:
                raise ValueError("%s" % infos)


# ************************************************************* Classe StrSpaceList *******************
class StrSpaceList(list):
    """"Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien du type défini dans "dataType"
    pour une chaine de plusieurs valeurs. test également si la valeur et comprise entre 2 bornes; et le nobre de valeur attendu.
    Returns list

    Example:

    >>> x=StrSpaceList(datatype = float, minValue = 0, maxValue = 1, nbInListMin = 4, nbInListMax = 4, sumOfValues = 1)
    >>> x.append("0.1 0.2 0.3 0.4")
    >>> x
    ['0.1 0.2 0.3 0.4']
    >>> x.append("0.5")
    Traceback (most recent call last):
        ...
    ValueError: You must enter values more than 1 ​​separated by a space ( 4 < nbValue < 4 )
    >>> x.append("2 0.2 0.3 0.4")
    Traceback (most recent call last):
        ...
    ValueError: "2" is greater than "1" ( 0 < value < 1 )
    >>> x.append("-1 0.2 0.3 0.4")
    Traceback (most recent call last):
        ...
    ValueError: "-1" is less than "0" ( 0 < value < 1 )
    """

    def __init__(self, datatype=None, minValue=None, maxValue=None, nbInListMin=None, nbInListMax=None,
                 sumOfValues=None, defaultvalue=None):
        """
        :rtype: list()
        :return: list check if not value of past list
        """
        self.dataType = datatype
        self.min = minValue
        self.max = maxValue
        self.nbInListMin = nbInListMin
        self.nbInListMax = nbInListMax
        self.sumOfValues = sumOfValues
        self.defaultValue = defaultvalue

    def append(self, value):
        """L.append(object) -- append object to end after do testItem"""
        self.testItem(value)
        super(StrSpaceList, self).append(value)

    def insert(self, at, value):
        self.testItem(value)
        return list.insert(self, at, value)

    def testItem(self, value):
        """Méthode qui test la valeur a ajouter"""
        if (value != None) and (value != ""):
            l = value.split(' ')

            for num in l:
                try:
                    self.dataType(num)
                except:
                    raise ValueError("\"%s\" is not type \"%s\"" % (num, self.dataType.__name__))
                try:
                    if self.min < self.dataType(num):
                        pass
                    else:
                        raise ValueError(
                                "\"%s\" is less than \"%s\" ( %s < value < %s )" % (num, self.min, self.min, self.max))
                    if self.dataType(num) < self.max:
                        pass
                    else:
                        raise ValueError(
                                "\"%s\" is greater than \"%s\" ( %s < value < %s )" % (
                                        num, self.max, self.min, self.max))
                except ValueError as infos:
                    raise ValueError("%s" % infos)
            try:
                if len(l) < self.nbInListMin:
                    raise ValueError(
                            "You must enter values more than %s ​​separated by a space ( %s < nbValue < %s )" % (
                                    len(l), self.nbInListMin, self.nbInListMax))
                if len(l) > self.nbInListMax:
                    raise ValueError(
                            "You must enter value less than %s ​​separated by a space ( %s < nbValue < %s )" % (
                                    len(l), self.nbInListMinn, self.nbInListMax))
            except ValueError as infos:
                raise ValueError("%s" % infos)

            if self.sumOfValues != None:
                s = 0
                for num in l:
                    s += self.dataType(num)
                if s != self.sumOfValues:
                    raise ValueError("The values must sum up to %s." % self.sumOfValues)


# ************************************************************* Classe StrSpace *******************
class LetterList(list):
    """Class which derives from list. Checks that the string added to the object (list) contains only the letters from the list "Letters"
    Returns list

    Example:

    >>> x=LetterList(letters=['A', 'C', 'T', 'G'])
    >>> x.append("cagtcgatgcatgctagctagtcagtcat")
    >>> x
    ['cagtcgatgcatgctagctagtcagtcat']
    >>> x.append("CGATCGATCGATCGT")
    >>> x
    ['cagtcgatgcatgctagctagtcagtcat', 'CGATCGATCGATCGT']
    >>> x.append("CGAGFBDCGATCGATCGT")
    Traceback (most recent call last):
        ...
    ValueError: "FBD" is not in "['A', 'C', 'T', 'G']"
    """

    def __init__(self, letters=None, defaultvalue=None):
        """
        :rtype: list()
        :return: list check if not value of past list
        """
        self.Letter = letters
        self.defaultValue = defaultvalue

    def append(self, value):
        """L.append(object) -- append object to end after do testItem"""
        self.testItem(value)
        super(LetterList, self).append(value)

    def insert(self, at, value):
        self.testItem(value)
        return list.insert(self, at, value)

    def testItem(self, value):
        """Méthode qui test la valeur a ajouter"""
        if (value != None) and (value != ""):
            value = value.upper()
            ln = ""
            for l in value:
                if not l in self.Letter:
                    ln += l
            if ln != "":
                raise ValueError("\"%s\" is not in \"%s\"" % (ln, self.Letter))
