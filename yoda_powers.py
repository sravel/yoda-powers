#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################
# Modules
##################################################
# Python modules
import os
import glob
import re
from pathlib import Path

##################################################
# Variables Globales

__doc__ = """
the toolbox module

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


def existant_file(path):
    """
    'Type' for argparse - checks that file exists and return the absolute path as PosixPath() with pathlib

    Notes:
        function need modules:

        - pathlib
        - argparse


    Arguments:
        path (str): a path to eixstant file

    Returns:
        :class:`PosixPath`: ``Path(path).resolve()``

    Raises:
         ArgumentTypeError: If file `path` does not exist.
         ArgumentTypeError: If `path` is not a valid file.

    Examples:
        >>> import argparse
        >>> parser = argparse.ArgumentParser(prog='test.py', description='''This is demo''')
        >>> parser.add_argument('-f', '--file', metavar="<path/to/file>",type=existant_file, required=True,
            dest='path_file', help='path to file')

    """
    from argparse import ArgumentTypeError
    from pathlib import Path

    if not Path(path).exists():
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise ArgumentTypeError(f'ERROR: file "{path}" does not exist')
    elif not Path(path).is_file():
        raise ArgumentTypeError(f'ERROR: "{path} " is not a valid file')

    return Path(path).resolve()




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


def sort_human(in_list, _nsre=None):
    """
    Sort a :class:`list` with alpha/digit on the way that humans expect,\n
    use list.sort(key=sort_human) or\n
    sorted(list, key=sort_human)).

    Arguments:
        in_list (:obj:`list`): a python :class:`list`
        _nsre (:obj:`re.compil`, optional): re expression use for compare , defaults re.compile('([0-9]+)'

    Returns:
        list: sorted with human sort number

    Example:
        >>> list_to_sorted = ["something1","something32","something17","something2","something29","something24"]
        >>> print(sorted(list_to_sorted, key=sort_human))
        ['something1', 'something17', 'something2', 'something25', 'something29', 'something32']
        >>> list_to_sorted.sort(key=sort_human)
        >>> print(list_to_sorted)
        ['something1', 'something17', 'something2', 'something25', 'something29', 'something32']

    """
    from warnings import warn
    import re
    if not _nsre:
        _nsre = re.compile('([0-9]+)')
    try:
        return [int(text) if text.isdigit() else f"{text}".lower() for text in re.split(_nsre, in_list)]
    except TypeError as e:
        if not isinstance(in_list, int):
            warn(f"Yoda_powers::sort_human : element '{in_list}' on the list not understand so don't sort this element\n",
                          SyntaxWarning, stacklevel=2)
            return in_list

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

#################################################
# CLASS
#################################################

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

