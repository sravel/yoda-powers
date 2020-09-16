#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################
# Modules
##################################################
# Python modules
from pathlib import PosixPath, Path
from datetime import datetime
import sys


##################################################
# Functions

def welcome_args(version_arg, parser_arg):
    """
    use this Decorator to add information to scripts with arguments

    Args:
        version_arg: the program version
        parser_arg: the function which return :class:`argparse.ArgumentParser`

    Returns:
        None:

    Notes:
        use at main() decorator for script with :class:`argparse.ArgumentParser`

    Examples:
        >>> from yoda_powers.toolbox import welcome_args
        >>> @welcome_args(version, build_parser())
        >>> def main():
        >>>     # some code
        >>> main()
        >>> ################################################################################
        >>> #                             prog_name and version                            #
        >>> ################################################################################
        >>> Start time: 16-09-2020 at 14:39:02
        >>> Commande line run: ./filter_mummer.py -l mummer/GUY0011.pp1.fasta.PH0014.pp1.fasta.mum
        >>>
        >>> - Intput Info:
        >>>         - debug: False
        >>>         - plot: False
        >>>         - scaff_min: 1000000
        >>>         - fragments_min: 5000
        >>>         - csv_file: blabla
        >>> PROGRAMME CODE HERE
        >>> Stop time: 16-09-2020 at 14:39:02       Run time: 0:00:00.139732
        >>> ################################################################################
        >>> #                               End of execution                               #
        >>> ################################################################################

    """
    def welcome(func):
        def wrapper():
            start_time = datetime.now()
            parser = parser_arg
            version = version_arg
            parse_args = parser.parse_args()
            # Welcome message
            print(
                    f"""{"#" * 80}\n#{Path(parser.prog).stem + " " + version:^78}#\n{"#" * 80}\nStart time: {start_time:%d-%m-%Y at %H:%M:%S}\nCommande line run: {" ".join(sys.argv)}\n""")
            # resume to user
            print(" - Intput Info:")
            for k, v in vars(parse_args).items():
                print(f"\t - {k}: {v}")
            print("\n")
            func()
            print(
                    f"""\nStop time: {datetime.now():%d-%m-%Y at %H:%M:%S}\tRun time: {datetime.now() - start_time}\n{"#" * 80}\n#{'End of execution':^78}#\n{"#" * 80}""")
        return wrapper
    return welcome


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
        path (str): a path to existent file

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
        raise ArgumentTypeError(f'ERROR: "{path}" is not a valid file')

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
        ['something1', 'something2', 'something17', 'something24', 'something29', 'something32']
        >>> list_to_sorted.sort(key=sort_human)
        >>> print(list_to_sorted)
        ['something1', 'something2', 'something17', 'something24', 'something29', 'something32']

    """
    from warnings import warn
    import re
    if not _nsre:
        _nsre = re.compile('([0-9]+)')
    try:
        return [int(text) if text.isdigit() else f"{text}".lower() for text in re.split(_nsre, in_list)]
    except TypeError:
        if not isinstance(in_list, int):
            warn(
                    f"Yoda_powers::sort_human : element '{in_list}' on the list not understand so don't sort this element\n",
                    SyntaxWarning, stacklevel=2)
            return in_list


def readable_dir(prospective_dir):
    """
    'Type' for argparse - checks that directory exists and if readable, then return the absolute path as PosixPath() with pathlib

    Notes:
        function need modules:

        - pathlib
        - argparse


    Arguments:
        prospective_dir (str): a path to existent path

    Returns:
        :class:`PosixPath`: ``Path(path).resolve()``

    Raises:
         ArgumentTypeError: If directory `path` does not exist.
         ArgumentTypeError: If `path` is not a valid directory.

    Examples:
        >>> import argparse
        >>> parser = argparse.ArgumentParser(prog='test.py', description='''This is demo''')
        >>> parser.add_argument('-f', '--file', metavar="<path/to/file>",type=readable_dir, required=True,
            dest='path_file', help='path to file')
    """
    from argparse import ArgumentTypeError
    from pathlib import Path
    import os

    if not Path(prospective_dir).exists():
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise ArgumentTypeError(f'ERROR: directory "{prospective_dir}" does not exist')
    elif not Path(prospective_dir).is_dir():
        raise ArgumentTypeError(f'ERROR: "{prospective_dir}" is not a valid directory')
    elif not os.access(prospective_dir, os.R_OK):
        raise ArgumentTypeError(f'ERROR: "{prospective_dir}" is not a readable dir')

    return Path(prospective_dir).resolve()


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
    import re
    return re.sub('|'.join(re.escape(key) for key in repls.keys()), lambda k: repls[k.group(0)], str)


def load_in_list(filename):
    """
    Check file exist and create generator with no line break
    file can be a gzip file with '.gz' extension

    Notes:
        function need modules:

        - pathlib

    Arguments:
        filename (str): a path to existent file

    Yields:
        :class:`str`: generator of rows without line break

    Raises:
         FileNotFoundError: If file `filename` does not exist or not valide file

    Example:
        >>> rows = load_in_list("filename")
        >>> list(rows)
        ["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
    """
    from pathlib import Path
    import gzip
    filename = Path(filename)
    if not filename.exists() or not filename.is_file():
        raise FileNotFoundError(
                f'ERROR: Yoda_powers.toolbox.load_in_list() file "{filename}" {"does not exist" if not filename.exists() else "is not a valid file"}')

    open_fn = gzip.open if filename.suffix == ".gz" else open
    with open_fn(filename, "rt") as file_in:
        return (line.rstrip() for line in file_in.readlines())


def load_in_list_col(filename, col=0, sep="\t"):
    """
    Check file exist and create generator with only selected column with no line break
    file can be a gzip file with '.gz' extension

    Arguments:
        filename (str): a path to existent file
        col (int, default): the selected column (python index). Default=0
        sep (str, default): the string separator. Default="\t"

    Yields:
        :class:`str`: generator of rows without line break

    Raises:
         FileNotFoundError: If file `filename` does not exist or not valide file

    Example:
        >>> rows = load_in_list_col("filename", col=1, sep=";")
        >>> list(rows)
        ["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
    """
    from pathlib import Path
    import gzip
    filename = Path(filename)
    if not filename.exists() or not filename.is_file():
        raise FileNotFoundError(
                f'ERROR: Yoda_powers.toolbox.load_in_list_col() file "{filename}" {"does not exist" if not filename.exists() else "is not a valid file"}')
    open_fn = gzip.open if filename.suffix == ".gz" else open
    with open_fn(filename, "rt") as file_in:
        # yield from (line.rstrip().split(sep)[col] for line in file_in.readlines())
        return (line.rstrip().split(sep)[col] for line in file_in.readlines())


def load_in_dict(filename, sep="\t"):
    """
    Check file exist and return a :class:`dict` with load rows first column is the key and value are other column.
    File can be a gzip file with '.gz' extension.

    Arguments:
        filename (str): a path to existent file
        sep (str, default): the string separator. Default="\t"

    Returns:
        :class:`dict`: a python :class:`dict` of file

    Raises:
         FileNotFoundError: If file `filename` does not exist or not valid file

    Example:
        >>> dico = load_in_dict(filename)
        >>> dico
        {
        "col1",["col2","col3"],
        "indiv1",["valeurcol2","valeurcol3"],
        "indiv2",["valeurcol2","valeurcol3"]
        }
    """
    from pathlib import Path
    import gzip
    filename = Path(filename)
    if not filename.exists() or not filename.is_file():
        raise FileNotFoundError(
                f'ERROR: Yoda_powers.toolbox.load_in_dict() file "{filename}" {"does not exist" if not filename.exists() else "is not a valid file"}')
    dico_out = {}
    open_fn = gzip.open if filename.suffix == ".gz" else open
    with open_fn(filename, "rt") as file_in:
        for line in file_in:
            tab_line = line.rstrip().split(sep)
            if tab_line[0] not in dico_out.keys():
                if len(tab_line[1:]) == 0:
                    dico_out[tab_line[0]] = None
                elif len(tab_line[1:]) == 1:
                    dico_out[tab_line[0]] = tab_line[1]
                else:
                    dico_out[tab_line[0]] = tab_line[1:]
    return dico_out


def load_in_dict_selected(filename, column_key=0, column_value=1, sep="\t"):
    """
    Check file exist and return a :class:`dict` with load rows first column is the key and value are other column.
    File can be a gzip file with '.gz' extension.

    Arguments:
        filename (str): a path to existent file
        column_key (int, default): the index for dict keys (python index). Default=0
        column_value (int, default): the index for dict value (python index). Default=1
        sep (str, default): the string separator. Default="\t"

    Returns:
        :class:`dict`: a python :class:`dict` of file

    Raises:
         FileNotFoundError: If file `filename` does not exist or not valid file
         IndexError: If missing data

    Example:
        >>> dico = load_in_dict(filename)
        >>> dico
        {
        "col1",["col2","col3"],
        "indiv1",["valeurcol2","valeurcol3"],
        "indiv2",["valeurcol2","valeurcol3"]
        }
    """
    from pathlib import Path
    import gzip
    filename = Path(filename)
    if not filename.exists() or not filename.is_file():
        raise FileNotFoundError(
                f'ERROR: Yoda_powers.toolbox.load_in_dict() file "{filename}" {"does not exist" if not filename.exists() else "is not a valid file"}')
    dico_out = {}
    open_fn = gzip.open if filename.suffix == ".gz" else open
    try:
        with open_fn(filename, "rt") as file_in:
            for num_line, line in enumerate(file_in):
                tab_line = line.rstrip().split(sep)
                if tab_line[column_key] not in dico_out.keys():
                    dico_out[tab_line[column_key]] = tab_line[column_value]
    except IndexError:
        raise IndexError(
                f'ERROR: Yoda_powers.toolbox.load_in_dict_selected() please check line {num_line + 1}, no value for column {column_value + 1}')

    return dico_out


def load_in_dict_dict(filename, sep="\t"):
    """
    Check file exist and return a :class:`dict` with load rows first column is the key and value are other column.
    File can be a gzip file with '.gz' extension.

    Arguments:
        filename (str): a path to existent file
        sep (str, default): the string separator. Default="\t"

    Returns:
        :class:`dict`: a python :class:`dict` of file

    Raises:
         FileNotFoundError: If file `filename` does not exist or not valid file
         IndexError: If missing data

    Example:
        >>> dico = load_in_dict_dict(filename)
        >>> dico
        {
        "indiv1",{"headerCol2":"toto","headerCol3":"tata"},
        "indiv2",{"headerCol2":"tutu","headerCol3":"titi"},
        "indiv3",{"headerCol2":"tete","headerCol3":"tyty"},
        }
    """
    from pathlib import Path
    import gzip
    from collections import defaultdict, OrderedDict
    filename = Path(filename)
    if not filename.exists() or not filename.is_file():
        raise FileNotFoundError(
                f'ERROR: Yoda_powers.toolbox.load_in_dict_dict() file "{filename}" {"does not exist" if not filename.exists() else "is not a valid file"}')
    dico_out = defaultdict(OrderedDict)
    open_fn = gzip.open if filename.suffix == ".gz" else open
    try:
        with open_fn(filename, "rt") as file_in:
            header = file_in.readline().rstrip().split(sep)
            for num_line, line in enumerate(file_in):
                tab_line = line.rstrip().split(sep)
                if tab_line[0] not in dico_out.keys():
                    for index, head in enumerate(header[1:]):
                        dico_out[tab_line[0]][head] = tab_line[index + 1]
    except IndexError:
        raise IndexError(
                f'ERROR: Yoda_powers.toolbox.load_in_dict_dict() please check line {num_line + 1}, no value for column {head}')
    return dico_out


#################################################
# CLASS
#################################################

class PrintCol:
    """
    Classe qui ajoute des méthodes à print pour afficher de la couleur

    Example:

    >>> PrintCol.red("j'affiche en rouge")
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
        print(f"{cls.__RED}{s}{cls.__END}")

    @classmethod
    def green(cls, s):
        print(f"{cls.__GREEN}{s}{cls.__END}")

    @classmethod
    def yellow(cls, s):
        print(f"{cls.__YELLOW}{s}{cls.__END}")

    @classmethod
    def lightPurple(cls, s):
        print(f"{cls.__LIGHT_PURPLE}{s}{cls.__END}")

    @classmethod
    def purple(cls, s):
        print(f"{cls.__PURPLE}{s}{cls.__END}")


class AutoVivification(dict):
    """
    Implementation of perl's autovivification feature.

    Example:

    >>> a = AutoVivification()
    >>> a[1][2][3] = 4
    >>> a[1][3][3] = 5
    >>> a[1][2]['test'] = 6
    >>> print(a)
    >>> {1: {2: {'test': 6, 3: 4}, 3: {3: 5}}}

    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class Directory(PosixPath):
    """
    Class which derives from PosixPath.
    Checks that the string is and path to valid directory
    add function like list all files/dirs

    Example:
        >>> dir = Directory("./")
        >>> print(dir)
        >>> print(dir.list_files)
        >>> for file in dir.list_files_ext([".py"]):
        >>>     print(file)
    """

    def __init__(self, path_directory=None):
        """
        Arguments:
            path_directory (str): a path to directory
        """
        from pathlib import Path

        if not Path(path_directory).exists() or not Path(path_directory).is_dir():
            raise NotADirectoryError(
                    f'ERROR: Yoda_powers.toolbox.Directory() directory "{path_directory}" {"does not exist" if not Path(path_directory).exists() else "is not a valid directory"}')

        self.path_directory = Path(path_directory).resolve()
        self.__sep = "\n"
        super().__init__()

    @property
    def list_path(self):
        """Generator of files/directory include on folder"""
        return self.path_directory.glob("*")

    @property
    def list_dir(self):
        """Generator of directory include on folder"""
        return (elm for elm in self.path_directory.glob("*") if elm.is_dir())

    @property
    def list_files(self):
        """Generator of files include on folder"""
        return (elm for elm in self.path_directory.glob("*") if elm.is_file())

    def list_files_ext(self, ext=None):
        """Generator of files with specify extension include on folder

        Arguments:
            ext (list): a list of extension like [".py"]
        Yields:
            :class:`PosixPath`: Generator of files with specify extension include on folder
        """
        if not isinstance(ext, list) or not ext:
            raise ValueError(f'ERROR: Yoda_powers.toolbox.directory.list_files_ext() "ext" must be a list not "{ext}"')
        return (elm for elm in self.path_directory.glob(f"**/*") if (elm.is_file() and elm.suffix in ext))

    def __repr__(self):
        return f"{self.__class__}({self.__dict__})"

    def __str__(self):
        """print format"""
        return f"""
path_directory={self.path_directory}
list_path:\n   - {"   - ".join([f'{elm.name}{self.__sep}' for elm in self.list_path])}
list_dir:\n   - {"   - ".join([f'{elm.name}{self.__sep}' for elm in self.list_dir])}
list_files:\n   - {"   - ".join([f'{elm.name}{self.__sep}' for elm in self.list_files])}
"""
