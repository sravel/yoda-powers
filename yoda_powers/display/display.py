#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################
# Modules
##################################################
# Python modules
from ..toolbox import sort_human

##################################################
# Functions

def dict_2_txt(dico, sep="\t"):
    """
    Function that takes a dictionary and returns a string with separator:

    Arguments:
        dico (:class:`dict`): the python dict to translate to formated string.
        sep (:class:`str`): the separator for join . Defaults to '\\t'.

    Returns:
        str: formated dict to string

    Examples:
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

    return "\n".join([f"{key}{sep}{dico[key]}" for key in sorted(dico.keys(), key=sort_human)])


def dict_dict_2_txt(dico, first="Info", sep="\t"):
    """
    Function that takes a dictionary and returns a tabular string with:

    Arguments:
        dico (:obj:`dict`): the python dict to translate to formated string.
        first (:obj:`str`): the first column header name. Default to 'Info'.
        sep (:obj:`str`): the separator for join. Defaults to '\\t'.

    Returns:
        str: formated dict to string

    Examples:
        >>> dico = {"Souche1":{"NUM":"171","MIN":"2042","MAX":"3133578","N50 BP":"938544","N50 NUM":"11"},
                    "Souche2":{"NUM":"182","MIN":"5004","MAX":"74254","N50 BP":"45245"}}
        >>> dict_dict_2_txt(dico,"souches")
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


def dict_list_2_txt(dico, sep="\t"):
    """
    Function that takes a dictionary of list and returns a tabular string with:

    Arguments:
        dico (:obj:`dict`): the python dict to translate to formated string.
        sep (:obj:`str`): the separator for join . Defaults to '\\t'.

    Returns:
        str: formated dict to string

    Examples:
        >>> dico = {"key1":["value1","value1"], "key2":["value2","value2"],"key3":["value3","value3"]}
        >>> dict_list_2_txt(dico)
        key1	value1	value1
        key2	value2	value2
        key3	value3	value3

    """

    txt_output = ""
    for key in sorted(dico.keys(), key=sort_human):
        txt_output += f"{key}{sep}{sep.join(sorted(dico[key], key=sort_human))}\n"
    return txt_output
