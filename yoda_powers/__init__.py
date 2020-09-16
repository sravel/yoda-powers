#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from . import bio
from . import display
from . import toolbox
from . import scripts

##################################################
# Variables Globals
__version__ = "1.0.0"

__doc__ = """
:author: Sebastien Ravel
:contact: sebastien.ravel@cirad.fr
:date: 22-07-2020
:version: """ + __version__ + """

Use it to import very handy functions.

Warnings:
    This module run with |PythonVersions|


Install
=======

Required Module install
-----------------------

Modules ``BioPython`` will be install with :class:``yoda_powers``:

Global install
^^^^^^^^^^^^^^

::

    # for all users (require root privilege)
    sudo pip3 install yoda_powers

    # for own
    pip3 install yoda_powers --user

Developing version
^^^^^^^^^^^^^^^^^^

If you want to use an **unofficial version** of the ``yoda_powers`` module, you need to work from
a clone of this ``git`` repository.
following actions:

1. Clone from github ::

    $ git clone https://github.com/sravel/yoda-powers.git

2. Go in the cloned directory ::

    $ cd yoda-powers

3. Install ``yoda-powers`` in editable mode ::

    $ sudo pip3 install -e .

Sub module information
======================

This module are split on three sub-module

- ``yoda_powers.toolbox``: handy functions common for many scripts to easy check file/directory, load file info, ...
- ``yoda_powers.display``: handy functions to display/write python object like dict of dict.
- ``yoda_powers.bio``: handy functions to manipulate biological data with Biopython.



Example of usage
================

Example:

>>> from yoda_powers.display import dict_2_txt
>>> dico = {'key1':'value1','key2':'value2','key3':'value3'}
>>> dict_2_txt(dico)
key1	value1
key2	value2
key3	value3


more info at https://yoda-powers.readthedocs.io/en/latest/

.. |PythonVersions| image:: https://img.shields.io/badge/python-3.6+-blue.svg
   :target: https://www.python.org/downloads
   :alt: Python /3.6+
"""
