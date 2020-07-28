#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from . import bio
from . import display
from . import toolbox

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

"""
