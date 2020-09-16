=====================================
Welcome to Yoda powers documentation!
=====================================

.. image:: ./docs/source/SupplementaryFiles/yoda.png
   :target: ./docs/source/SupplementaryFiles/yoda.png
   :alt: yoda Logo

|PythonVersions| |PypiPackage|


Table of Contents
=================

.. contents::
   :depth: 2


About this package
==================

Use it to import very handy functions.

**Documentation available at :** `<https://yoda-powers.readthedocs.io/en/latest>`_


Quick manual
============

Install
-------

Global install
^^^^^^^^^^^^^^
.. warning::
    not working yet, module in development

To compile and install yoda_powers you should do this::

    # for all users (requiring super-user rights)
    sudo pip install yoda_powers
    # for the current user only
    pip install yoda_powers --user

Install for a specific version of python for example python3.7::

    # for all users (requiring super-user rights)
    sudo python3.7 -m pip install yoda_powers
    # for the current user only
    python3.7 -m pip install yoda_powers --user

Developing version
^^^^^^^^^^^^^^^^^^
If you want to use an **unofficial version** of the ``yoda_powers`` module, you need to work from
a clone of this ``git`` repository.
following actions:

1. Clone from github ::

    git clone https://github.com/sravel/yoda-powers.git
    
2. Go in the cloned directory ::

    cd yoda_powers

3. Install in editable mode ::

    pip install -e . --user


Citation
========

:author: Sebastien Ravel
:contact: sebastien.ravel@cirad.fr

Please cite the github module



License
=======

Licencied under GPLv3
Intellectual property belongs to CIRAD and authors.


.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7+-blue.svg
   :target: https://www.python.org/downloads
   :alt: Python /3.7+

.. |PypiPackage| image:: https://badge.fury.io/py/Yoda_powers.svg
   :target: https://pypi.org/project/Yoda_powers
   :alt: PyPi package
