#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from io import open
import os
from setuptools import setup, find_packages

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, HERE + "/yoda_powers/")
import yoda_powers


########################################################################
# main function
def main():
    setup(
            name=yoda_powers.__name__,
            version=yoda_powers.__version__,
            description=yoda_powers.__doc__,
            long_description=open(os.path.join(HERE, 'README.rst'), encoding='utf-8').read(),
            long_description_content_type='text/x-rst',
            # these are optional and override conf.py settings
            command_options={
                    'build_sphinx': {
                            'project'   : ('setup.py', yoda_powers.__name__),
                            'version'   : ('setup.py', yoda_powers.__version__),
                            'release'   : ('setup.py', yoda_powers.__version__),
                            'source_dir': ('setup.py', './docs/source'),
                            'build_dir' : ('setup.py', './docs/build'),

                    }},
            classifiers=[
                    'Development Status :: 5 - Production/Stable',
                    'Environment :: Other Environment',
                    'Intended Audience :: Developers',
                    'Intended Audience :: End Users/Desktop',
                    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                    'Operating System :: OS Independent',
                    'Programming Language :: Python :: 3.7',
                    'Natural Language :: English',
            ],
            author="Sébastien Ravel",
            url="https://github.com/sravel/yoda-powers",
            download_url="https://github.com/sravel/yoda-powers/archive/{}.tar.gz".format(yoda_powers.__version__),
            license='LGPL license',
            platforms=['cross-platform', 'Windows i686', 'MacOSX AMD64'],
            keywords=[
                    'python',
            ],
            package_dir=["yoda_powers","yoda_powers"],
            packages=find_packages(),
            package_data={
                    'yoda_powers': ['*.ini'],
            },
            include_package_data=True,
            install_requires=[
                    'BioPython',
            ],
            options={
                    'bdist_wheel':
                        {'universal': True}
            },
            # scripts=['./yoda_powers/scripts/'],
            zip_safe=False,  # Don't install the lib as an .egg zipfile
            entry_points={
                    'yoda_powers'    : ["yoda_powers = __init__"],
                    'console_scripts': ["cli.py = yoda_powers.scripts.cli:main", ]
            },
    )


if __name__ == '__main__':
    main()
