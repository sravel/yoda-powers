#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys
from setuptools import setup, find_packages

HERE = Path(__file__).resolve().parent
sys.path.insert(0, HERE.stem)
sys.path.insert(0, HERE.joinpath("/yoda_powers/").stem)
import yoda_powers

# auto build console_scripts entry
from yoda_powers import scripts

console_scripts_list = []
for script in scripts.list_scripts:
    if "__init__" not in script.name:
        txt = f"{script.name} = yoda_powers.scripts.{script.stem}:main"
        console_scripts_list.append(txt)


########################################################################
# main function
def main():
    setup(
            name=yoda_powers.__name__,
            version=yoda_powers.__version__,
            description=yoda_powers.__doc__,
            long_description_content_type='text/x-rst',
            long_description=HERE.joinpath('README.rst').open(encoding='utf-8').read(),
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
            platforms=['cross-platform', 'Windows i686', 'MacOSX AMD64'],
            keywords=[
                    'python',
            ],
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
                    yoda_powers.__name__: [yoda_powers.__name__ + " = __init__"],
                    'console_scripts'   : console_scripts_list
            },
    )


if __name__ == '__main__':
    main()
