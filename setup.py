import os
import re
from warnings import warn
from subprocess import check_output, CalledProcessError
from setuptools import setup

_packages = ['fitsstamp']
_scripts = ['scripts/fitsstamp_cut', 'scripts/fitsstamp_paste']
_requires = ['numpy', 'astropy', 'pyregion']

def _version_from_init():
    main_init = os.path.join(_packages[0], '__init__.py')
    ver_pat = re.compile('__version__\s*=\s*\'([^\']*)\'')
    with open(main_init) as initfile:
        for line in initfile.readlines():
            match = ver_pat.match(line)
            if match is not None:
                return match.group(1)
    return ''

_version = _version_from_init()

if os.path.exists('.git'):
    try:
        _version = check_output(['git', 'describe'])
    except CalledProcessError:
        warn('git versioning information found but call to git describe '
             'failed. Version number will not include commit number and hash '
             'for development builds.')

setup(name=_packages[0],
      version=_version,
      author='Matt Mechtley',
      author_email='matt@cherem.org',
      url='https://github.com/mmechtley/fitsstamp',
      license='BSD 2-clause',
      description='Utilities for cutting and pasting small "stamp" images from '
                  'astronomical FITS image files.',
      packages=_packages,
      scripts=_scripts,
      requires=_requires,
      zip_safe=True)