fitsstamp
=========
Python module and command-line utilities for cutting and pasting small "stamp" 
images from astronomical FITS image files. Stamp cutting locations can be 
defined in one of two ways:

- SAOImage DS9 region files, including composite regions consisting of multiple 
  included/excluded parts by giving multiple regions a common text label
- A FITS-format label image where each individual object's pixels are flagged 
  with a unique non-zero integer (e.g. sextractor segmentation map)

Several additional options are available, including the ability to blank out 
non-object pixels by setting them to a fixed value.

A module function and command-line utility are also included for pasting stamps 
back into their source image. This can be used, for example, to cut out 
individual objects, model them out with software such as galfit, then paste the 
"cleaned" sections back into the original image.

Dependencies and Installation
-----------------------------
The module relies upon:

- numpy
- [astropy](http://www.astropy.org/), for FITS file I/O and WCS transforms
- [pyregion](http://leejjoon.github.io/pyregion/users/overview.html), for 
  parsing and manipulating SAOImage DS9 region files

To install, simply run:
    
    python setup.py install

To install in a non-standard location (e.g. a Dropbox folder)
    
    python setup.py install --prefix=~/Dropbox/Python
