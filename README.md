Overview
========

NEAT is a program for calculating chemical abundances in photoionised nebulae.  It can propagate uncertainties and compensate for measurement biases affecting weak lines.

The code includes the necessary atomic data files, but you might also wish to use some different atomic data sets which are available separately in the Atomic-data repository.

http://github.com/rwesson/Atomic-data

Installation
============

To compile the code, type

  make

in the directory the files are in.  Our makefile assumes that you have gfortran installed.  Other compilers should work as well but may require some tweaks to the code.  Assuming that worked fine, then install the code by typing

  sudo make install

If you don't have root access, or just want to install the code to somewhere other than the standard linux directories, then you should edit the makefile and changed the value of the variable PREFIX accordingly.

Manual
======

Full documentation is available on the NEAT website:

http://www.nebulousresearch.org/codes/neat

License
=======

NEAT is released under the GNU General Public License, v3.  The full text is in the file LICENSE.
