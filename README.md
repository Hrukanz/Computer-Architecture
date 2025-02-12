ENCE461 Assignment 2 Template
=============================

See assignment instructions [here](doc/instructions/instructions.pdf)

Contents
 - `doc/` - assignment instructions, lab notes, report template.
 - `reference/` - correct output for test comparison.
 - `poisson.c` - basic template to work from. Write your solution here.
 - `threads.c` - example on how to use POSIX thread library.
 - `test.sh` - automatic testing script.

Building
--------

Build instructions are provided in the Makefile.  Use make to build.


Testing
-------

Run `make test`

It will automatically run your solution for three cube sizes and compare the
output against some correct reference files. **Do not edit these reference
files!**
