#pyBDSIM#

A python package containing both utilities for preparing and analysing BDSIM input and output as well as controlling BDSIM.

## Authors ##

L. Nevay
S. Boogert

## Setup ##

From within the pybdsim root directory:

$ make install

or for development:

$ make develop


```
#!python

$>python
$>>> import pybdsim
$>>> a = pybdsim.Data.Load("run1_output.txt")
$>>> hist(a.Xp())
```
