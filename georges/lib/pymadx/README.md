#pymadx#

A python package containing both utilities for processing and analysing MADX output.

## Authors ##

L. Nevay
S. Boogert
S. Walker
A. Abramov
W. Shields
J. Snuverink

## Setup ##

From within the pymadx root directory:

$ make install

or for development:

$ make develop


```
#!python

$>python
$>>> import pymadx
$>>> t = pymadx.Data.Tfs("twiss.tfs")
```

## Dependencies ##

 * matplotlib
 * numpy