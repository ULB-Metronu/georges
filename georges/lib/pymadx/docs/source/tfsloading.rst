===============================
TFS File Loading & Manipulation
===============================

MADX outputs Twiss information as well as PTC tracking data in their own Table
File System (TFS). This is the only format used by MADX. pymadx includes a class
called Tfs for the purpose of loading and manipulating this data.

The TFS format is described in the MADX manual available from `the madx website <http://madx.web.cern.ch>`_.
The format
roughly is described as a text file. The file contains first a header with key and
value pairs for one-off definitions. This is proceeded by a line
with column names and a line with the data type of each column. After this each line
typically represents the values of the lattice for a particular element with each new
line containing the values at a subsequent element in the lattice. We maintain the
concept of this table and refer to 'rows' and 'columns'.

Tfs Class Features
------------------

* Loading of TFS files.
* Loading of TFS files compressed and archived with .tar.gz suffix without decompressing.
* Report a count of all different element types.
* Get a particular column.
* Get a particular row.
* Get elements of a particular type.
* Get a numerical index from the name of the element.
* Find the curvilinear S coordinate of an element by name.
* Find the name of the nearest element at a given S coordinate.
* Plot an optics diagram.
* Roll a lattice to start from a different point.
* Calculate a beam size given the Twiss parameters, dispersion and emittance (in the header).
* Determining whether a given component perturbs the beam.
* Extract a 'segment' if PTC data is present.
* Slice a lattice (in the Python sense) with new S coordinates.


Loading
-------

A file may be loading by constructing a Tfs instance from a file name.

>>> import pymadx
>>> a = pymadx.Data.Tfs("myTwissFile.tfs")

.. note:: The import will be assumed from now on in examples.

A file compressed using tar and gzip may also be loaded without first uncompressing
without any difference in functionality. Not temporary files are created::

  tar -czf myTwissFile.tar.gz myTwissFile.fs
  
>>> import pymadx
>>> a = pymadx.Data.Tfs("myTwissFile.tar.gz")

.. note:: The detection of a compressed file is based on 'tar' or 'gz' existing
	  in the file name.

Twiss File Preparation
----------------------

You may export twiss data from MADX with a choice of columns. We often find it beneficial
to not specify any columns at all, which results in all available columns being written.
This large number (~70) makes the file less human-readable but ensures no information is
omitted. Such an export will also increase the file size, however, we recommend compressing
the file with tar and gzip as the ASCII files compress very well with a typically compression
ratio of over 10:1.

The following MADX syntax in a MADX input script will prepare a Tfs file with all columns where
"SEQUENCENAME" is the name of the sequence in MADX.::

  select,flag=twiss, clear; 
  twiss,sequence=SEQUENCENAME, file=outputfilename.tfs;


Querying
--------

The Tfs class can be used to query the data in various ways.

Basic Information
*****************

 * All data is stored in the **data** object inside the class
 * The header information is stored in **header**.
 * The names of all elements in order is stored in **sequence**.
 * The names of all columns in the file is stored in **columns**

Generally, members beginning with small letters are objects and capital letters are functions.

A nice summary of the file can be provided with the `ReportPopulations` function.::

  a = pymadx.Data.Tfs("mytwissfile.tar.gz")
  a.ReportPopulations()

  Filename > twiss_v5.2.tfs
  Total number of items > 1032
  Type........... Population
  MULTIPOLE...... 516
  DRIFT.......... 201
  QUADRUPOLE..... 102
  MARKER......... 78
  MONITOR........ 64
  SBEND.......... 24
  SEXTUPOLE...... 18
  HKICKER........ 15
  VKICKER........ 14

Indexing and Slicing
********************

The instance may be indexed like a normal Python iterable structure such as a list or a tuple.
Square brackets with a number *i* will return the *ith* element in the sequence. A Python 'slice'
may also be used where a range of elements may be selected. If only one element is indexed a
Python dictionary is returned for that element. If a range is required, another Tfs instance
is returned::

  a = pymadx.Data.Tfs("mytwissfile.tar.gz")
  a[3]         # 4th element in sequence (0,1,2,3!)
  a[3:10]      # 4th to 11th elements (tfs instance returned)
  a[3:10:-1]   # similarly but in steps on -1 ie reversed
  a['IP1':300] # find element named IP1 (exactly) and start from that until the #301th element
  a['IP3':]    # find element named IP3 (exactly) and take from there to the end of the file
  a['L230A']   # returns a Python dictionary for element named L230A

If you know the name of an element you can search for it and get the index from that.::

  a.IndexFromName('L230A')
  >>> 995

You can also search by nearest curvilinear S coordinate along the beam line.::

  a.IndexFromNearestS(34.4)
  >>> 225
  a[225]['NAME']

Row or Element
**************

A row of data is an entry for a particular element. The Tfs class is conceptually a list of
elements. Each element is represented by a Python dictionary that has a key for each column.
The list of acceptable keys (ie names of columns) can be found in the member named 'colums'.::

  a.columns #prints out list of column names

If a single element is indexed, a dictionary is returned and can be accessed - even in one step.::

  d = a[123]
  d['NAME']
  >>> 'MQD8X'
  a[123]['NAME'] # equivalent


Looping & Iterating
*******************

The Tfs class may be iterated over like a list in Python. For each iteration a dictionary
for that element is returned.::
  
  for el in a:
      print(el['NAME'])

Beam Sizes
**********

For convenience the beam size is calculated from the Beta amplitude functions, the emittance
and dispersion if they are present. The emittance is defined by 'EX' and 'EY' in the header.
These are calculated according to

.. math::

   \sigma_x &= \sqrt{ \beta_x \epsilon_x + D(S)^2 \frac{\sigma_{E}^{2}}{\beta_{\mathrm{Lorentz}}^{2}}} \\
   \sigma_y &= \sqrt{ \beta_y \epsilon_y + D(S)^2 \frac{\sigma_{E}^{2}}{\beta_{\mathrm{Lorentz}}^{2}}}

:math:`\sigma_E` in MADX is fractional. Here we use the relation

.. math::

   \sigma_E = \frac{\Delta E}{E} = \beta_{\mathrm{Lorentz}}^{2} \frac{\Delta p}{p}

.. note:: MADX input files often don't have a sensible emittance defined as it is not always
	  required. Ensure the emittance is what you intended it to be in the Tfs file.


Modification
************

It is not recommended to modify the data structures inside the Tfs class. Of course one can,
but one must be careful of Python's copying behaviour. Often a 'deep copy' is required or
care must be taken to modify the original and not a reference to a particular variable.


