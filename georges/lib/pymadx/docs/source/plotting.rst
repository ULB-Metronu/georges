========
Plotting
========

The `pymadx.Plot` module provides various plotting utilities.

Plotting Features
-----------------

* Make default optics plots.
* Add a machine lattice to any pre-existing plot.
* Interactive plots with the machine diagram following the mouse zooming.
* Interactive plots with searching for nearest element.

Optics Plots
------------  

A simple optics plot may be made with the following syntax::

  a = pymadx.Data.Tfs("mytwissfile.tar.gz")
  a.Plot()

This creates a plot of the Beta amplitude functions against curvilinear S position. A
colour diagram representing the machine is also produced above the graph as shown below.

.. figure:: figures/optics.pdf
   :width: 90%
   :align: center

The command has optional arguments such as a title string to be put at the top of the graph
and whether to also plot the horizontal dispersion function. This function is provided as
a quick utility and not the ultimate plotting script. The user can make their own plot and
then append a machine diagram at the end if they wish.::

  f = matplotlib.pyplot.figure()
  # user plotting commands here
  pymadx.Plot.AddMachineLatticeToFigure(f, "mytwissfile.tar.gz")

`gcf()` is a matplotlib.pyplot function to get a reference to the current matplotlib
figure and can be used as the first argument.::

  pymadx.Plot.AddMachineLatticeToFigure(gcf(), "mytwissfile.tar.gz")

.. note:: It becomes difficult to adjust the axes and layout of the graph after adding the
	  machine description. It is therefore strongly recommended to do this last.


Colour Coding
-------------

Each magnet is colour coded an positioned depending on its type and strength.

+--------------+------------------+--------------+-----------------------------------------------+
| **Type**     | **Shape**        | **Colour**   | **Vertical Position**                         |
+==============+==================+==============+===============================================+
| drift        | N/A              | Not shown    | N/A                                           |
+--------------+------------------+--------------+-----------------------------------------------+
| sbend        | Rectangle        | Blue         | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| rbend        | Rectangle        | Blue         | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| hkicker      | Rectangle        | Purple       | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| vkicker      | Rectangle        | Pink         | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| quadrupole   | Rectangle        | Red          | Top half for K1L > 0; Bottom half for K1L < 0 |
+--------------+------------------+--------------+-----------------------------------------------+
| sextupole    | Hexagon          | Yellow       | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| octupole     | Hexagon          | Green        | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| multiple     | Hexagon          | Light grey   | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| rcollimator  | Rectangle        | Black        | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| ecollimator  | Rectangle        | Black        | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+
| *any other*  | Rectangle / Line | Light Grey   | Central always                                |
+--------------+------------------+--------------+-----------------------------------------------+

.. note:: In all cases if the element is a magnet and the appropriate strength is zero, it is
	  shown as a grey line.

Plot Interactivity
------------------

With the adition of the machine axes, some extra interactivity is included to the matplotlib
figures.

 * zooming - if the 'right-click and drag' zoom feature is used on the machine diagram, the
   graph will automatically update and follow the machine diagram.
 * xlim - setting the horizontal graph limits with the 'xlim' command will update both the
   machine diagram and the graph.
 * querying - right-clicking anywhere on the graph will print out the name of the nearest element
   in the terminal.

