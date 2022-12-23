.. kernels:

Kernels
-------

The file `kernels.py` contains the loops that are the core of the particles propagation based on
their coordinates. Different batches are available, to allow a matrix (order 1) propagation,
a tensor (order 2) propagation or a matrix followed by a tensor (orders 1+2) propagations.