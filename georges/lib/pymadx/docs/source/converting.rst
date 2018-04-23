==========
Conversion
==========

pymadx provides the ability to convert a MADX model into another format. These
are detailed in the sections below.

For conversion of MADX to BDSIM GMAD format, please see the pybdsim documentation
`<http://www.pp.rhul.ac.uk/bdsim/pybdsim/convert.html#madxtfs2gmad>`_.

Mad8ToMadx
==========

A MAD8 model can be converted to MADX. This relies on the pymad8 package
and the conversion is performed there.


TfsToPtc
========

A MADX model as described by a full twiss table in a TFS file can be
prepared into a new MADX / PTC model suitable for tracking with PTC
immediately.
