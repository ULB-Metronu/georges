#! /usr/bin/env python2.7

import pybdsim

a = pybdsim.Builder.Machine()

a.AddDrift('d1',1.2,apertureType="rectangular", aper1=(3,"cm"), aper2=(1.5,"cm"))
a.AddQuadrupole("q1",0.2,k1=0.034515)
a.AddDrift('d2',0.5)
a.AddQuadrupole("q2",0.2,k1=-0.034515)
a.AddDrift('d3',0.76)
a.AddDipole('dp','sbend',4.5,0.0789)
a.AddDrift('d4',0.6)

a.AddSampler("all")

a.Write("testlattice")
