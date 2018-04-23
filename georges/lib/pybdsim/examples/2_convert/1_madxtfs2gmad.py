#! /usr/bin/env python2.7

import pybdsim

def Main():
    pybdsim.Convert.MadxTfs2Gmad("twiss35tevb1_short.tar.gz", "lhcb1")

if __name__ == "__main__":
    Main()
