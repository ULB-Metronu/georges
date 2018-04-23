import platform as _platform
import ROOT as _ROOT

# Check for which libraries to load
if _platform.platform().find("Darwin") != -1 :
    _ROOT.gSystem.Load("libbdsimRootEvent.dylib")
elif _platform.platform().find("Linux") != -1:
    _ROOT.gSystem.Load("libbdsimRootEvent.so")
else :
    pass
    # TODO issue warning on unknown system

from Options   import *
from Model     import *
from Run       import *
from Event     import *
from Processed import *
from Root      import *
# from Plots     import *