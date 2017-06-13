import os
import sys
print(sys.path)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import beamline
import beam
import madx