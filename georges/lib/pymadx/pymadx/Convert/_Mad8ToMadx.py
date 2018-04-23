#wrapper for clean interface
#wrapper functions should be placed in here and long complicated 

def Mad8ToMadx(inputfilename):
    """
    Convert a Mad8 file to the syntax of a MadX file.

    The output file name will be the input file name with
    the file extension replaced by .xsifx
    """
    _Mad8.Mad8ToMadX(inputfilename)
