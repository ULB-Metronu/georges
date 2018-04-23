import pybdsim
import subprocess


def adjust_collimator(colName, colParam, colVal, colFile, fileToBeConverted,
                      outputFile):
    collimator_obj = pybdsim.Convert.Mad8SavelineCollimatorDatabase(colFile)
    collimator_dict = collimator_obj.GetDict()
    
    if colName in collimator_dict:
        collimator_dict[colParam] = colVal
        outputFile = '%s_%s_%s_%s.gmad' % (outputFile, colName, colParam,  colVal)
    
    convert_collimator(fileToBeConverted, outputFile, collimator_dict)
            

def convert_collimator(fileToBeConverted, outputFile, collimatorSettings):
    pybdsim.Convert.Mad8Saveline2Gmad(fileToBeConverted, outputFile, 
                                       collimator_dict=collimatorSettings, 
                                       beam=True, optics=True, loss=True)
    
    subprocess.call(['bdsim', '--file=%s' % outputFile, '--output=root',
                     '--outfile=%s%s' % (outputFile.split('.')[0], outputFile.split('.')[1]),
                     '--seed=1', '--batch'])
    
