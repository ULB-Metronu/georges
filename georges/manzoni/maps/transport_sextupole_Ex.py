from numba import njit 
import numpy as np 
from numpy import cos, sin, cosh, sinh, sqrt 
from numba.typed import List as nList

def compute_transport_sextupole_Ex_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0] 
    K2: float = element_parameters[1] 
    d: float = element_parameters[2] 
    K2 + K2/(1+d) 
    R = np.zeros(6,6)
    
    R[0,0] = 1 
    R[0,1] = L 
    R[1,1] = 1 
    R[2,2] = 1 
    R[2,3] = L 
    R[3,3] = 1 
    R[4,4] = 1
    R[5,5] = 1
    return R 

def compute_transport_sextupole_Ex_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0] 
    K2: float = element_parameters[1] 
    d: float = element_parameters[2] 
    K2 + K2/(1+d) 
    T = np.zeros(6,6,6)
    T[0,0,0] = -K2*L**2/2 
    T[0,0,1] = -K2*L**3/3 
    T[0,1,1] = -K2*L**4/12 
    T[0,2,2] = K2*L**2/2 
    T[0,2,3] = K2*L**3/3 
    T[0,3,3] = K2*L**4/12 
    T[1,0,0] = -K2*L 
    T[1,0,1] = -K2*L**2 
    T[1,1,1] = -K2*L**3/3 
    T[1,2,2] = K2*L 
    T[1,2,3] = K2*L**2 
    T[1,3,3] = K2*L**3/3 
    T[2,0,2] = K2*L**2 
    T[2,0,3] = K2*L**3/3 
    T[2,1,2] = K2*L**3/3 
    T[2,1,3] = K2*L**4/6 
    T[3,0,2] = 2*K2*L 
    T[3,0,3] = K2*L**2 
    T[3,1,2] = K2*L**2 
    T[3,1,3] = 2*K2*L**3/3 
    return T 
