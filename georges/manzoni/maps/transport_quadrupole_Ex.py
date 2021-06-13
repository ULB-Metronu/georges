from numba import njit 
import numpy as np 
from numpy import cos, sin, cosh, sinh, sqrt 

def compute_transport_quadrupole_Ex_matrix(element_parameters: list) -> np.ndarray:
    L: float = element_parameters[0] 
    K1: float = element_parameters[1] 
    d: float = element_parameters[2] 
    K1 = K1/(1+d) 
    R = np.zeros((6,6))
    R[4,4] = 1
    R[5,5] = 1
    if (K1 == 0):
        R[0,0] = 1 
        R[0,1] = L 
        R[1,1] = 1 
        R[2,2] = 1 
        R[2,3] = L 
        R[3,3] = 1 
    if (K1 > 0):
        R[0,0] = cos(sqrt(K1)*L) 
        R[0,1] = sin(sqrt(K1)*L)/sqrt(K1) 
        R[1,0] = -sqrt(K1)*sin(sqrt(K1)*L) 
        R[1,1] = cos(sqrt(K1)*L) 
        R[2,2] = cosh(sqrt(K1)*L) 
        R[2,3] = sinh(sqrt(K1)*L)/sqrt(K1) 
        R[3,2] = sqrt(K1)*sinh(sqrt(K1)*L) 
        R[3,3] = cosh(sqrt(K1)*L) 
    if (K1 < 0):
        R[0,0] = cosh(L*sqrt(-K1)) 
        R[0,1] = sinh(L*sqrt(-K1))/sqrt(-K1) 
        R[1,0] = sqrt(-K1)*sinh(L*sqrt(-K1)) 
        R[1,1] = cosh(L*sqrt(-K1)) 
        R[2,2] = cos(L*sqrt(-K1)) 
        R[2,3] = sin(L*sqrt(-K1))/sqrt(-K1) 
        R[3,2] = -sqrt(-K1)*sin(L*sqrt(-K1)) 
        R[3,3] = cos(L*sqrt(-K1)) 
    return R 

def compute_transport_quadrupole_Ex_tensor(element_parameters: list) -> np.ndarray:
    L: float = element_parameters[0] 
    K1: float = element_parameters[1] 
    d: float = element_parameters[2] 
    K1 = K1/(1+d) 
    T = np.zeros((6,6,6))
    if (K1 > 0):
        T[0,0,4] = sqrt(K1)*L*sin(sqrt(K1)*L)/2 
        T[0,1,4] = -L*cos(sqrt(K1)*L)/2 + sin(sqrt(K1)*L)/(2*sqrt(K1)) 
        T[1,0,4] = sqrt(K1)*sin(sqrt(K1)*L)/2 + K1*L*cos(sqrt(K1)*L)/2 
        T[1,1,4] = sqrt(K1)*L*sin(sqrt(K1)*L)/2 
        T[2,2,4] = -sqrt(K1)*L*sinh(sqrt(K1)*L)/2 
        T[2,3,4] = -L*cosh(sqrt(K1)*L)/2 + sinh(sqrt(K1)*L)/(2*sqrt(K1)) 
        T[3,2,4] = -sqrt(K1)*sinh(sqrt(K1)*L)/2 - K1*L*cosh(sqrt(K1)*L)/2 
        T[3,3,4] = -sqrt(K1)*L*sinh(sqrt(K1)*L)/2 
    if (K1 < 0):
        T[0,0,4] = -L*sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[0,1,4] = -(-L*sqrt(-K1)*cosh(L*sqrt(-K1)) + sinh(L*sqrt(-K1)))/(2*sqrt(-K1)) 
        T[1,0,4] = K1*L*cosh(L*sqrt(-K1))/2 - sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[1,1,4] = L*sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[2,2,4] = L*sqrt(-K1)*sin(L*sqrt(-K1))/2 
        T[2,3,4] = L*cos(L*sqrt(-K1))/2 - sin(L*sqrt(-K1))/(2*sqrt(-K1)) 
        T[3,2,4] = -K1*L*cos(L*sqrt(-K1))/2 + sqrt(-K1)*sin(L*sqrt(-K1))/2 
        T[3,3,4] = -L*sqrt(-K1)*sin(L*sqrt(-K1))/2 
    return T 

