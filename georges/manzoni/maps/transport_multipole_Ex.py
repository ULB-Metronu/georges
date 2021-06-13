from numba import njit 
import numpy as np 
from numpy import cos, sin, cosh, sinh, sqrt
from numba.typed import List as nList

def compute_transport_multipole_Ex_matrix(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0] 
    K1: float = element_parameters[1] 
    d: float = element_parameters[3] 
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

def compute_transport_multipole_Ex_tensor(element_parameters: nList) -> np.ndarray:
    L: float = element_parameters[0] 
    K1: float = element_parameters[1] 
    K2: float = element_parameters[2] 
    d: float = element_parameters[3] 
    K1 = K1/(1+d) 
    K2 = K2/(1+d) 
    T = np.zeros((6,6,6))
    if (K1 == 0):
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
    if (K1 > 0):
        T[0,0,0] = K2*(cos(sqrt(K1)*L)**2 + cos(sqrt(K1)*L) - 2)/(3*K1) 
        T[0,0,1] = 2*K2*(cos(sqrt(K1)*L) - 1)*sin(sqrt(K1)*L)/(3*K1**(3/2)) 
        T[0,0,4] = sqrt(K1)*L*sin(sqrt(K1)*L)/2 
        T[0,1,1] = K2*(sin(sqrt(K1)*L)**2 + 2*cos(sqrt(K1)*L) - 2)/(3*K1**2) 
        T[0,1,4] = -L*cos(sqrt(K1)*L)/2 + sin(sqrt(K1)*L)/(2*sqrt(K1)) 
        T[0,2,2] = K2*(-3*cos(sqrt(K1)*L) + cosh(sqrt(K1)*L)**2 + 2)/(5*K1) 
        T[0,2,3] = -2*K2*(sin(sqrt(K1)*L) - sinh(2*sqrt(K1)*L)/2)/(5*K1**(3/2)) 
        T[0,3,3] = K2*(2*cos(sqrt(K1)*L) + cosh(sqrt(K1)*L)**2 - 3)/(5*K1**2) 
        T[1,0,0] = -K2*(sin(sqrt(K1)*L) + sin(2*sqrt(K1)*L))/(3*sqrt(K1)) 
        T[1,0,1] = 2*K2*(-cos(sqrt(K1)*L) + cos(2*sqrt(K1)*L))/(3*K1) 
        T[1,0,4] = sqrt(K1)*sin(sqrt(K1)*L)/2 + K1*L*cos(sqrt(K1)*L)/2 
        T[1,1,1] = 2*K2*(cos(sqrt(K1)*L) - 1)*sin(sqrt(K1)*L)/(3*K1**(3/2)) 
        T[1,1,4] = sqrt(K1)*L*sin(sqrt(K1)*L)/2 
        T[1,2,2] = K2*(3*sin(sqrt(K1)*L) + sinh(2*sqrt(K1)*L))/(5*sqrt(K1)) 
        T[1,2,3] = -2*K2*(cos(sqrt(K1)*L) - cosh(2*sqrt(K1)*L))/(5*K1) 
        T[1,3,3] = -2*K2*(sin(sqrt(K1)*L) - sinh(2*sqrt(K1)*L)/2)/(5*K1**(3/2)) 
        T[2,0,2] = 2*K2*(2*sin(sqrt(K1)*L)*sinh(sqrt(K1)*L) - cos(sqrt(K1)*L)*cosh(sqrt(K1)*L) + cosh(sqrt(K1)*L))/(5*K1) 
        T[2,0,3] = -2*K2*(-2*sin(sqrt(K1)*L)*cosh(sqrt(K1)*L) + cos(sqrt(K1)*L)*sinh(sqrt(K1)*L) + sinh(sqrt(K1)*L))/(5*K1**(3/2)) 
        T[2,1,2] = -2*K2*(sin(sqrt(K1)*L)*cosh(sqrt(K1)*L) + 2*cos(sqrt(K1)*L)*sinh(sqrt(K1)*L) - 3*sinh(sqrt(K1)*L))/(5*K1**(3/2)) 
        T[2,1,3] = -2*K2*(sin(sqrt(K1)*L)*sinh(sqrt(K1)*L) + 2*cos(sqrt(K1)*L)*cosh(sqrt(K1)*L) - 2*cosh(sqrt(K1)*L))/(5*K1**2) 
        T[2,2,4] = -sqrt(K1)*L*sinh(sqrt(K1)*L)/2 
        T[2,3,4] = -L*cosh(sqrt(K1)*L)/2 + sinh(sqrt(K1)*L)/(2*sqrt(K1)) 
        T[3,0,2] = 2*K2*(3*sin(sqrt(K1)*L)*cosh(sqrt(K1)*L) + cos(sqrt(K1)*L)*sinh(sqrt(K1)*L) + sinh(sqrt(K1)*L))/(5*sqrt(K1)) 
        T[3,0,3] = 2*K2*(3*sin(sqrt(K1)*L)*sinh(sqrt(K1)*L) + cos(sqrt(K1)*L)*cosh(sqrt(K1)*L) - cosh(sqrt(K1)*L))/(5*K1) 
        T[3,1,2] = 2*K2*(sin(sqrt(K1)*L)*sinh(sqrt(K1)*L) - 3*cos(sqrt(K1)*L)*cosh(sqrt(K1)*L) + 3*cosh(sqrt(K1)*L))/(5*K1) 
        T[3,1,3] = 2*K2*(sin(sqrt(K1)*L)*cosh(sqrt(K1)*L) - 3*cos(sqrt(K1)*L)*sinh(sqrt(K1)*L) + 2*sinh(sqrt(K1)*L))/(5*K1**(3/2)) 
        T[3,2,4] = -sqrt(K1)*sinh(sqrt(K1)*L)/2 - K1*L*cosh(sqrt(K1)*L)/2 
        T[3,3,4] = -sqrt(K1)*L*sinh(sqrt(K1)*L)/2 
    if (K1 < 0):
        T[0,0,0] = K2*(cos(sqrt(K1)*L)**2 + cosh(L*sqrt(-K1)) - 2)/(3*K1) 
        T[0,0,1] = -2*K2*(1 - cosh(L*sqrt(-K1)))*sinh(L*sqrt(-K1))/(3*(-K1)**(3/2)) 
        T[0,0,4] = -L*sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[0,1,1] = K2*(-sinh(L*sqrt(-K1))**2 + 2*cosh(L*sqrt(-K1)) - 2)/(3*K1**2) 
        T[0,1,4] = -(-L*sqrt(-K1)*cosh(L*sqrt(-K1)) + sinh(L*sqrt(-K1)))/(2*sqrt(-K1)) 
        T[0,2,2] = K2*(cosh(sqrt(K1)*L)**2 - 3*cosh(L*sqrt(-K1)) + 2)/(5*K1) 
        T[0,2,3] = -2*K2*(-sin(2*L*sqrt(-K1))/2 + sinh(L*sqrt(-K1)))/(5*(-K1)**(3/2)) 
        T[0,3,3] = K2*(cosh(sqrt(K1)*L)**2 + 2*cosh(L*sqrt(-K1)) - 3)/(5*K1**2) 
        T[1,0,0] = -K2*(-sqrt(-K1)*sinh(L*sqrt(-K1)) - sqrt(-K1)*sinh(2*L*sqrt(-K1)))/(3*K1) 
        T[1,0,1] = 2*K2*(cosh(L*sqrt(-K1)) - cosh(2*L*sqrt(-K1)))/(3*K1) 
        T[1,0,4] = K1*L*cosh(L*sqrt(-K1))/2 - sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[1,1,1] = -2*K2*(cosh(L*sqrt(-K1)) - 1)*sinh(L*sqrt(-K1))/(3*(-K1)**(3/2)) 
        T[1,1,4] = L*sqrt(-K1)*sinh(L*sqrt(-K1))/2 
        T[1,2,2] = K2*(-sqrt(-K1)*sin(2*L*sqrt(-K1)) - 3*sqrt(-K1)*sinh(L*sqrt(-K1)))/(5*K1) 
        T[1,2,3] = -2*K2*(cos(2*L*sqrt(-K1)) - cosh(L*sqrt(-K1)))/(5*K1) 
        T[1,3,3] = 2*K2*(-sqrt(-K1)*sin(2*L*sqrt(-K1))/2 + sqrt(-K1)*sinh(L*sqrt(-K1)))/(5*K1**2) 
        T[2,0,2] = 2*K2*(-2*sin(L*sqrt(-K1))*sinh(L*sqrt(-K1)) - cos(L*sqrt(-K1))*cosh(L*sqrt(-K1)) + cos(L*sqrt(-K1)))/(5*K1) 
        T[2,0,3] = -2*K2*(sin(L*sqrt(-K1))*cosh(L*sqrt(-K1)) + sin(L*sqrt(-K1)) - 2*cos(L*sqrt(-K1))*sinh(L*sqrt(-K1)))/(5*(-K1)**(3/2)) 
        T[2,1,2] = -2*K2*(2*sin(L*sqrt(-K1))*cosh(L*sqrt(-K1)) - 3*sin(L*sqrt(-K1)) + cos(L*sqrt(-K1))*sinh(L*sqrt(-K1)))/(5*(-K1)**(3/2)) 
        T[2,1,3] = -2*K2*(-sin(L*sqrt(-K1))*sinh(L*sqrt(-K1)) + 2*cos(L*sqrt(-K1))*cosh(L*sqrt(-K1)) - 2*cos(L*sqrt(-K1)))/(5*K1**2) 
        T[2,2,4] = L*sqrt(-K1)*sin(L*sqrt(-K1))/2 
        T[2,3,4] = L*cos(L*sqrt(-K1))/2 - sin(L*sqrt(-K1))/(2*sqrt(-K1)) 
        T[3,0,2] = 2*K2*(sin(L*sqrt(-K1))*cosh(L*sqrt(-K1)) + sin(L*sqrt(-K1)) + 3*cos(L*sqrt(-K1))*sinh(L*sqrt(-K1)))/(5*sqrt(-K1)) 
        T[3,0,3] = 2*K2*(3*sin(L*sqrt(-K1))*sinh(L*sqrt(-K1)) - cos(L*sqrt(-K1))*cosh(L*sqrt(-K1)) + cos(L*sqrt(-K1)))/(5*K1) 
        T[3,1,2] = 2*K2*(sin(L*sqrt(-K1))*sinh(L*sqrt(-K1)) + 3*cos(L*sqrt(-K1))*cosh(L*sqrt(-K1)) - 3*cos(L*sqrt(-K1)))/(5*K1) 
        T[3,1,3] = -2*K2*(-3*sin(L*sqrt(-K1))*cosh(L*sqrt(-K1)) + 2*sin(L*sqrt(-K1)) + cos(L*sqrt(-K1))*sinh(L*sqrt(-K1)))/(5*(-K1)**(3/2)) 
        T[3,2,4] = -K1*L*cos(L*sqrt(-K1))/2 + sqrt(-K1)*sin(L*sqrt(-K1))/2 
        T[3,3,4] = -L*sqrt(-K1)*sin(L*sqrt(-K1))/2 
    return T 

