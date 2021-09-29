from numba import njit
import numpy as np
from numpy import cos, sin, cosh, sinh, sqrt, tan


@njit
def compute_transport_fringe_in_Ex_matrix(element_parameters: list) -> np.ndarray:
    h: float = element_parameters[0] 
#     alpha: float = element_parameters[1] 
    beta1: float = element_parameters[3] 
    R1: float = element_parameters[7] 
    FINT: float = element_parameters[6] 
    Y_APPERTURE: float = element_parameters[5] 
    d: float = element_parameters[8] 
    h = h/(1+d)     
    phi1 = FINT*h*Y_APPERTURE*cos(beta1)**(-1)*(1+sin(beta1)**2)
    R = np.zeros((6,6))
    R[0,0] = 1 
    R[1,0] = h*tan(beta1) 
    R[1,1] = 1 
    R[2,2] = 1 
    R[3,2] = -h*tan(beta1 - phi1) 
    R[3,3] = 1 
    R[4,4] = 1
    R[5,5] = 1
    return R 
#             h,  # 0
#             self.K1.m_as('m**-2'),  # 1
#             self.K2.m_as('m**-3'),  # 2
#             face_angle,  # 3    
#             self.TILT.m_as('radian'),  # 4
#             HGAP, # 5
#             fringe_x,  # 6
#             self.R2 #7
def compute_transport_fringe_in_Ex_tensor(element_parameters: list) -> np.ndarray:
    h: float = element_parameters[0] 
#     alpha: float = element_parameters[1] 
    K1: float = element_parameters[1] 
    beta1: float = element_parameters[3] 
    R1: float = element_parameters[7] 
    FINT: float = element_parameters[6] 
    Y_APPERTURE: float = element_parameters[5] 
    d: float = element_parameters[8] 
    h = h/(1+d) 
    K1 = K1/(1+d) 
    T = np.zeros((6,6,6))
    phi1 = FINT*h*Y_APPERTURE*cos(beta1)**(-1)*(1+sin(beta1)**2) 
    T[0,0,0] = -h*tan(beta1)**2/2 
    T[0,2,2] = h*cos(beta1)**(-1)**2/2 
    T[1,0,0] = K1*tan(beta1) + (h*cos(beta1)**(-1))**3/(2*R1) 
    T[1,0,1] = h*tan(beta1)**2 
    T[1,0,4] = -h*tan(beta1) 
    T[1,2,2] = -K1*tan(beta1) + 0.5*h**2*tan(beta1) + h**2 + tan(beta1)**3 - h*cos(beta1)**(-1)**3/(2*R1) 
    T[1,2,3] = -h*tan(beta1)**2 
    T[2,0,2] = h*tan(beta1)**2 
    T[3,0,2] = -2*K1*tan(beta1) - h*cos(beta1)**(-1)**3/R1 
    T[3,0,3] = -h*tan(beta1)**2 
    T[3,1,2] = -h*cos(beta1)**(-1)**2 
    T[3,2,4] = -h*phi1*cos(beta1 - phi1)**(-2) + h*tan(beta1) 
    return T 



