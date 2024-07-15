import numpy as np

def elastic_2d(E,nu):
    '''
    constitutive matric of plane stress
    '''
    
    C = np.zeros((3, 3))
    enu = E/(1 - nu**2)
    mnu = (1 - nu)/2
    C[0, 0] = enu
    C[0, 1] = nu*enu
    C[1, 0] = C[0, 1]
    C[1, 1] = enu
    C[2, 2] = enu*mnu
    return C