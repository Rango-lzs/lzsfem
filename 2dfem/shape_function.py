
import numpy as np

def jacoper(dNdr, coord):
    '''
    calculate the jacobe matrix , only depend on the coord
    for quad4 element : 2*4 @ 4 *2 = 2*2
    '''

    jaco = dNdr @ coord 
    det = np.linalg.det(jaco)

    jaco_rev = np.linalg.inv(jaco)

    return det,jaco_rev