import numpy as np

'''
this module  calculate the shape function at giving coords at param space, 
actually at the gauss point as gauss quadrature numerical method is used
'''

def shape_quad4(r,s):
    #1*4
    N = np.array(
    [(1-r)*(1-s),
    (1+r)*(1-s),
    (1+r)*(1+s),
    (1-r)(1+s)
    ])

    #dNdr and dNds 2*4
    dNdr = np.array([
        [-1+s,1-s,1+s,-1-s],
        [-1+r,-1-r,1+r,1-r]
    ])

    return N,dNdr
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