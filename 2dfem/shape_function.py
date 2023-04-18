import numpy as np

'''
this module  calculate the shape function at giving coords at param space, 
actually at the gauss point as gauss quadrature numerical method is used
'''

def shape_quad4(r,s):
    #1*4
    N = 0.25*np.array(
    [(1-r)*(1-s),
    (1+r)*(1-s),
    (1+r)*(1+s),
    (1-r)*(1+s)
    ])

    #dNdr and dNds 2*4
    dNdr = 0.25*np.array([
        [-1+s,1-s,1+s,-1-s],
        [-1+r,-1-r,1+r,1-r]
    ])

    return N,dNdr

def elas_diff_2d(r,s,coord,shape_func):
    N,dNdr = shape_func(r,s)
    det, jaco_inv = jacoper(dNdr,coord)

    #[2 * 4]
    dNdx = jaco_inv @ dNdr

    H = np.zeros((2,8))
    B = np.zeros((3,8))

    # 2 * 8
    H[0,0::2] = N
    H[1,1::2] = N

    # 3 * 8
    B[0,0::2] = dNdx[0,:]
    B[1,1::2] = dNdx[1,:]
    B[2,0::2] = dNdx[1,:]
    B[2,1::2] = dNdx[0,:]

    return H, B, det


def jacoper(dNdr, coord):
    '''
    calculate the jacobe matrix , only depend on the coord
    for quad4 element : 2*4 @ 4 *2 = 2*2
    '''

    jaco = dNdr @ coord 
    det = np.linalg.det(jaco)

    jaco_rev = np.linalg.inv(jaco)

    return det,jaco_rev