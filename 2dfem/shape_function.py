import numpy as np

'''
this module  calculate the shape function at giving coords at param space, 
actually at the gauss point as gauss quadrature numerical method is used
'''

#查找形函数要点
#1.插值节点的位置
#2.插值函数的形式，多维的拉格朗日插值

#三节点三角形
def shape_tri3(r,s):
    """
    s
    3
    |
    |\
    | \
    |__\ __ r
    1   2
    """
    N = np.array([1-r-s, r ,s])
    dNdr = np.array([
        [-1,1,0]
        [-1,0,1]
    ]}
    return N, dNdr

#3节点 二阶三角形单元
def shape_tri6(r,s):
    N = np.array(
        [(1 - r - s) - 2*r*(1 - r - s) - 2*s*(1 - r - s),
         r - 2*r*(1 - r - s) - 2*r*s,
         s - 2*r*s - 2*s*(1-r-s),
         4*r*(1 - r - s),
         4*r*s,
         4*s*(1 - r - s)])
    dNdr = np.array([
        [4*r + 4*s - 3, 4*r - 1, 0, -8*r - 4*s + 4, 4*s, -4*s],
        [4*r + 4*s - 3, 0, 4*s - 1, -4*r, 4*r, -4*r - 8*s + 4]])
    return N, dNdr
    

#四节点四边形等参元
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

## 8节点二阶四边形
def shape_quad8(r, s):
    """
    Shape functions and derivatives for a 8-noded serendipity element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s).
    """
    N = 0.25*np.array(
        [(1.0 - r)*(1.0 - s) - (1.0 - r)*(1.0 - s**2) - (1.0 - s)*(1.0 - r**2),
         (1.0 + r)*(1.0 - s) - (1.0 + r)*(1.0 - s**2) - (1.0 - s)*(1.0 - r**2),
         (1.0 + r)*(1.0 + s) - (1.0 + r)*(1.0 - s**2) - (1.0 + s)*(1.0 - r**2),
         (1.0 - r)*(1.0 + s) - (1.0 - r)*(1.0 - s**2) - (1.0 + s)*(1.0 - r**2),
         2.0*(1.0 - s)*(1.0 - r**2),
         2.0*(1.0 + r)*(1.0 - s**2),
         2.0*(1.0 + s)*(1.0 - r**2),
         2.0*(1.0 - r)*(1.0 - s**2)])
    dNdr = 0.25*np.array([
        [-2.0*r*(s - 1.0) - s**2 + s,
         -2.0*r*(s - 1.0) + s**2 - s,
         -2.0*r*(-s - 1.0) + s**2 + s,
         -2.0*r*(-s - 1.0) - s**2 - s,
         -2.0*r*(2.0 - 2.0*s),
         2.0 - 2.0*s**2,
         -2.0*r*(2.0*s + 2.0),
         2.0*s**2 - 2.0],
        [-r**2 + r - 2.0*s*(r - 1.0),
         -r**2 - r - 2.0*s*(-r - 1.0),
         r**2 + r - 2.0*s*(-r - 1.0),
         r**2 - r - 2.0*s*(r - 1.0),
         2.0*r**2 - 2.0,
         -2.0*s*(2.0*r + 2.0),
         2.0 - 2.0*r**2,
         -2.0*s*(2.0 - 2.0*r)]])
    return N, dNdr


## 3D elements  四面体单元C3D4
def shape_tet4(r, s, t):
    """
    Shape functions and derivatives for a linear tetrahedron

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    N = np.array([1 - r - s - t, r, s, t])
    dNdr = np.array([
        [-1, 1, 0, 0],
        [-1, 0, 1, 0],
        [-1, 0, 0, 1]])
    return N, dNdr

## 3D elements  四面体单元C3D10
def shape_tet10(r, s, t):
    """
    Shape functions and derivatives for a linear tetrahedron

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    pass


#8节点六面体 C3D8
def shape_hex8(r, s, t):
    """
    Shape functions and derivatives for a trilinear element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    N = np.array([
        (1 - r)*(1 - s)*(1 - t), (1 - s)*(1 - t)*(r + 1),
        (1 - t)*(r + 1)*(s + 1), (1 - r)*(1 - t)*(s + 1),
        (1 - r)*(1 - s)*(t + 1), (1 - s)*(r + 1)*(t + 1),
        (r + 1)*(s + 1)*(t + 1), (1 - r)*(s + 1)*(t + 1)])
    dNdr = np.array([
        [(1 - t)*(s - 1), (1 - s)*(1 - t),
         (1 - t)*(s + 1), (1 - t)*(-s - 1),
         (1 - s)*(-t - 1), (1 - s)*(t + 1),
         (s + 1)*(t + 1), -(s + 1)*(t + 1)],
        [(1 - t)*(r - 1), (1 - t)*(-r - 1),
         (1 - t)*(r + 1), (1 - r)*(1 - t),
         -(1 - r)*(t + 1), -(r + 1)*(t + 1),
         (r + 1)*(t + 1), (1 - r)*(t + 1)],
        [-(1 - r)*(1 - s), -(1 - s)*(r + 1),
         -(r + 1)*(s + 1), -(1 - r)*(s + 1),
         (1 - r)*(1 - s), (1 - s)*(r + 1),
         (r + 1)*(s + 1), (1 - r)*(s + 1)]])
    return 0.125*N, 0.125*dNdr


#8节点六面体 C3D20
def shape_hex20(r, s, t):
    """
    Shape functions and derivatives for a trilinear element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    N = np.array([
        (1 - r)*(1 - s)*(1 - t), (1 - s)*(1 - t)*(r + 1),
        (1 - t)*(r + 1)*(s + 1), (1 - r)*(1 - t)*(s + 1),
        (1 - r)*(1 - s)*(t + 1), (1 - s)*(r + 1)*(t + 1),
        (r + 1)*(s + 1)*(t + 1), (1 - r)*(s + 1)*(t + 1)])
    dNdr = np.array([
        [(1 - t)*(s - 1), (1 - s)*(1 - t),
         (1 - t)*(s + 1), (1 - t)*(-s - 1),
         (1 - s)*(-t - 1), (1 - s)*(t + 1),
         (s + 1)*(t + 1), -(s + 1)*(t + 1)],
        [(1 - t)*(r - 1), (1 - t)*(-r - 1),
         (1 - t)*(r + 1), (1 - r)*(1 - t),
         -(1 - r)*(t + 1), -(r + 1)*(t + 1),
         (r + 1)*(t + 1), (1 - r)*(t + 1)],
        [-(1 - r)*(1 - s), -(1 - s)*(r + 1),
         -(r + 1)*(s + 1), -(1 - r)*(s + 1),
         (1 - r)*(1 - s), (1 - s)*(r + 1),
         (r + 1)*(s + 1), (1 - r)*(s + 1)]])
    return 0.125*N, 0.125*dNdr


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