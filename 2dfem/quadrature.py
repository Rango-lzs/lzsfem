
import numpy as np
from itertools import product

'''
this module  implement the gauss quadrature 
give the gauss points and weights for different kind if parametric element
'''

#%% General
def gauss_1d(npts):
    """Return Gauss points and weights for Gauss quadrature in 1D

    Parameters
    ----------
    npts : int
      Number of quadrature points.

    Returns
    -------
    wts : ndarray
      Weights for the Gauss-Legendre quadrature.
    pts : ndarray
      Points for the Gauss-Legendre quadrature.
    """
    if npts == 2:
        pts = [-0.577350269189625764, 0.577350269189625764]
        wts = [1.00000000000000000, 1.00000000000000000]
    elif npts == 3:
        pts = [-0.774596669241483377, 0, 0.774596669241483377]
        wts = [0.555555555555555556, 0.888888888888888889,
               0.555555555555555556]
    elif npts == 4:
        pts = [-0.861136311594052575, -0.339981043584856265,
               0.339981043584856265, 0.861136311594052575]
        wts = [0.347854845137453857, 0.652145154862546143,
               0.652145154862546143, 0.347854845137453857]  
    else:
        msg = "The number of points should be in [2, 10]"
        raise ValueError(msg)

    return pts, wts


def gauss_nd(npts, ndim=2):
    """
    Return Gauss points and weights for Gauss quadrature in
    an ND hypercube.

    Parameters
    ----------
    npts : int
      Number of quadrature points.

    Returns
    -------
    nd_wts : ndarray
      Weights for the Gauss-Legendre quadrature.
    nd_pts : ndarray
      Points for the Gauss-Legendre quadrature.
    """
    pts, wts = gauss_1d(npts)
    nd_pts = np.array(list(product(pts, repeat=ndim)))
    nd_wts = product(wts, repeat=ndim)
    nd_wts = [np.prod(nd_wt) for nd_wt in nd_wts]
    return nd_pts, nd_wts