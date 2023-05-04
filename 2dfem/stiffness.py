import numpy as np
import material as mat
import quadrature as gauss
import shape_function as shape


def stiff_func(eleType):
    ele_func = {
        1 : elas_tri3,
        2 : elas_quad4,
        3 : elas_tet4,
        4 : elas_hex8
    }

    try:
        return ele_func[eleType]
    except:
        raise ValueError("undefined element type!")


def spring(coord, params):
    pass

def truss(coord, params):
    pass

def elas_tri3(coord, params):
    pass

# 计算单元的刚度矩阵，跟单元类型，插值型函数，材料模型相关
def elas_quad4(coord, params):
    """Quadrilateral element with 4 nodes
    """
    stiff_mat = np.zeros([8, 8])
    mass_mat = np.zeros([8, 8])
    C = mat.elastic_2d(params[0],params[1])
    if len(params) == 2:
        dens = 1
    else:
        dens = params[-1]
    gpts, gwts = gauss.gauss_nd(2)
    for cont in range(gpts.shape[0]): 
        r, s = gpts[cont, :]
        H, B, det = shape.elas_diff_2d(r, s, coord, shape.shape_quad4)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += dens*factor* (H.T @ H)
    return stiff_mat, mass_mat


def elas_tet4(coord, params):
    pass

def elas_hex8(coord, params):
    pass