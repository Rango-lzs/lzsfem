import numpy as np

'''
此模块用于单元刚度矩阵组装，应包含如何功能
1、节点自由度计数
2、局部自由度和全局自由度的映射关系
3、总体刚度矩阵的组装，稀疏矩阵的存储
4、总体载荷矩阵的组装

'''

def eq_counter(nodes):
    '''
    给节点自由度进行编号，记录节点自由度对应的全局自由度
    bc[node_size][2]
    '''
    size = len(nodes)
    node_dofs = 2
    bc = np.zeros((size,2))
    equ = -1
    for node_idx in range(size):
        dof_value = nodes[node_idx][-2:]
        for dof in range(node_dofs):
            value = dof_value[dof]
            if value == 0:
                equ += 1
                bc[node_idx][dof] = equ   
    return bc, equ  


def dof_mapping(elems, nodes):
    '''
    记录单元自由度对应的全局自由度

    返回值： dof_map[elem][ndof] = bc[elem_nodes: -1][2] 
    '''  
    elem_num  = len(elems)
    ndof=2
    dof_map = np.zeros((elem_num, ndof*4))

    bc, equ_num = eq_counter(nodes)

    for elem_idx in range(elem_num):
        elem_nodes = elems[elem_idx][3:]

        dof_map[elem_idx] = bc[elem_nodes].flatten()

    return dof_map
