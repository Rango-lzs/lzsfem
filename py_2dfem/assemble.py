import numpy as np
import stiffness as stiff

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
            else:
                bc[node_idx][dof] = -1
    return bc, equ+1  


def dof_mapping(elems, nodes):
    '''
    记录单元自由度对应的全局自由度

    返回值： dof_map[elem][ndof] = bc[elem_nodes: -1][2] 
    '''  
    elem_num  = len(elems)
    ndof=2
    dof_map = np.zeros((elem_num, ndof*4),dtype = int)

    bc, equ_num = eq_counter(nodes)

    for elem_idx in range(elem_num):
        elem_nodes = elems[elem_idx][3:]

        dof_map[elem_idx] = bc[elem_nodes].flatten()

    return dof_map, bc ,equ_num

def assemble(elems, nodes):
    '''
    组装单刚成总刚stiff_matrix
    '''
    #dof_map[i][0:8]
    # 总刚的大小
   
    dof_map, bc, eqn = dof_mapping(elems, nodes)
    stiff_matrix = np.zeros((eqn,eqn))

    elems_num = len(elems)

    for ele_idx in range(elems_num):
        # 获取单元的单元类型，材料参数，然后计算
        # coord 2*4
        coord = nodes[elems[ele_idx, 3:], 1:3]
        param = [1,0.3]
        eleType = 2 #elems[ele_idx, 4]
        kloc,mloc = stiff.stiff_func(eleType)(coord, param)

        ele_dofs = kloc.shape[0]
        ele_map = dof_map[ele_idx]

        for i in range(ele_dofs):
            for j in range(ele_dofs):
                global_i = ele_map[i]
                global_j = ele_map[j]
                if global_i != -1 and global_j != -1:
                    stiff_matrix[global_i,global_j] = stiff_matrix[global_i,global_j]+ kloc[i,j]

    return stiff_matrix


def loadasem(loads, bc_array, neq, ndof_node=2):
    '''
    组装载荷矩阵
    '''
    # 遍历载荷数组
    load_count = loads.shape[0]
    rhs_vec = np.zeros(neq)

    for idx in range(load_count):
        node = int(loads[idx,0])
        for dof in range(ndof_node):
            global_dof = int(bc_array[node,dof])
            if global_dof != -1:
                rhs_vec[global_dof] = loads[idx,dof+1]
    return rhs_vec



        
