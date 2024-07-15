import readfile as io
import assemble as asm
import solver
import postprocessor as pos

'''
the main routine of fem analysis steps
'''

def analysis():
    elems ,nodes,loads,mats = io.read_input("")


# 读取模型信息
elems ,nodes,loads,mats = io.read_input("square-4_elements\\")

# 计算组装矩阵
dof_map, bc, eqn= asm.dof_mapping(elems, nodes)

# 组装刚度举证和载荷
stiff_matrix = asm.assemble(elems, nodes)
load_vec = asm.loadasem(loads, bc, eqn)

# 求解
disp = solver.linear_solver(stiff_matrix, load_vec)

disp_complete = pos.complete_disp(bc, nodes, disp)

pos.plot_node_field(disp_complete, nodes, elems, plt_type="contourf", levels=12,
                    savefigs=True, title=["Ux","Uy"], figtitle=["Ux","Uy"],
                    filename=["Ux","Uy"])

print(disp)

# 后处理