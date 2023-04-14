import readfile as io
import assemble as asm

'''
the main routine of fem analysis steps
'''

def analysis():
    elems ,nodes,loads,bcs = io.read_input("")



elems ,nodes,loads,bcs = io.read_input("square-4_elements\\")

dof_map = asm.dof_mapping(elems, nodes)
