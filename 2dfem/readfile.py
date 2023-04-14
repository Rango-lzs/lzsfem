import numpy as np

'''file format
element : id, eletype ，matType, node...
1 1 1 2 3 4
2 1 2 3 4 5

node: id, x,y , dof bc (0 for free, -1 for restrainted) 
1 0.5 0.5 -1 0

loads:
node_id, fx, fy

mats：
E1,nu1
E2,nu2
...

'''

def read_input(file_path = ""):
    elems = np.loadtxt(file_path+"eles.txt", dtype= int, ndmin= 2)
    nodes = np.loadtxt(file_path+"nodes.txt",dtype= float,ndmin= 2)
    loads = np.loadtxt(file_path+"loads.txt",dtype= float, ndmin= 2)
    mats = np.loadtxt(file_path+"mater.txt", dtype=  float, ndmin= 2)
    return elems, nodes,loads,mats


#nodes, elems,loads,mats = read_input("lzsfem\square-4_elements\\")
#print(nodes)
#print(elems)
#print(loads)
#print(mats)