import numpy as np

def read_input(file_path = ""):
    elems = np.loadtxt(file_path+"eles.txt", dtype= int, ndmin= 2)
    nodes = np.loadtxt(file_path+"nodes.txt",dtype= float,ndmin= 2)
    loads = np.loadtxt(file_path+"loads.txt",dtype= float, ndmin= 2)
    mats = np.loadtxt(file_path+"mater.txt", dtype=  float, ndmin= 2)
    return nodes, elems,loads,mats


#nodes, elems,loads,mats = read_input("lzsfem\square-4_elements\\")
#print(nodes)
#print(elems)
#print(loads)
#print(mats)