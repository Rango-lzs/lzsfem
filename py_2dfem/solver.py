
import numpy as np

def linear_solver(A,b):
    return np.linalg.inv(A) @ b