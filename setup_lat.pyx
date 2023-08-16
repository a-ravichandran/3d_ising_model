# cython: language_level=3
import numpy as np

def Lattice_3D(n, start):
    # n: dimensions of 3D Lattice n X n X n
    # returns n X n x n lattice array and the nearest neighbour dictionary

    if start == "cold":
        lat = np.ones(shape=(n,n,n), dtype=np.int8)
    
    elif start == "hot":
        lat = np.random.choice([-1, 1], size=(n, n, n), p=[.5,.5])
    
    cdef dict NN = {} # nearest neighbours

    for i in range(n):
        for j in range(n):
            for k in range(n):
                # ~ [up, right, down, left, top, bottom]     
                NN[(i,j,k)] = [((i-1)%n, j, k),
                            (i,(j+1)%n, k),
                            ((i+1)%n, j, k),
                            (i,(j-1)%n, k),
                            (i,j,(k+1)%n),
                            (i,j,(k-1)%n)]

    return lat, NN