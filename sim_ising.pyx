# cython: language_level=3
import numpy as np
from comp_quan import calc_H, calc_M

def check_flip(lat,NN,i,j,k,B,J,beta):
    # lat: current lattice state
    # NN: nearest neighbor information
    # i,j,k: position of spin in focus
    # B: magnetic field
    # J: interaction energy parameter
    # beta: inverse temperature
    # returns boolean for whether spin should be flip and delta_H

    cdef float del_H = 0
    del_H += 2*B*lat[i,j,k]
    
    for pos in NN[i,j,k]:
        del_H += 2*J*lat[pos]*lat[i,j,k]

    cdef float prob 
    cdef float temp

    if del_H<=0:
        return 1, del_H
    else:
        prob = np.exp(-beta*del_H)
        temp = np.random.rand()
        
        if temp<prob:
            return 1, del_H
        else:
            del_H = 0
            return 0, del_H
        
def sim_ising(lat, NN, beta, B, J, N):
    # N: number of flip steps
    # lat: initial lattice state
    # NN: nearest neighbours information
    # beta: inverse temperature
    # B: magnetic field * magnetic moment
    # J: exchange/interaction energy
    # N: number of time steps
    # returns history of Energy, Net Magnetization and Lattice States

    H_hist = np.zeros(shape=(N,), dtype=np.float32)
    del_H_hist = np.zeros(shape=(N,), dtype=np.float32)
    M_hist = np.zeros(shape=(N,), dtype=np.int16)
    lat_hist = np.zeros(shape=(N,lat.shape[0],lat.shape[1],lat.shape[2]), dtype=np.int8)

    lat_hist[0] = np.copy(lat)
    H_hist[0] = calc_H(lat=lat_hist[0], NN=NN, B=B, J=J) # initial H
    M_hist[0] = calc_M(lat=lat_hist[0], NN=NN, B=B, J=J) # initial M

    cdef float del_H
    cdef int flip

    for a in range(1,N):
        lat_hist[a] = np.copy(lat_hist[a-1])
        i = np.random.randint(low=0, high=lat.shape[0])
        j = np.random.randint(low=0, high=lat.shape[1])
        k = np.random.randint(low=0, high=lat.shape[2])
        # print(i,j)

        flip, del_H = check_flip(lat=lat_hist[a],NN=NN,i=i,j=j,k=k,B=B,J=J,beta=beta)

        if flip == 1:
            H_hist[a] = H_hist[a-1] + del_H
            M_hist[a] = M_hist[a-1] - 2*lat_hist[a][i,j,k]
            lat_hist[a,i,j,k] *= -1

        else:
            H_hist[a] = np.copy(H_hist[a-1])
            M_hist[a] = np.copy(M_hist[a-1])

        del_H_hist[a] = del_H
    return H_hist, M_hist, lat_hist, del_H_hist

def sim_ising_no_lat(lat, NN, beta, B, J, N):
    # N: number of flip steps
    # lat: initial lattice state
    # NN: nearest neighbours information
    # beta: inverse temperature
    # B: magnetic field * magnetic moment
    # J: exchange/interaction energy
    # N: number of time steps
    # returns history of Energy, Net Magnetization and Lattice States

    H_hist = np.zeros(shape=(N,), dtype=np.float32)
    del_H_hist = np.zeros(shape=(N,), dtype=np.float32)
    M_hist = np.zeros(shape=(N,), dtype=np.int16)

    H_hist[0] = calc_H(lat=lat, NN=NN, B=B, J=J) # initial H
    M_hist[0] = calc_M(lat=lat, NN=NN, B=B, J=J) # initial M
    lat_new = np.copy(lat)
    cdef float del_H
    cdef int flip

    for a in range(1,N):
        i = np.random.randint(low=0, high=lat.shape[0])
        j = np.random.randint(low=0, high=lat.shape[1])
        k = np.random.randint(low=0, high=lat.shape[2])
        # print(i,j)

        flip, del_H = check_flip(lat=lat_new,NN=NN,i=i,j=j,k=k,B=B,J=J,beta=beta)

        if flip == 1:
            H_hist[a] = H_hist[a-1] + del_H
            M_hist[a] = M_hist[a-1] - 2*lat_new[i,j,k]
            lat_new[i,j,k] *= -1

        else:
            H_hist[a] = np.copy(H_hist[a-1])
            M_hist[a] = np.copy(M_hist[a-1])

        del_H_hist[a] = del_H
    return H_hist, M_hist, lat_new, del_H_hist
