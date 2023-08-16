# cython: language_level=3
import numpy as np

def calc_H(lat, NN, B, J):
    # lat: current lattice state
    # NN: nearest neighbors of lattice points
    # B: magnetic field
    # J: interaction energy parameter
    # returns the Hamiltonian of the system
    cdef float H = 0
    for i in range(lat.shape[0]):
        for j in range(lat.shape[1]):
            for k in range(lat.shape[2]):
                H += -B*lat[i,j,k]
                temp = 0
                for pos in NN[(i,j,k)]:
                    temp += -J*lat[i,j,k]*lat[pos]
                H += temp*0.5
    return H

def calc_M(lat, NN, B, J):
    # lat: current lattice state
    # NN: nearest neighbors of lattice points
    # B: magnetic field
    # J: interaction energy parameter
    # returns the Magnetization of the system
    cdef float M = 0
    M = np.sum(lat)

    return M

def calc_avg_H(H_hist, n_therm, n_corr, size):
    cdef float H_avg = 0
    cdef int temp = 0
    for i in range(n_therm, H_hist.shape[0], n_corr):
        H_avg += H_hist[i]
        temp+=1
    H_avg /= temp

    return H_avg/(size**3)

def calc_avg_M(M_hist, n_therm, n_corr, size):
    cdef float M_avg = 0
    cdef int temp=0
    for i in range(n_therm, M_hist.shape[0], n_corr):
        M_avg += M_hist[i]
        temp+=1
    M_avg /= temp

    return M_avg/(size**3)

def calc_avg_chi(M_hist, n_therm, n_corr, beta, size):
    cdef float M_avg = 0
    cdef float M2_avg = 0
    cdef int temp = 0
    for i in range(n_therm, M_hist.shape[0], n_corr):
        M_avg += M_hist[i]
        M2_avg += (M_hist[i]**2)
        temp+=1
    M_avg /= temp
    M2_avg /= temp

    return ((M2_avg-(M_avg**2))*beta*2)/(size**2)

def calc_avg_C(H_hist, n_therm, n_corr, beta, size):
    cdef float H_avg = 0
    cdef float H2_avg = 0
    cdef int temp=0
    for i in range(n_therm, H_hist.shape[0], n_corr):
        H_avg += H_hist[i]
        H2_avg += (H_hist[i]**2)
        temp+=1
    H_avg /= temp
    H2_avg /= temp

    return ((H2_avg-(H_avg**2))*beta**2)/(size**2)