 
from comp_quan import *
import h5py
import json
import numpy as np

dir_path = f"data/N10"
T = np.loadtxt(dir_path+"/T.txt")
with open(dir_path+"/MC_pars.json","r") as file:
    mc_pars = json.load(file)
del(file)

N_proc = mc_pars["N_proc"]
N_chain = mc_pars["N_chain"]
N_lat = mc_pars["N_lat"]

n_corr_H = 1000
n_corr_M = 1000
n_therm_H = int(5e4)
n_therm_M = int(5e4)

quant = {"chi": [], "C": [], "H": [], "M": [], "T":T.tolist()}

for i in range(len(T)):
    with h5py.File(dir_path+f"/sim-{i:02d}/chain-avg.h5", "r") as f:
        quant["C"].append(calc_avg_C(H_hist=f["H_hist"][:], n_therm=n_therm_H, n_corr=n_corr_H, size=N_lat, beta=1.0/T[i]))
        quant["chi"].append(calc_avg_chi(M_hist=f["M_hist"][:], n_therm=n_therm_M, n_corr=n_corr_M, size=N_lat, beta=1.0/T[i]))
        quant["H"].append(calc_avg_H(H_hist=f["H_hist"][:], n_therm=n_therm_H, n_corr=n_corr_H, size=N_lat))
        quant["M"].append(abs(calc_avg_M(M_hist=f["M_hist"][:], n_therm=n_therm_M, n_corr=n_corr_M, size=N_lat)))
        
with open(dir_path+"/quant.json","w") as file:
     json.dump(quant,file)
del(file)
