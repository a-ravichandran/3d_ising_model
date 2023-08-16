import h5py
import multiprocessing as mp
from configparser import ConfigParser
import os
import numpy as np
import time
import json

from setup_lat import Lattice_3D
from sim_ising import sim_ising_no_lat


# MULTIPROCESSING WORKER FUNC ---------------------------
def worker(lat,NN,beta,B,J,N_chain,i,seed,dir_path, save_lat_hist):
    np.random.seed(seed=seed)
    H_hist, M_hist, _ , _ = sim_ising_no_lat(lat=lat, NN=NN, beta=beta, B=B, J=J, N=N_chain)

    with h5py.File(f'{dir_path}/chain-{i}.h5', 'a') as f:

        dataset = f.create_dataset('H_hist', data=H_hist)
        dataset.attrs['descr'] = 'Hamiltonian History'
        del(H_hist)

        dataset = f.create_dataset('M_hist', data=M_hist)
        dataset.attrs['descr'] = 'Magnetization History'
        del(M_hist)

        if save_lat_hist:
            dataset = f.create_dataset('lat_hist', data=lat_hist)
            dataset.attrs['descr'] = 'Lattice History'
            del(lat_hist)


def master_func(conf_i):
    # READING CONFIGURATION FILE ----------------------------
    config = ConfigParser()
    config.optionxform = str
    config.read(f"configs/control_ising-{conf_i}.ini")

    control = {}
    for key in config["ising"]:
        control[key] = config["ising"][key]

    N_lat = int(control["N_lat"])
    N_chain = int(control["N_chain"])
    N_proc = int(control["N_proc"])
    B = float(control["B"])
    J = float(control["J"])
    kB = float(control["kB"])
    save_lat_hist = True if int(control["save_lat_hist"])==1 else False
    T_start = int(control["T_start"])
    T_end = int(control["T_end"])
    T = np.linspace(float(control["T_start"]), float(control["T_end"]),int(control["N_T"]))
    start = control["start"]

    # DIRECTORIES TO STORE DATA --------------------------
    dir_path = f"data/N{N_lat}T{T_start}_{T_end}"
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    par_dict = {}
    par_dict["N_chain"] = N_chain
    par_dict["N_lat"] = N_lat
    par_dict["N_proc"] = N_proc
    par_dict["T"] = T.tolist()
    with open(dir_path+"/MC_pars.json", "w") as fp:
        json.dump(par_dict, fp)

    np.savetxt(dir_path+"/T.txt", T)

    for i in range(len(T)):
        path = dir_path + f"/sim-{i:02d}"
        if not os.path.exists(path):
            os.mkdir(path)

    for j in range(len(T)):
        # SETTING UP LATTICE AND DATA --------------------------
        lat, NN = Lattice_3D(n=N_lat, start=start)
        path_T = dir_path+f"/sim-{j:02d}"
        with h5py.File(f"{path_T}/sim_pars.h5", 'a') as f:
            # g = f.create_group("pars")
            f.attrs["decr"] = "Simulation Parameters"

            d = f.create_dataset("kB", data=kB)
            d = f.create_dataset("B", data=B)
            d = f.create_dataset("J", data=J)
            d = f.create_dataset("N_lat", data=N_lat)
            d = f.create_dataset("N_chain", data=N_chain)
            d = f.create_dataset("N_proc", data=N_proc)
            d = f.create_dataset("start", data=start)
            d = f.create_dataset("T", data=T[j])

        t = time.time()

        # SIMULATING ISING -------------------------------------
        beta = 1.0/T[j]
        p_list = []
        seeds = np.random.randint(low=0, high=256, size=N_proc, dtype=int)
        print(f"\nT = {T[j]}:")
        for i in range(N_proc):
            p = mp.Process(target=worker, args=[lat, NN, beta, B, J, N_chain, i, seeds[i], path_T, save_lat_hist])
            print(f"Starting {p.name}")
            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

        for p in p_list:
            p.close()
            print(f"Finished {p.name} in ", time.time()-t, f" seconds)")
        del(p_list)
        #  -------------------------------------

pm_list = []
for i in range(1,5):
    pm = mp.Process(target=master_func, args=[i])
    print(f"Starting Master Process {i}")
    pm.start()
    pm_list.append(pm)

for pm in pm_list:
    pm.join()

for pm_i in range(len(pm_list)):
    pm_list[pm_i].close()
    print(f"Finished Master Process {pm_i}")
del(pm_list)
