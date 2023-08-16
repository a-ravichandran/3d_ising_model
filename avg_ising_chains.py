import h5py
from matplotlib import pyplot as plt
import json
import numpy as np

dir_path = f"data/N10"
T = np.loadtxt(dir_path+"/T.txt")
with open(dir_path+"/MC_pars.json","r") as file:
    mc_pars = json.load(file)
del(file)

N_proc = mc_pars["N_proc"]
N_chain = mc_pars["N_chain"]
N_lat = 10 #mc_pars["N_lat"]

for i in range(len(T)):
    fig = plt.figure(figsize=[8,3])
    axs = fig.subplots(nrows=1, ncols=2)
    H_hist = np.zeros(shape=(N_chain,), dtype=np.float32)
    M_hist = np.zeros(shape=(N_chain,), dtype=np.float32)
    for idx in range(N_proc):
        with h5py.File(dir_path+f"/sim-{i:02d}/chain-{idx}.h5", "r") as f:
            H_hist += np.array(f["H_hist"][:])
            M_hist += np.array(abs(f["M_hist"][:]))

    with h5py.File(dir_path+f"/sim-{i:02d}/chain-avg.h5", "w") as fw:
        dataset = fw.create_dataset('H_hist', data=H_hist)
        dataset.attrs['descr'] = 'Hamiltonian History'

        dataset = fw.create_dataset('M_hist', data=M_hist)
        dataset.attrs['descr'] = 'Magnetization History'

    axs[0].plot(H_hist/(N_proc*(N_lat**3)), label=f"chain-{idx}", alpha=.9, ls="--")
    axs[1].plot(M_hist/(N_proc*(N_lat**3)), label=f"chain-{idx}", alpha=.9, ls="--")

    fig.suptitle(f"T={T[i]}")
    axs[0].set_title("Avg. Hamiltonian History")
    axs[1].set_title("Avg. |Magnetization| History")
    axs[0].set_xlabel("MC Steps")
    axs[1].set_xlabel("MC Steps")
    axs[0].grid()
    axs[1].grid()

    del(H_hist)
    del(M_hist)

    plt.tight_layout()
    fig.savefig(dir_path+f"/sim-{i:02d}/mc_history.pdf", facecolor='w', transparent=False)
    plt.show()

plt.close("all")
