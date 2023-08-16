import matplotlib.pyplot as plt
import json
import numpy as np

fig = plt.figure(figsize=[6,6])
axs = fig.subplots(nrows=2, ncols=2)

with open(dir_path+"/quant.json","r") as file:
    quant = json.load(file)
del(file)

start = 0
end = None

ax = axs[0,0]
ax.set_title("H")
ax.plot(quant["T"][start:end], quant["H"][start:end], marker="x", ls="--", color="blue")
ax.grid()

ax = axs[1,0]
ax.set_title("C")
ax.plot(quant["T"][start:end], quant["C"][start:end], marker="x", ls="--", color="blue")
ax.grid()

ax = axs[0,1]
ax.set_title("M")
ax.plot(quant["T"][start:end], quant["M"][start:end], marker="x", ls="--", color="red")
ax.grid()

ax = axs[1,1]
ax.set_title("chi")
ax.plot(quant["T"][start:end], quant["chi"][start:end], marker="x", ls="--", color="red")
ax.grid()

plt.show()
fig.savefig(dir_path+f"/quant.pdf", facecolor='w', transparent=False)
plt.close("all")
