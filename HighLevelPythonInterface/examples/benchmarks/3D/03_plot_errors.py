import numpy as np
import matplotlib.pyplot as plt

err_idx = 3
labels = []
data = []
for diff in [1,3,5]:
    raw = np.loadtxt("results/errors_and_times_dtmm_diff=%d.dat" % (diff,))
    data.append({"dy":raw[:,1],"t":raw[:,2],"err":raw[:,err_idx]})
    labels.append("DTMM(%d)" % (diff,))

raw = np.loadtxt("results/errors_and_times_bpm.dat")
data.append({"dy":raw[:,1],"t":raw[:,2],"err":raw[:,err_idx]})
labels.append("BPM")

Nd = len(labels)
x = np.arange(Nd)
w = 0.35

i = 1

plt.subplot(2,1,1)
plt.bar(x, [d["err"][i] for d in data], w)
plt.gca().set_xticks(x)
plt.gca().set_xticklabels(labels)
plt.gca().set_yscale('log')

plt.subplot(2,1,2)
plt.bar(x, [d["t"][i] for d in data], w)
plt.gca().set_xticks(x)
plt.gca().set_xticklabels(labels)

plt.show()
