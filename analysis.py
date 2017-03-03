from sys import argv
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

filename_beg = 'n-100k_p-'

bcast_times = []; peval_times = []
bcast_t_all = np.array([])
peval_t_all = np.array([])

p = range(2, 13);

for i in p:
    filename = filename_beg + i + '.txt'
    with open(filename) as f:
        for line in f:
            x, result, bcast_t, peval_t = line.split()
            bcast_times.append(bcast_t)
            peval_times.append(peval_times)

    bcast_t_all = append(bcast_t_all, np.mean(np.array(bcast_times)))
    peval_t_all = append(peval_t_all, np.mean(np.array(peval_times)))

f1 = plt.figure()
f2 = plt.figure()

f1.plot(p, bcast_t_all)
f1.savefig('bcast-t_vs_p.pdf')

f2.plot(p, peval_t_all)
f2.savefig('peval-t_vs_p.pdf')

