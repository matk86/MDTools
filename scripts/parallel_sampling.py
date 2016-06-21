import numpy as np
# uncomment the following on matgen
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import multiprocessing
from multiprocessing import Pool, Manager


moves = {1: [1, 0, 0], 
         2: [0, 1, 0], 
         3: [0, 0, 1], 
         4: [-1, 0, 0], 
         5: [0, -1, 0], 
         6: [0, 0, -1]}

def sampler(nmoves):
    prev_move = 1
    total_move = np.array(moves[prev_move])
    for i in range(nmoves):
        prev_move = np.random.randint(1, 7)
        total_move += np.array(moves[prev_move])
    return np.linalg.norm(total_move)

def normal_distribution(x, mu, sigma):
    return np.exp(-(x-mu)**2/(2*sigma**2)) / np.sqrt(2*np.pi) / sigma

def main(args):
    (samples, nmoves, nsamples) = args
    for i in range(nsamples):
        samples.append(sampler(nmoves))


if __name__=='__main__':
    nsamples = 100000
    nmoves = 100
    manager = Manager()
    samples = manager.list()
    # number of process to run simultaneously.
    nprocs = multiprocessing.cpu_count()
    p = Pool(nprocs)
    p.map(main, [(samples, nmoves, nsamples/nprocs)]*nprocs)
    avg = np.average(samples)
    std = np.std(samples)
    print avg, std  
    n, bins, _ = plt.hist(np.array(samples), 50)
    normal = [normal_distribution(x, avg, std) for x in bins]
    plt.plot(bins, np.array(normal)*max(n)*10, 'r')
    plt.savefig("distribution.png")
    plt.show()
