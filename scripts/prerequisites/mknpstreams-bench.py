from itertools import repeat, count
import npstreams as ns
import numpy as np
from tqdm import tqdm
from timeit import default_timer as timer
from pathlib import Path

OUTPUT = Path("data") / "introduction" / "npstreams-benchmark"
OUTPUT.mkdir(exist_ok=True)


def benchmark_arrsize(size):
    """
    Run an average of 10 arrays of size `size`. Monitor the running time.

    Parameters
    ----------
    size : int
        Array size.

    Returns
    -------
    nptime : float
        NumPy computation time in seconds.
    nstime : float
        npstreams computation time in seconds.
    """
    stream = [np.random.random((size, size)) for _ in range(10)]

    start = timer()
    nssum = ns.average(stream)
    nstime = timer() - start

    start = timer()
    npsum = np.average(np.stack(list(stream), axis=-1), axis=2)
    nptime = timer() - start

    return nptime, nstime


def benchmark_seqsize(size):
    """
    Run an average of ``size`` arrays of size (512, 512). Monitor the running time.

    Parameters
    ----------
    size : int
        Array sequence size.

    Returns
    -------
    nptime : float
        NumPy computation time in seconds.
    nstime : float
        npstreams computation time in seconds.
    """
    start = timer()
    nssum = ns.average(np.random.random((256, 256)) for _ in range(size))
    nstime = timer() - start

    stack = np.empty(shape=(256, 256, size), dtype=float)
    start = timer()
    np.stack([np.random.random((256, 256)) for _ in range(size)], axis=-1, out=stack)
    npsum = np.average(stack, axis=2)
    nptime = timer() - start

    return nptime, nstime


def memory_profile(size):
    """
    Run an average of 10 arrays of size `size`.

    Parameters
    ----------
    size : int
        Array size.

    Returns
    -------
    npmem : float
        Maximum memory usage during NumPy computation in MB.
    nsmem : float
        Maximum memory usage during npstreams computation in MB.
    """
    stream = [np.random.random((size, size)) for _ in range(10)]
    stack = np.stack(stream, axis=-1)
    return stack.nbytes / 1e6, stream[0].nbytes / 1e6


sizes = list(range(512, 4096, 256)) + list(range(4096, 8192, 512))
seqsizes = range(64, 2000, 50)

np_times = list()
np_seq_times = list()
np_seq_par_times = list()
np_memories = list()

ns_times = list()
ns_seq_times = list()
ns_seq_par_times = list()
ns_memories = list()

for size in tqdm(sizes):
    nptime, nstime = benchmark_arrsize(size)
    np_times.append(nptime)
    ns_times.append(nstime)

    npmem, nsmem = memory_profile(size)
    np_memories.append(npmem)
    ns_memories.append(nsmem)

for seqsize in tqdm(seqsizes):
    nptime, nstime = benchmark_seqsize(seqsize)
    np_seq_times.append(nptime)
    ns_seq_times.append(nstime)

np_times = np.array(np_times)
ns_times = np.array(ns_times)

np_seq_times = np.array(np_seq_times)
ns_seq_times = np.array(ns_seq_times)

np_memories = np.array(np_memories)
ns_memories = np.array(ns_memories)

times = np.empty(shape=(len(sizes), 3), dtype=float)
times[:, 0] = sizes
times[:, 1] = np_times
times[:, 2] = ns_times
np.save(OUTPUT / "times.npy", times)

memory = np.empty(shape=(len(sizes), 3), dtype=float)
memory[:, 0] = sizes
memory[:, 1] = np_memories
memory[:, 2] = ns_memories
np.save(OUTPUT / "memory.npy", memory)

seqtimes = np.empty(shape=(len(seqsizes), 3), dtype=float)
seqtimes[:, 0] = seqsizes
seqtimes[:, 1] = np_seq_times
seqtimes[:, 2] = ns_seq_times
np.save(OUTPUT / "seqtimes.npy", seqtimes)
