"Bernstein-Vazirani"
import sys

import numpy as np
import qiskit
from qiskit import IBMQ
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer.noise import NoiseModel
from qiskit.quantum_info.operators import Operator

USE_IBMQ = True
SIMULATOR_NOISE = False
VERBOSE = True
NUM_SHOTS = 1024

if USE_IBMQ or SIMULATOR_NOISE:
    IBMQ.save_account("YOUR_KEY_HERE")


def multiply(x, y):
    """Bitwise multiplication of two ints"""
    ret = 0
    while x or y:
        ret += (x & 1) * (y & 1)
        x >>= 1
        y >>= 1
    return bool(ret % 2)


def add(x, y):
    """The sum of two bits, mod 2."""
    return (x + y) % 2


def get_U_f_matrix(n_inp, n_anc, f):
    side = 2 ** (n_inp + n_anc)
    ret = np.zeros((side, side))

    for i in range(side):
        ret[bool(f(i >> n_anc)) ^ i, i] = 1

    return ret


def create_bv_circuit(n, f):
    circuit = qiskit.QuantumCircuit(n + 1, n)

    # Initialize last qubit to 1
    circuit.x(n)

    circuit.barrier()

    for i in range(n + 1):
        circuit.h(i)

    circuit.barrier()

    U_f = Operator(get_U_f_matrix(n, 1, f))
    circuit.unitary(U_f, reversed(range(n + 1)), "U_f")

    circuit.barrier()

    for i in range(n):
        circuit.h(i)

    for i in range(n):
        circuit.measure(i, n - i - 1)

    return circuit


def run(n, f):
    b = f(0)

    circuit = create_bv_circuit(n, f)
    # print("Circuit")
    # print(circuit.draw())

    if USE_IBMQ:
        provider = IBMQ.load_account()
        backend = provider.backend.ibmq_quito
    elif SIMULATOR_NOISE:
        provider = IBMQ.load_account()
        noise_backend = provider.backend.ibmq_quito
        noise_model = NoiseModel.from_backend(noise_backend)
        backend = QasmSimulator(noise_model=noise_model)
    else:
        backend = QasmSimulator()

    # Execute the circuit on the qasm simulator.
    # We've set the number of repeats of the circuit
    # to be 1024, which is the default.
    transpiled = qiskit.transpile(circuit, backend)
    job_sim = backend.run(transpiled, shots=NUM_SHOTS)

    # Grab the results from the job.
    result_sim = job_sim.result()

    if VERBOSE:
        print(result_sim)

    counts = result_sim.get_counts()
    # print(counts)

    a = int(max(counts.keys(), key=lambda x: counts[x]), base=2)

    return a, b

    # breakpoint()
    # assert len(counts) == 1
    # for a in counts.keys():
    #     return int(a[::-1], base=2), b


if __name__ == "__main__":
    # n = 5
    # a = 12
    # b = 1

    try:
        n = int(sys.argv[1])
        a = int(sys.argv[2])
        b = int(bool(int(sys.argv[3])))
        print(f"{n} qubits, a={a}, b={b}")
    except Exception:
        print("Usage: python bv.py [n] [a] [b]")
        exit(1)

    actual_a, actual_b = run(n, lambda x: add(multiply(a, x), b))
    print(f"\tActual: a={actual_a}, b={actual_b}")

    # assert (actual_a, actual_b) == (a, b)
