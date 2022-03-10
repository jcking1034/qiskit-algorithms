"Bernstein-Vazirani"
import itertools
import sys
import time

import numpy as np
import qiskit
from qiskit import IBMQ
from qiskit.providers.aer import QasmSimulator
from qiskit.quantum_info.operators import Operator

USE_IBMQ = True
NUM_SHOTS = 1024


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


def get_U_f_matrix(n, f):
    side = 2 ** (n + 1)
    ret = np.zeros((side, side))

    for xb in range(2 ** (n + 1)):
        x = xb >> 1
        b = xb & 1
        y = f(x)
        ret[xb, x << 1 | (b ^ y)] = 1

    # print(ret)
    return ret


def create_bv_circuit(n, f):
    circuit = qiskit.QuantumCircuit(n + 1, n)

    # Initialize last qubit to 1
    circuit.x(n)

    circuit.barrier()

    for i in range(n + 1):
        circuit.h(i)

    circuit.barrier()

    U_f = Operator(get_U_f_matrix(n, f))
    circuit.unitary(U_f, reversed(range(n + 1)), "U_f")

    circuit.barrier()

    for i in range(n):
        circuit.h(i)

    circuit.measure(range(n), range(n))

    return circuit


def run(n, f):
    b = f(0)

    circuit = create_bv_circuit(n, f)
    # print("Circuit")
    # print(circuit.draw())

    if USE_IBMQ:
        IBMQ.save_account("YOUR_KEY_HERE")
        provider = IBMQ.load_account()
        backend = provider.backend.ibmq_quito
    else:
        backend = QasmSimulator()

    # Execute the circuit on the qasm simulator.
    # We've set the number of repeats of the circuit
    # to be 1024, which is the default.
    transpiled = qiskit.transpile(circuit, backend)
    job_sim = backend.run(transpiled, shots=NUM_SHOTS)

    # Grab the results from the job.
    result_sim = job_sim.result()

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
        print("Usage: python bv.py [num qubits] [a] [b]")
        exit(1)

    actual_a, actual_b = run(n, lambda x: add(multiply(a, x), b))
    print(f"\tActual: a={actual_a}, b={actual_b}")

    # assert (actual_a, actual_b) == (a, b)
