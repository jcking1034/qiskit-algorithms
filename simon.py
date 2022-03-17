"""Simons"""

import copy
import sys

import cirq
import numpy as np
import qiskit
import sympy
from bitstring import BitArray
from qiskit import IBMQ
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer.noise import NoiseModel
from qiskit.quantum_info.operators import Operator

# n is the size input of the function (number of non ancilla bits)

USE_IBMQ = False
SIMULATOR_NOISE = False
VERBOSE = True
NUM_SHOTS = 1

if USE_IBMQ or SIMULATOR_NOISE:
    IBMQ.save_account("YOUR_KEY_HERE")


def result_to_int(result, n):
    return int(max(result.keys(), key=lambda x: result[x]), base=2)


def result_to_bitlist(result, n):
    """ "qbit 0 (MSB) is first"""
    s = max(result.keys(), key=lambda x: result[x])

    return [int(i) for i in s[::-1]]


def bitstring_to_int(bitstring_list: list):
    assert type(bitstring_list) == list

    num = 0
    for i in range(len(bitstring_list)):
        num += bitstring_list[len(bitstring_list) - 1 - i] * 2**i
    return num


def get_U_f_matrix(n_inp, n_anc, f):
    side = 2 ** (n_inp + n_anc)
    ret = np.zeros((side, side))

    for i in range(side):
        ret[f(i >> n_anc) ^ i, i] = 1

    return ret


def simon(n: int, f) -> BitArray:
    constraints = []
    circuitRuntimes = []

    while len(constraints) < n - 1:
        print("Trying to get another constraint")
        circuit = create_simon_circuit(f, n)
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

        transpiled = qiskit.transpile(circuit, backend)
        job_sim = backend.run(transpiled, shots=NUM_SHOTS)

        result_sim = job_sim.result()
        if VERBOSE:
            print(result_sim)

        counts = result_sim.get_counts()
        # print(counts)

        bitstring_list = result_to_bitlist(counts, n)
        bitstring_int = result_to_int(counts, n)

        # end_time = time.perf_counter()
        # circuitRuntimes.append(end_time - start_time)

        # we don't care about 0 solns
        ret_matrix, _ = rref_f2(np.array(constraints + [bitstring_list]))
        if bitstring_int != 0 and np.any(ret_matrix[-1]):
            constraints.append(bitstring_list)

    # we have n linearly independent constraints and can solve the matrix equation
    print("constraints", constraints)
    _, sol = rref_f2(np.array(constraints))
    ret_s = bitstring_to_int(list(sol))
    if f(0) != f(ret_s):
        ret_s = 0

    return ret_s, sum(circuitRuntimes)


def rref_f2(matrix_orig: np.array):
    if len(matrix_orig) <= 0:
        raise RuntimeError("I am mad!")

    matrix = matrix_orig.copy()

    num_rows = len(matrix)
    num_cols = len(matrix[0])

    free_cols = np.ones(num_cols, dtype=int)

    for i in range(num_rows):

        leftmost1_row = None
        one_pos = len(matrix[0])

        for j in range(i, num_rows):
            np_nonz_tmp = np.nonzero(matrix[j])[0]
            if len(np_nonz_tmp) == 0:
                pass
            elif np_nonz_tmp[0] < one_pos:
                leftmost1_row = j
                one_pos = np_nonz_tmp[0]

        if leftmost1_row == None:
            break

        free_cols[one_pos] = 0

        # print(i, leftmost1_row)
        matrix[[i, leftmost1_row]] = matrix[[leftmost1_row, i]]

        for j in range(num_rows):
            if j == i:
                continue

            matrix[j] = (matrix[j] + matrix[i] * matrix[j][one_pos]) % 2

    # DO THE OTHER STUFF
    tmp1 = matrix.transpose().dot(matrix.dot(free_cols)) % 2
    return matrix, tmp1 | free_cols


def is_linearly_indp(list_of_bitstrings_orig: list, new_bitstring: list):
    list_of_bitstrings = list_of_bitstrings_orig[:]
    list_of_bitstrings.append(new_bitstring)
    M = sympy.Matrix(list_of_bitstrings)
    M, _ = M.rref()
    lastRow = list(M[-1, :])
    zeroVec = [0] * len(new_bitstring)
    return zeroVec != lastRow


def is_linearly_indp_matan(list_of_bitstrings_orig: list, new_bitstring: list):
    """
    do gausiann elimination
    find some vector with leftmost 1, use that to elminate all other vectors that have one in that poisiton
    dont use that row again, and look at columns further to the right
    """
    list_of_bitstrings = copy.deepcopy(list_of_bitstrings_orig)
    list_of_bitstrings.append(new_bitstring)
    n = len(list_of_bitstrings[0])
    counter = 0
    need_resorting = True

    while counter < len(list_of_bitstrings):
        # if we need to resort after doing some operations in the last step (or not re-sorting at all yet) do it!
        if need_resorting:
            list_of_bitstrings.sort(key=lambda x: bitstring_to_int(x), reverse=True)
            need_resorting = False

        if bitstring_to_int(list_of_bitstrings[len(list_of_bitstrings) - 1]) == 0:
            return False

        # make sure that there is only one bitstring that has this column poisiton on
        if list_of_bitstrings[counter][counter] == 1:
            need_resorting = True
            for i in range(len(list_of_bitstrings)):
                if i != counter:
                    list_of_bitstrings[i] = subtruct_bitstrings(
                        list_of_bitstrings[i], list_of_bitstrings[counter]
                    )
        counter += 1

    # do a last check after the last subtruction
    if [0] * n in list_of_bitstrings:
        return False
    else:
        return True


def subtruct_bitstrings(subtruct_from: list, subtructor: list):
    assert len(subtruct_from) == len(subtructor)
    return [(subtruct_from[i] - subtructor[i]) % 2 for i in range(len(subtructor))]


def run_simon_quantum_subroutine(circuit: cirq.Circuit, qbits, n) -> BitArray:
    """simulate the circuit and return the s value which was returned"""

    simulator = cirq.Simulator()
    result = simulator.run(circuit)
    return result_to_bitlist(result, n), result_to_int(result, n)


def create_simon_circuit(f, n: int):
    circuit = qiskit.QuantumCircuit(2 * n, n)

    for i in range(n):
        circuit.h(i)

    circuit.barrier()

    U_f = Operator(get_U_f_matrix(n, n, f))
    circuit.unitary(U_f, reversed(range(2 * n)), "U_f")

    circuit.barrier()

    for i in range(n):
        circuit.h(i)

    circuit.measure(range(n), range(n))

    return circuit


def constructSimonDict(n: int, s: int):
    rand_perm = np.random.permutation(2**n)
    simon_array = np.ones(2**n, dtype=int) - 2
    for i in range(2**n):
        if simon_array[i] == -1:
            simon_array[i] = simon_array[i ^ s] = rand_perm[i]
    return simon_array


if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        s = int(sys.argv[2])
        if n < 2:
            print("We need n >= 2")
            exit(1)
        if s < 0 or s > 2**n - 1:
            print("s should be in [0, 2**n - 1]")
            exit(1)
        print(f"{n} qubits, s={s}")
    except Exception:
        print("Usage: python simon.py [n] [s]")
        exit(1)

    simon_array = list(constructSimonDict(n, s))

    ret_s, average_circuit_runtime_per_circ = simon(
        n, lambda x: simon_array[x]
    )  # run simon routine

    print(f"\tReturned s is {ret_s}, actual s is {s}")

    # # start by constructing a function which meets the simon criterions
    # MAX_S = 10
    # MAX_N = 6

    # circ_avg_runtimes = []
    # tot_runtimes = []
    # for n in range(2, MAX_N):
    #     runtime_per_n = []
    #     circ_runtime_per_n = []
    #     for s_itr in range(0, MAX_S):
    #         # make simon function
    #         s = randint(0, 2**n - 1)

    #         simon_array = list(constructSimonDict(n, s))
    #         # build circuit make qbits

    #         start_time = time.perf_counter()  # start time

    #         ret_s, average_circuit_runtime_per_circ = simon(
    #             n, lambda x: simon_array[x]
    #         )  # run simon routine

    #         if ret_s != s:  # catch if simon failed
    #             print("error! ret_s:", ret_s, "s:", s, "simonDict", simon_array, n)
    #             exit(1)

    #         end_time = time.perf_counter()

    #         # add meassurement
    #         circ_runtime_per_n.append(average_circuit_runtime_per_circ)
    #         runtime_per_n.append(end_time - start_time)

    #     # add round messurements
    #     circ_avg_runtimes.append(circ_runtime_per_n)
    #     tot_runtimes.append(runtime_per_n)

    # print(
    #     "complete algorithm total runtimes (seconds)",
    #     [average(row) for row in tot_runtimes],
    #     "total_circ_run (seconds)",
    #     [average(row) for row in circ_avg_runtimes],
    # )
