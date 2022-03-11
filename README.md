# README

To use these programs, ensure that all dependencies (`qiskit`, `bitstring`, etc.) are up to date.

The following sections contain instructions on how to use each program. Note that the programs themselves have a few adjustable constants. For example, the `USE_IBMQ` constant can be set to `True` or `False`, depending on whether circuits should be run on a simulator, or on a quantum computer.

## Deutsch-Jozsa

```sh
python dj.py [num qubits] [type: 'constant'|'balanced']
```

Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, the Result of the executed circuit, as provided by qiskit, will be printed. Finally, the result of the algorithm is printed. This will be a string stating that the selected function is "constant" or "balanced".

## Bernstein-Vazirani

```sh
python bv.py [num qubits] [a] [b]
```

Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, the Result of the executed circuit, as provided by qiskit, will be printed. Finally, the result of the algorithm is printed. This will be a string stating the calculated values of a and b, based on the results of the circuit execution.

## Simon

```sh
python simon.py [num qubits] [s]
```
Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, for each new attempt at obtaining a constraint equation, the Result of the executed circuit, as provided by qiskit, will be printed. Once enough constraints are obtained, they will be output, as well. Finally, the result of the algorithm is printed. This will be a string stating the calculated values of s, as well as the actual value of s.

## Grover

```sh
python grover.py [NUM BITS] [NUM CASES]
```
Upon running the program, several pieces of information will be output. For each `n` between 1 and `NUM BITS` inclusive, information about the current case is output (including the number of bits, number of iteration, expected number of correct runs, actual number of correct nums, total run time, and average simulation run time). In addition, the run information of the circuit simulation is output.

## QAOA

## Shor
