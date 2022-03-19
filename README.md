# README

To use these programs, ensure that all dependencies (`qiskit`, `bitstring`, `cirq`, etc.) are up to date.

The following sections contain instructions on how to use each program. Note that the programs themselves have a few adjustable constants. For example, the `USE_IBMQ` constant can be set to `True` or `False`, depending on whether circuits should be run on a simulator, or on a quantum computer. Additionally, the `SIMULATOR_NOISE` constant can be set to `True` or `False` - if `USE_IBMQ` is `False` and `SIMULATOR_NOISE` is `True`, then all simulations will be performed using noise models based on an actual IBM quantum computer's noise. Furthermore, the `VERBOSE` constant can be set to `True` or `False` to adjust some of the extra logging.

Please note that in some cases, problems with large `n` may be extremely slow to execute. In addition, note that running on actual IBM quantum computers requires that circuits use 5 or less qubits. Therefore, for some cases with circuits using too many qubits, it is not possible to run circuits on IBM quantum computers.

## Deutsch-Jozsa

```sh
python dj.py [n] [type: constant|balanced]
```

To run the program, `n` should be provided, as well as the function is constant or balanced (input exactly `constant` or `balanced`, with no quotes). Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, if verbose logging is enabled, the Result of the executed circuit, as provided by qiskit, will be printed. Finally, the result of the algorithm is printed. This will be a string stating that the selected function is "constant" or "balanced".

## Bernstein-Vazirani

```sh
python bv.py [n] [a] [b]
```

Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, if verbose logging is enabled, the Result of the executed circuit, as provided by qiskit, will be printed. Finally, the result of the algorithm is printed. This will be a string stating the calculated values of a and b, based on the results of the circuit execution.

## Simon

```sh
python simon.py [n] [s]
```
Upon running the program, several pieces of information will be output. First, the given parameters will be printed back to the user. Next, if verbose logging is enabled, for each new attempt at obtaining a constraint equation, the Result of the executed circuit, as provided by qiskit, will be printed. Once enough constraints are obtained, they will be output, as well. Finally, the result of the algorithm is printed. This will be a string stating the calculated values of s, as well as the actual value of s.

## Grover

```sh
python grover.py [NUM BITS] [NUM CASES]
```

Upon running the program, several pieces of information will be output. For each `n` between 1 and `NUM BITS` inclusive, information about the current case is output (including the number of bits, number of iteration, expected number of correct runs, actual number of correct nums, total run time, and average simulation run time). In addition, the Result information of the execution of the circuit is output, if verbose logging is enabled.

## QAOA

```sh
python qaoa.py [MIN SIZE] [MAX SIZE]
```

When the program runs, it will run the QAOA algorithm on the graphs described in the files contained in `adjacency_matrices`, specifically the graphs with MIN SIZE to MAX SIZE nodes. For each problem, the program prints out some information about the solution obtained, as well as some diagnostic information. In addition, the Result information of the execution of the circuit is output, if verbose logging is enabled.

## Shor

```sh
python shor.py [NUMBER TO FACTOR]
```

When the program runs, it will print out the steps it performs to transform this result into a factor of the number we wish to factor. When a factor is found, the program prints some information concerning the performance of the circuit before terminating successfully. In addition, the Result information of the execution of the circuit is output, if verbose logging is enabled.
