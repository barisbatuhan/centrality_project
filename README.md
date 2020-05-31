# Centrality Project

**Description:** CS 406 Project for Parallelizing Centrality Measurements: Degree 1&amp;2, Closeness and Betweenness Centralities

## Authors

- [Barış Batuhan Topal](https://github.com/barisbatuhan)

- [Çağhan Köksal](https://github.com/caghankoksal)

- [Furkan Ergün](https://github.com/furkaneergun)

- [Hakan Ogan Alpar](https://github.com/oalpar)

## Requirements

- CUDA 10.0
- gcc 5.3.0

## Folders and Descriptions

- **CPU:** Includes a coarse-grained parallelization implementation. Number of threads are set as stated in the `Makefile`. 
- **GPU:** Includes a hybrid parallelization approach based on both coarse grained and fine grained parallelization levels. 

For both implementations, the matrix file to run should be specified inside the main functions of the source codes. 

## Introduction

*TO BE FILLED*

## Implementation

### First Steps

*TO BE FILLED*

### CPU Level Parallelization

*TO BE FILLED*

### GPU Level Parallelization

*TO BE FILLED*

## Results & Discussion

*TO BE FILLED*

| Methods | 494_bus.mtx | c-43.mtx | shuttle_eddy.mtx | wing_nodal.mtx | wave.mtx |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Brandes's Algorithm | 0 | 0 | 0 | 0 | 0 |
| CPU 1 Thread | 0 | 0 | 0 | 0 | 0 |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| CPU 2 Threads | 0 | 0 | 0 | 0 | 0 |
| CPU 4 Threads | 0 | 0 | 0 | 0 | 0 |
| CPU 8 Threads | 0 | 0 | 0 | 0 | 0 |
| CPU 16 Threads | 0 | 0 | 0 | 0 | 0 |
| CPU 32 Threads | 0 | 0 | 0 | 0 | 0 |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Fine Grained GPU | 0 | 0 | 0 | 0 | 0 |
| Edge Parallel GPU | 0 | 0 | 0 | 0 | 0 |
| Vertex Parallel GPU | 0 | 0 | 0 | 0 | 0 |
| Fine & Coarse Grained GPU | 0 | 0 | 0 | 0 | 0 |
| Our GPU Implementation | 0 | 0 | 0 | 0 | 0 |

## Further Improvements

*TO BE FILLED*

## References

*TO BE FILLED*


