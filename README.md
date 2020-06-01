# Centrality Project

**Description:** CS 406 Project for Parallelizing Centrality Measurements: Degree 1&amp;2, Closeness and Betweenness Centralities

## Authors

- [Barış Batuhan Topal](https://github.com/barisbatuhan)

- [Çağhan Köksal](https://github.com/caghankoksal)

- [Furkan Ergün](https://github.com/furkaneergun)

- [Hakan Ogan Alpar](https://github.com/oalpar)

## Requirements

For GPU Code:

- CUDA 10.0
- gcc 5.3.0

For CPU Code:

- gcc 8.2.0

This implementation is tested on a HPC cluster running CentOS 6.5 and having 60 Intel(R) Xeon(R) E7-4870 v2 @ 2.30GHz as CPU and 12 GB Nvidia Tesla K40c as GPU with 2880 CUDA cores.

## Compilation & Run

- **CPU:** Includes a coarse-grained parallelization implementation using OpenMP. To run, either the `Makefile` inside the folder can be run or the command below can be executed: 

```
g++ ./ugraph.h main.cpp -O3 -fopenmp
./a.out <number_of_threads>
``` 

- **GPU:** Includes a hybrid parallelization approach based on both coarse grained and fine grained parallelization levels. To run, either the `Makefile` inside the folder can be run or the command below can be executed: 

```
nvcc -o cent_calc cent_gpu.cu
./cent_calc
```

**For both implementations**, the matrix file to run should be specified inside the main functions of the source codes. 

## Problem Statement

In today's world, graphs are being used widely in many different areas. For instance, while Google uses this data structure in “Google Maps” by representing the roads that connect different places as edges, Facebook analyzes their social network by representing each person as a vertex and their relations as edges [(Geeks For Geeks)](https://www.geeksforgeeks.org/applications-of-graph-data-structure/). Since graphs have a variety of use cases, measures and heuristics chosen to analyze these graph structures also vary. Programs and algorithms to analyze these graphs require lots of computation power. Moreover, because of the increasing data sizes, one may gain a significant time difference by parallelizing these heuristics.

We recognize that in the heart of most these network analysis are Breadth First Search (BFS), and if an end user decides to use different heuristics he or she would have to run multiple BFS algorithms to get the results for analysis. Our purpose is to unify multiple centrality measures under one program such that our program can give the results faster for cases where multiple centrality measures are queried and parallelize the overall ​BFS algorithm with the hope that it can give faster results for even a single centrality measure query. For this project, we specifically focus on the following centrality measures: degree 1&2 centrality, closeness centrality and betweenness centrality.

## Implementation

### Graph Representation

To narrow our problem space, in this project we only worked on undirected and unweighted sparse graphs. All these graphs are read and stored in [Compressed Sparse Row Representation (CSR)](https://www.researchgate.net/publication/324640550_A_survey_on_NoSQL_stores) to save space in sparse graphs. An example for CSR structure can be seen below [(Science Direct)](https://www.sciencedirect.com/topics/computer-science/graph-representation):

![CSR Example](images/csr_example.png)

### Centrality Measurements

We chose 4 different centrality measurements, which are using BFS structure during the calculation. By using the same BFS result for each of the solution, calculating all these 4 heuristics in a shorter time is aimed. 

#### Degree One Centrality

For every single node, the number of links held by this node is calculated.

#### Degree Two Centrality

For all vertices, the number of their neighbors and their neighbors’ neighbors is found. In other words, for every single vertex, Breadth-First Search (BFS) algorithm is run and the total count of nodes is computed, which has a distance closer than or equal to 2.

![Degree 2 Centrality Equation](images/deg2_formula.png)

#### Closeness Centrality

The closeness of each node to other nodes is found by running BFS and adding all the distances found as a result.

![Closeness Centrality Equation](images/cc_formula.png)

#### Betweenness Centrality

Measures the number of shortest paths between 2 different nodes (s and t), in which our vertex (v) lies on. The more our node is between these 2 other nodes, the higher will be the value for our node. The shortest paths between 2 nodes are determined by BFS algorithms.

![Betweenness Centrality Equation](images/bc_formula.png)

### CPU Level Parallelization

*TO BE FILLED*

### GPU Level Parallelization

*TO BE FILLED*

## Results & Discussion

To measure our performance, four graphs with different sizes are being used, which are retrieved from [Suite-Sparse Matrix Collection](http://faculty.cse.tamu.edu/davis/suitesparse.html):

|  | **494_bus.mtx** | **c-43.mtx** | **wing_nodal.mtx** | **wave.mtx** |
| :--- | :---: | :---: | :---: | :---: |
| **Number Of Nodes** | 494 | 11125 | 10937 | 156317 |
| **Number of Edges** | 1080 | 67400 | 75488 | 1059331 |

Since there is no work calculating all of these 4 centrality measurements at the same time, our results are compared with a [betweenness centrality implementation having CUDA GPU parallelization](https://github.com/pvgupta24/Graph-Betweenness-Centrality), since betweenness centrality is the most time consuming measurement among all these selected. The implementation we have chosen for comparison is directly implementing the [methods presented by NVIDIA](https://devblogs.nvidia.com/accelerating-graph-betweenness-centrality-cuda/). Our results are stated below (The ones stated with italic fonts are from the repository we have used for comparison):

| **Methods** | **494_bus.mtx** | **c-43.mtx** | **wing_nodal.mtx** | **wave.mtx** |
| :--- | ---: | ---: | ---: | ---: |
| ***Brandes without Optimization*** | 0.06 s | 57.550 s | 66.500 s | - |
| ***Brandes with -O3*** | 0.01 s | 11.370 s | 20.770 s | 4915.840 s |
| **CPU 1 Thread** | 0.013 s | 9.882 s | 25.568 s | 4578.670 s |
| **CPU 2 Threads** | 0.011 s | 5.240 s | 13.867 s | 2343.540 s |
| **CPU 4 Threads** | 0.012 s | 2.927 s | 7.812 s | 1198.080 s |
| **CPU 8 Threads** | 0.010 s | 1.640 s | 4.105 s | 631.664 s |
| **CPU 16 Threads** | 0.009 s | 0.926 s | 2.154 s | 332.555 s |
| **CPU 32 Threads** | 0.008 s | 0.525 s | 1.129 s | 164.977 s |
| ***Edge Parallel GPU*** | 0.112 s | 22.757 s | 99.278 s | - |
| ***Vertex Parallel GPU*** | 0.281 s | 222.499 s | 128.733 s | - |
| ***Fine Grained GPU*** | 0.354 s | 258.352 s | 79.210 s | - |
| ***Fine & Coarse Grained GPU*** | 0.061 s | 17.081 s | 5.448 s | 1482.232 s |
| **Our GPU Implementation** | 0.036 s | 25.807 s | 20.890 s | - |

*TO BE FILLED*

## Further Improvements

*TO BE FILLED*



