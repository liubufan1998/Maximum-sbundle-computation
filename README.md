# Efficient Maximum s-bundle search via local vertex connectivity

This repository implements the maximum *s*-bundle computation algorithm **SymBD** proposed in our paper. If you are using the code, please cite our paper.


## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "SymBD", which corresponds to our SymBD algorithm.

## Get datasets
The three datasets used in our paper can be found from the links below.

1. 2nd DIMACS graphs collection: http://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/clique/

2. 10th DIMACS graphs collection: https://networkrepository.com/dimacs10.php

3. Real-world graphs collection: https://networkrepository.com/index.php


## Data format
Our program only accepts a **binary input data file** (i.e., with a **.bin** file extension). If the data is not in binary form, you need to use our provided data transformation tool `translate.cpp` which is contained in the directory `data_transform_tool` to transform the data. Additionally, we have included an **original data file** `bio-celegans.mtx` within the same file directory, which can assist in understanding the native format of the graph data to be read.

We provide an example of converting a '.mtx' data file into binary format below. 

    g++ -O3 -std=c++11 translate.cpp -o a.out -w
    ./a.out bio-celegans.mtx bio-celegans.bin

## Run the code
The usage procedure for our program is as follows.
```sh
$ ./SymBD -g {path_to_graph} -s {s_value}
```

An example of computing the exact maximum 6-bundle for our provided example graph 'bio-pdb1HYS.bin' is as follows
```sh
$ ./SymBD -g datasets/bio-pdb1HYS.bin -s 6
```
## Result Analysis
If you are using our provided example above, then you will see the following outputs (may be there are some slight differences):
```sh
**** SymBK_H (Release) build at 16:11:27 Jul 11 2024 ***
File: bio-pdb1HYS n= 36417 m= 2154174 s = 6
Graph init ok
Now we find a larger solution, and its size is 59
*** Degeneracy s-bundle size: 58, max_core: 74, UB: 80, Time: 32,257 (microseconds)
*** After core shrink: n = 33,929, m = 2,048,332 (undirected)
*** After core_truss_copruning: n = 3279, m = 127259 (undirected)
*** 1-hop subgraph degen s-bundle size: 60, Time: 3,119 (microseconds)
*** After core_truss_copruning: n = 3279, m = 127259 (undirected)
2hop-Degen find a larger solution of size 61
2hop-Degen find a larger solution of size 62
2hop-Degen find a larger solution of size 63
*** 2-hop-Degen sbundle size: 63, Time: 13,944 (microseconds)
*** After core_truss_copruning: n = 345, m = 11289 (undirected)
*** The time for preprocessing is 1,286,236 (microseconds)
search_cnt: 10, ave_density: 0.97234, min_density: 0.93862
Branch_nodes: 518
*** Search time: 75,990
	Maximum 6-bundle Size: 63, Total Time: 1,362,248 (microseconds)
```
From the obtained empirical results, we can extract the following key information:

- **Graph Data Characteristics:**
  - **File Name**: `bio-pdb1HYS`
  - **Vertices (n)**: 36,417
  - **Edges (m)**: 2,154,174

- **s-Bundle Sizes Obtained via Our Three-stage Heuristic Approach (The details can be found in Section 5.2 of our paper.):**
  - **Stage 1: Degeneracy s-Bundle Size**: 58
  - **Stage 2: 1-hop Subgraph Degen s-Bundle Size**: 60
  - **Stage 3: 2-hop-Degen s-Bundle Size**: 63

- **Performance Metrics:**
  - **Preprocessing Time**: 1,286,236 microseconds
  - **Branches Generated**: 518
  - **Maximum 6-bundle Size**: 63
  - **Total Runtime**: 1,362,248 microseconds
